//! Parallel clade accumulation, partial-taxa mapping, and O(N log N) Euler-RMQ CA heights.

use std::path::PathBuf;
use std::sync::Arc;

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::helper::progress::{self, ProgressReporter};
use crate::parser::tree::{
    AnnotationTarget, NodeId, ParseOptions, Rootedness, TaxonTable, Tree, TreeFmt,
    read_trees_chunked_par,
};
use crate::report::stats::TreeSummaryFileStats;

use super::bitset::{
    CladeBitsArena, CladeId, ZobristTable, bit_is_set, canonicalize_for_mask,
    canonicalize_unrooted, first_set_bit, popcount, set_bit,
};
use super::consensus::{TreeLocation, hash_clade_list, map_tree_nodes_to_clades};
use super::stats_acc::{AnnotationKeyTable, CladeAnnotationsAcc, MetricAccumulator};
use super::{BranchLengthStrategy, BurninSpec, RootedMode, SummaryMethod, TreeSummaryOptions};

const CHUNK_SIZE: usize = 512;

/// Aggregated statistics for one unique clade or bipartition across input trees.
#[derive(Debug, Clone)]
pub(crate) struct CladeStats {
    pub(crate) count: u64,
    pub(crate) decisive_count: u64,
    pub(crate) branch_lengths: MetricAccumulator,
    pub(crate) heights: MetricAccumulator,
    pub(crate) ca_heights: MetricAccumulator,
    pub(crate) annotations: CladeAnnotationsAcc,
}

impl CladeStats {
    fn new(store_samples: bool) -> Self {
        Self {
            count: 0,
            decisive_count: 0,
            branch_lengths: MetricAccumulator::new(store_samples),
            heights: MetricAccumulator::new(store_samples),
            ca_heights: MetricAccumulator::new(store_samples),
            annotations: CladeAnnotationsAcc::default(),
        }
    }

    fn merge(
        &mut self,
        other: Self,
        other_keys: &AnnotationKeyTable,
        global_keys: &mut AnnotationKeyTable,
    ) {
        self.count += other.count;
        self.decisive_count += other.decisive_count;
        self.branch_lengths.merge(other.branch_lengths);
        self.heights.merge(other.heights);
        self.ca_heights.merge(other.ca_heights);
        self.annotations
            .merge(other.annotations, other_keys, global_keys);
    }
}

/// Global interned clade/bipartition table plus tip and root summaries.
#[derive(Debug, Clone)]
pub(crate) struct GlobalCladeTable {
    pub(crate) taxa: Arc<TaxonTable>,
    pub(crate) zobrist: ZobristTable,
    pub(crate) arena: CladeBitsArena,
    pub(crate) by_hash: AHashMap<(u64, u64), Vec<CladeId>>,
    pub(crate) clades: Vec<CladeStats>,
    pub(crate) tip_stats: Vec<CladeStats>,
    pub(crate) root_stats: CladeStats,
    pub(crate) keys: AnnotationKeyTable,
    pub(crate) store_samples: bool,
}

impl GlobalCladeTable {
    pub(crate) fn new(taxa: Arc<TaxonTable>, store_samples: bool) -> Self {
        let n_taxa = taxa.len();
        let zobrist = ZobristTable::new(n_taxa);
        let arena = CladeBitsArena::new(zobrist.words());
        let tip_stats = (0..n_taxa)
            .map(|_| CladeStats::new(store_samples))
            .collect();
        Self {
            taxa,
            zobrist,
            arena,
            by_hash: AHashMap::new(),
            clades: Vec::new(),
            tip_stats,
            root_stats: CladeStats::new(store_samples),
            keys: AnnotationKeyTable::default(),
            store_samples,
        }
    }

    pub(crate) fn lookup(&self, bits: &[u64], hash: (u64, u64)) -> Option<CladeId> {
        let bucket = self.by_hash.get(&hash)?;
        bucket
            .iter()
            .find(|&&cid| self.arena.get(cid) == bits)
            .copied()
    }

    pub(crate) fn insert_zero_clade(&mut self, bits: &[u64], hash: (u64, u64)) -> CladeId {
        if let Some(cid) = self.lookup(bits, hash) {
            return cid;
        }
        let cid = self.arena.push(bits);
        self.clades.push(CladeStats::new(self.store_samples));
        self.by_hash.entry(hash).or_default().push(cid);
        cid
    }

    fn lookup_or_insert(&mut self, bits: &[u64], hash: (u64, u64)) -> CladeId {
        if let Some(cid) = self.lookup(bits, hash) {
            return cid;
        }
        self.insert_zero_clade(bits, hash)
    }
}

/// Unique tree topology entry tracked for MAP selection and credible-set CSV output.
#[derive(Debug, Clone)]
pub(crate) struct TopologyRecord {
    pub(crate) hash: (u64, u64),
    pub(crate) clade_ids: Vec<CladeId>,
    pub(crate) count: u64,
    pub(crate) first_location: TreeLocation,
}

/// Per-tree clade membership list kept when scoring MCC trees.
#[derive(Debug, Clone)]
pub(crate) struct TreeCladeRecord {
    pub(crate) location: TreeLocation,
    pub(crate) clade_ids: Vec<CladeId>,
}

/// Complete Pass 1 output ready for consensus/MCC/MAP/Target selection and Pass 2 CA heights.
pub(crate) struct SummaryRunData {
    pub(crate) table: GlobalCladeTable,
    pub(crate) topologies: Vec<TopologyRecord>,
    pub(crate) tree_records: Vec<TreeCladeRecord>,
    pub(crate) total_read: usize,
    pub(crate) total_burned_in: usize,
    pub(crate) total_summarized: usize,
    pub(crate) rooted: bool,
    pub(crate) union_ntax: usize,
    pub(crate) is_partial_target: bool,
}

/// Per-file resolved burn-in plan so Pass 1 and Pass 2 skip identical tree prefixes.
#[derive(Debug, Clone)]
pub(crate) struct FilePlan {
    pub(crate) path: PathBuf,
    pub(crate) total_trees: usize,
    pub(crate) burnin_skip: usize,
}

pub(crate) fn compute_burnin_skip(total_trees: usize, spec: BurninSpec) -> Result<usize> {
    let skip = match spec {
        BurninSpec::Count(n) => {
            if n >= total_trees {
                bail!(
                    "Burn-in count ({}) leaves 0 trees out of {} tree(s)",
                    n,
                    total_trees
                );
            }
            n
        }
        BurninSpec::Fraction(frac) => {
            if !(0.0..1.0).contains(&frac) {
                bail!("Burn-in fraction must be in [0.0, 1.0), got {}", frac);
            }
            let n = ((total_trees as f64) * frac).floor() as usize;
            if n >= total_trees {
                bail!(
                    "Burn-in fraction ({:.1}%) leaves 0 trees out of {}",
                    frac * 100.0,
                    total_trees
                );
            }
            n
        }
    };
    Ok(skip)
}

pub(crate) fn needs_sample_storage(options: &TreeSummaryOptions) -> bool {
    options.hpd_level > 0.0
        || options.branch_lengths == BranchLengthStrategy::Median
        || options.summarize_annotations
}

/// Run Pass 1 directly on an in-memory slice of trees.
pub(crate) fn summarize_in_memory_pass1(
    trees: &[Tree],
    target: Option<&Tree>,
    options: &TreeSummaryOptions,
) -> Result<SummaryRunData> {
    let total_read = trees.len();
    let skip = compute_burnin_skip(total_read, options.burnin)?;
    let retained = &trees[skip..];
    if retained.is_empty() {
        bail!("No trees remain after burn-in");
    }

    let first_tree = target.unwrap_or(&retained[0]);
    let rooted = resolve_rootedness(first_tree, options.rooted_mode)?;
    let store_samples = needs_sample_storage(options);

    if options.allow_partial_taxa {
        if options.method == SummaryMethod::Target {
            let tgt = target.context("Target method requires a target tree")?;
            return run_partial_target_in_memory(retained, tgt, rooted, options, total_read, skip);
        }
        let (shared_taxa, union_ntax) = discover_shared_taxa_in_memory(retained, target)?;
        let pruned: Vec<Tree> = retained
            .par_iter()
            .map(|t| prune_tree_to_taxa(t, &shared_taxa))
            .collect::<Result<Vec<_>>>()?;
        let mut res = run_identical_taxa_in_memory(
            &pruned,
            shared_taxa,
            rooted,
            options,
            store_samples,
            total_read,
            skip,
        )?;
        res.union_ntax = union_ntax;
        return Ok(res);
    }

    let taxa = build_taxon_table_from_tree(first_tree)?;
    let union_ntax = taxa.len();
    run_identical_taxa_in_memory(
        retained,
        taxa,
        rooted,
        options,
        store_samples,
        total_read,
        skip,
    )
    .map(|mut r| {
        r.union_ntax = union_ntax;
        r
    })
}

/// Run Pass 1 across tree files in bounded parallel chunks.
pub(crate) fn summarize_files_pass1(
    plans: &[FilePlan],
    input_fmt: TreeFmt,
    target: Option<&Tree>,
    options: &TreeSummaryOptions,
    progress: Option<&dyn ProgressReporter>,
) -> Result<(SummaryRunData, Vec<TreeSummaryFileStats>)> {
    let parse_opts = ParseOptions {
        parse_annotations: options.summarize_annotations,
        ..Default::default()
    };

    let first_tree = if let Some(tgt) = target {
        tgt.clone()
    } else {
        let p = &plans[0];
        crate::parser::tree::read_single_tree_at(&p.path, input_fmt, &parse_opts, p.burnin_skip)?
    };
    let rooted = resolve_rootedness(&first_tree, options.rooted_mode)?;
    let store_samples = needs_sample_storage(options);

    if options.allow_partial_taxa {
        if options.method == SummaryMethod::Target {
            let tgt = target.context("Target method requires a target tree")?;
            return run_partial_target_files(
                plans,
                input_fmt,
                &parse_opts,
                tgt,
                rooted,
                options,
                progress,
            );
        }
        let (shared_taxa, union_ntax) =
            discover_shared_taxa_files(plans, input_fmt, &parse_opts, target, progress)?;
        return run_identical_taxa_files(
            plans,
            input_fmt,
            &parse_opts,
            shared_taxa,
            union_ntax,
            true,
            rooted,
            options,
            store_samples,
            progress,
        );
    }

    let taxa = build_taxon_table_from_tree(&first_tree)?;
    let union_ntax = taxa.len();
    run_identical_taxa_files(
        plans,
        input_fmt,
        &parse_opts,
        taxa,
        union_ntax,
        false,
        rooted,
        options,
        store_samples,
        progress,
    )
}

fn run_identical_taxa_in_memory(
    retained: &[Tree],
    taxa: Arc<TaxonTable>,
    rooted: bool,
    options: &TreeSummaryOptions,
    store_samples: bool,
    total_read: usize,
    total_burned_in: usize,
) -> Result<SummaryRunData> {
    let taxon_index = build_taxon_index(&taxa);
    let n_taxa = taxa.len();
    let zobrist = ZobristTable::new(n_taxa);
    let need_clade_list = options.method == SummaryMethod::Mcc;
    let need_heights = rooted;

    let per_chunk: Vec<ChunkAccum> = retained
        .par_chunks(CHUNK_SIZE)
        .enumerate()
        .map(|(chunk_idx, chunk)| {
            let base = chunk_idx * CHUNK_SIZE;
            process_tree_slice(
                chunk,
                |local_i| TreeLocation::InMemory(base + local_i + total_burned_in),
                &taxa,
                &taxon_index,
                &zobrist,
                rooted,
                store_samples,
                options.summarize_annotations,
                need_heights,
                need_clade_list,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let mut global = GlobalCladeTable::new(taxa, store_samples);
    let mut topo_map = AHashMap::new();
    let mut topologies = Vec::new();
    let mut tree_records = Vec::new();

    for chunk in per_chunk {
        merge_chunk_accum(
            chunk,
            &mut global,
            &mut topo_map,
            &mut topologies,
            &mut tree_records,
            need_clade_list,
        );
    }

    Ok(SummaryRunData {
        table: global,
        topologies,
        tree_records,
        total_read,
        total_burned_in,
        total_summarized: retained.len(),
        rooted,
        union_ntax: n_taxa,
        is_partial_target: false,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_identical_taxa_files(
    plans: &[FilePlan],
    input_fmt: TreeFmt,
    parse_opts: &ParseOptions,
    taxa: Arc<TaxonTable>,
    union_ntax: usize,
    prune_to_shared: bool,
    rooted: bool,
    options: &TreeSummaryOptions,
    store_samples: bool,
    progress: Option<&dyn ProgressReporter>,
) -> Result<(SummaryRunData, Vec<TreeSummaryFileStats>)> {
    let taxon_index = build_taxon_index(&taxa);
    let n_taxa = taxa.len();
    let zobrist = ZobristTable::new(n_taxa);
    let need_clade_list = options.method == SummaryMethod::Mcc;
    let need_heights = rooted;

    let mut global = GlobalCladeTable::new(Arc::clone(&taxa), store_samples);
    let mut topo_map = AHashMap::new();
    let mut topologies = Vec::new();
    let mut tree_records = Vec::new();
    let mut total_read = 0usize;
    let mut total_burned_in = 0usize;
    let mut total_summarized = 0usize;
    let mut per_file = Vec::with_capacity(plans.len());

    for (file_idx, plan) in plans.iter().enumerate() {
        total_read += plan.total_trees;
        total_burned_in += plan.burnin_skip;

        let count = progress::track(progress, &plan.path, || {
            read_trees_chunked_par(
                &plan.path,
                input_fmt,
                parse_opts,
                plan.burnin_skip,
                CHUNK_SIZE,
                progress,
                |start_idx, mut trees| {
                    if prune_to_shared {
                        trees = trees
                            .into_par_iter()
                            .map(|t| prune_tree_to_taxa(&t, &taxa))
                            .collect::<Result<Vec<_>>>()?;
                    }
                    let sub_chunks: Vec<ChunkAccum> = trees
                        .par_chunks(64)
                        .enumerate()
                        .map(|(sub_idx, slice)| {
                            let base = start_idx + sub_idx * 64 + plan.burnin_skip;
                            process_tree_slice(
                                slice,
                                |local_i| TreeLocation::File {
                                    file_idx,
                                    tree_idx: base + local_i,
                                },
                                &taxa,
                                &taxon_index,
                                &zobrist,
                                rooted,
                                store_samples,
                                options.summarize_annotations,
                                need_heights,
                                need_clade_list,
                            )
                        })
                        .collect::<Result<Vec<_>>>()?;

                    for chunk in sub_chunks {
                        merge_chunk_accum(
                            chunk,
                            &mut global,
                            &mut topo_map,
                            &mut topologies,
                            &mut tree_records,
                            need_clade_list,
                        );
                    }
                    Ok(())
                },
            )
        })?;

        total_summarized += count;
        per_file.push(TreeSummaryFileStats {
            path: plan.path.display().to_string(),
            trees_read: plan.total_trees,
            trees_burned_in: plan.burnin_skip,
            trees_retained: count,
        });
    }

    Ok((
        SummaryRunData {
            table: global,
            topologies,
            tree_records,
            total_read,
            total_burned_in,
            total_summarized,
            rooted,
            union_ntax,
            is_partial_target: false,
        },
        per_file,
    ))
}

struct LocalCladeEntry {
    bits_start: usize,
    hash: (u64, u64),
    stats: CladeStats,
}

struct ChunkAccum {
    arena: CladeBitsArena,
    by_hash: AHashMap<(u64, u64), Vec<u32>>,
    clades: Vec<LocalCladeEntry>,
    tip_stats: Vec<CladeStats>,
    root_stats: CladeStats,
    keys: AnnotationKeyTable,
    tree_clade_lists: Vec<(TreeLocation, Vec<u32>)>,
}

impl ChunkAccum {
    fn new(n_taxa: usize, words: usize, store_samples: bool) -> Self {
        Self {
            arena: CladeBitsArena::new(words),
            by_hash: AHashMap::new(),
            clades: Vec::new(),
            tip_stats: (0..n_taxa)
                .map(|_| CladeStats::new(store_samples))
                .collect(),
            root_stats: CladeStats::new(store_samples),
            keys: AnnotationKeyTable::default(),
            tree_clade_lists: Vec::new(),
        }
    }

    fn get_or_insert(&mut self, bits: &[u64], hash: (u64, u64), store_samples: bool) -> u32 {
        if let Some(bucket) = self.by_hash.get(&hash) {
            for &id in bucket {
                if self.arena.get(CladeId(id)) == bits {
                    return id;
                }
            }
        }
        let cid = self.arena.push(bits);
        self.clades.push(LocalCladeEntry {
            bits_start: cid.index(),
            hash,
            stats: CladeStats::new(store_samples),
        });
        self.by_hash.entry(hash).or_default().push(cid.0);
        cid.0
    }
}

#[allow(clippy::too_many_arguments)]
fn process_tree_slice<F>(
    trees: &[Tree],
    loc_fn: F,
    expected_taxa: &Arc<TaxonTable>,
    taxon_index: &AHashMap<&str, usize>,
    zobrist: &ZobristTable,
    rooted: bool,
    store_samples: bool,
    summarize_annotations: bool,
    need_heights: bool,
    need_clade_list: bool,
) -> Result<ChunkAccum>
where
    F: Fn(usize) -> TreeLocation,
{
    let n_taxa = expected_taxa.len();
    let words = zobrist.words();
    let mut accum = ChunkAccum::new(n_taxa, words, store_samples);
    let mut cached_ptr: Option<*const TaxonTable> = None;
    let mut cached_remap: Vec<Option<usize>> = Vec::new();

    let mut node_bits = Vec::new();
    let mut node_hashes = Vec::new();
    let mut node_heights = Vec::new();
    let mut seen = vec![0u32; n_taxa];
    let mut seen_gen = 0u32;
    let mut key_bits = vec![0u64; words];

    for (local_idx, tree) in trees.iter().enumerate() {
        if tree.is_empty() {
            bail!("Input contains an empty tree");
        }
        seen_gen = seen_gen.wrapping_add(1);
        if seen_gen == 0 {
            seen.fill(0);
            seen_gen = 1;
        }

        // Cache O(1) TaxonId -> canonical index array across trees sharing one TaxonTable.
        let remap = if let Some(arc_taxa) = tree.taxa() {
            let ptr = Arc::as_ptr(arc_taxa);
            if cached_ptr != Some(ptr) {
                cached_remap.clear();
                cached_remap.reserve(arc_taxa.len());
                for name in arc_taxa.iter() {
                    cached_remap.push(taxon_index.get(name).copied());
                }
                cached_ptr = Some(ptr);
            }
            Some(cached_remap.as_slice())
        } else {
            None
        };

        let n_nodes = tree.len();
        node_bits.clear();
        node_bits.resize(n_nodes * words, 0u64);
        node_hashes.clear();
        node_hashes.resize(n_nodes, (0u64, 0u64));

        if need_heights {
            compute_node_heights(tree, &mut node_heights);
        }

        let root = tree.root();
        let root_is_binary = !rooted && tree.child_count(root) == 2;
        let first_rc = if root_is_binary {
            tree.first_child(root)
        } else {
            None
        };
        let second_rc = first_rc.and_then(|c| tree.next_sibling(c));
        let root_combined_bl = if root_is_binary && let (Some(c1), Some(c2)) = (first_rc, second_rc)
        {
            combine_opt_len(tree.branch_length(c1), tree.branch_length(c2))
        } else {
            None
        };
        let one_rc_leaf = root_is_binary
            && first_rc.is_some_and(|c| tree.is_leaf(c))
                != second_rc.is_some_and(|c| tree.is_leaf(c));

        let mut tree_clades = Vec::with_capacity(n_taxa);
        let mut tip_count_seen = 0usize;

        for id in tree.postorder() {
            let idx = id.index();
            if tree.is_leaf(id) {
                let tid = resolve_leaf_taxon_id(tree, id, remap, taxon_index)?;
                if seen[tid] == seen_gen {
                    let name = expected_taxa
                        .get(crate::parser::tree::TaxonId::from_raw(tid as u32))
                        .unwrap_or("<unknown>");
                    bail!("Duplicate taxon '{}' in input tree", name);
                }
                seen[tid] = seen_gen;
                tip_count_seen += 1;

                let base = idx * words;
                set_bit(&mut node_bits[base..base + words], tid);
                node_hashes[idx] = zobrist.taxon_hash(tid);

                let tip = &mut accum.tip_stats[tid];
                tip.count += 1;
                let bl = if (Some(id) == first_rc || Some(id) == second_rc) && one_rc_leaf {
                    root_combined_bl
                } else {
                    tree.branch_length(id)
                };
                if let Some(len) = bl {
                    tip.branch_lengths.push(len);
                }
                if need_heights {
                    tip.heights.push(node_heights[idx]);
                }
                if summarize_annotations {
                    record_node_annots(
                        tree,
                        id,
                        &mut tip.annotations,
                        &mut accum.keys,
                        store_samples,
                    );
                }
                continue;
            }

            let base = idx * words;
            let mut h = (0u64, 0u64);
            for child in tree.children(id) {
                let c_idx = child.index();
                let c_base = c_idx * words;
                for w in 0..words {
                    node_bits[base + w] |= node_bits[c_base + w];
                }
                h.0 ^= node_hashes[c_idx].0;
                h.1 ^= node_hashes[c_idx].1;
            }
            node_hashes[idx] = h;

            if tree.is_root(id) {
                accum.root_stats.count += 1;
                if need_heights {
                    accum.root_stats.heights.push(node_heights[idx]);
                }
                if summarize_annotations {
                    record_node_annots(
                        tree,
                        id,
                        &mut accum.root_stats.annotations,
                        &mut accum.keys,
                        store_samples,
                    );
                }
                continue;
            }

            if root_is_binary && Some(id) == second_rc {
                continue;
            }

            let pop = popcount(&node_bits[base..base + words]);
            let nontrivial = if rooted {
                pop >= 2 && pop < n_taxa
            } else {
                pop >= 2 && pop <= n_taxa.saturating_sub(2)
            };
            if !nontrivial {
                continue;
            }

            key_bits.copy_from_slice(&node_bits[base..base + words]);
            let mut key_hash = h;
            if !rooted {
                canonicalize_unrooted(&mut key_bits, &mut key_hash, n_taxa, zobrist);
            }

            let local_cid = accum.get_or_insert(&key_bits, key_hash, store_samples);
            let entry = &mut accum.clades[local_cid as usize];
            entry.stats.count += 1;
            let bl = if root_is_binary && Some(id) == first_rc {
                root_combined_bl
            } else {
                tree.branch_length(id)
            };
            if let Some(len) = bl {
                entry.stats.branch_lengths.push(len);
            }
            if need_heights {
                entry.stats.heights.push(node_heights[idx]);
            }
            if summarize_annotations {
                record_node_annots(
                    tree,
                    id,
                    &mut entry.stats.annotations,
                    &mut accum.keys,
                    store_samples,
                );
            }
            tree_clades.push(local_cid);
        }

        if tip_count_seen != n_taxa {
            bail!(
                "Taxon set mismatch: tree has {} taxa, expected {} (use --allow-partial-taxa for heterogeneous taxon sets)",
                tip_count_seen,
                n_taxa
            );
        }

        if need_clade_list || !tree_clades.is_empty() || n_taxa <= 3 {
            let loc = loc_fn(local_idx);
            accum.tree_clade_lists.push((loc, tree_clades));
        }
    }

    Ok(accum)
}

fn merge_chunk_accum(
    chunk: ChunkAccum,
    global: &mut GlobalCladeTable,
    topo_map: &mut AHashMap<(u64, u64), Vec<usize>>,
    topologies: &mut Vec<TopologyRecord>,
    tree_records: &mut Vec<TreeCladeRecord>,
    keep_tree_records: bool,
) {
    for (tid, tip) in chunk.tip_stats.into_iter().enumerate() {
        global.tip_stats[tid].merge(tip, &chunk.keys, &mut global.keys);
    }
    global
        .root_stats
        .merge(chunk.root_stats, &chunk.keys, &mut global.keys);

    let mut local_to_global = Vec::with_capacity(chunk.clades.len());
    for entry in chunk.clades {
        let bits = chunk.arena.get(CladeId(entry.bits_start as u32));
        let gid = global.lookup_or_insert(bits, entry.hash);
        global.clades[gid.index()].merge(entry.stats, &chunk.keys, &mut global.keys);
        local_to_global.push(gid);
    }

    for (loc, local_list) in chunk.tree_clade_lists {
        let mut global_list: Vec<CladeId> = local_list
            .into_iter()
            .map(|lid| local_to_global[lid as usize])
            .collect();
        global_list.sort_unstable();
        global_list.dedup();
        let topo_hash = hash_clade_list(&global_list);

        let bucket = topo_map.entry(topo_hash).or_default();
        let mut matched = None;
        for &topo_idx in bucket.iter() {
            if topologies[topo_idx].clade_ids == global_list {
                matched = Some(topo_idx);
                break;
            }
        }
        if let Some(topo_idx) = matched {
            topologies[topo_idx].count += 1;
            if loc < topologies[topo_idx].first_location {
                topologies[topo_idx].first_location = loc;
            }
        } else {
            let new_idx = topologies.len();
            bucket.push(new_idx);
            topologies.push(TopologyRecord {
                hash: topo_hash,
                clade_ids: global_list.clone(),
                count: 1,
                first_location: loc,
            });
        }

        if keep_tree_records {
            tree_records.push(TreeCladeRecord {
                location: loc,
                clade_ids: global_list,
            });
        }
    }
}

#[inline]
fn resolve_leaf_taxon_id(
    tree: &Tree,
    id: NodeId,
    remap: Option<&[Option<usize>]>,
    taxon_index: &AHashMap<&str, usize>,
) -> Result<usize> {
    if let (Some(map), Some(raw_tid)) = (remap, tree.taxon_id(id))
        && let Some(Some(canon)) = map.get(raw_tid.index())
    {
        return Ok(*canon);
    }
    let label = tree
        .label(id)
        .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
    taxon_index.get(label).copied().ok_or_else(|| {
        anyhow::anyhow!(
            "Taxon '{}' was not present in the reference taxon set (use --allow-partial-taxa for heterogeneous taxon sets)",
            label
        )
    })
}

fn build_taxon_table_from_tree(tree: &Tree) -> Result<Arc<TaxonTable>> {
    let mut table = TaxonTable::with_capacity(tree.tip_count());
    let mut seen = AHashSet::with_capacity(tree.tip_count());
    for leaf in tree.leaves() {
        let lbl = tree
            .label(leaf)
            .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
        if !seen.insert(lbl) {
            bail!("Duplicate taxon '{}' in tree", lbl);
        }
        table.intern(lbl);
    }
    if table.len() < 3 {
        bail!(
            "At least 3 taxa are required for tree summary (found {})",
            table.len()
        );
    }
    Ok(Arc::new(table))
}

fn build_taxon_index(taxa: &TaxonTable) -> AHashMap<&str, usize> {
    let mut map = AHashMap::with_capacity(taxa.len());
    for (idx, name) in taxa.iter().enumerate() {
        map.insert(name, idx);
    }
    map
}

fn resolve_rootedness(tree: &Tree, mode: RootedMode) -> Result<bool> {
    match mode {
        RootedMode::Rooted => Ok(true),
        RootedMode::Unrooted => Ok(false),
        RootedMode::Auto => {
            let rooted = match tree.metadata.rootedness {
                Rootedness::Rooted => true,
                Rootedness::Unrooted => false,
                Rootedness::Unknown => tree.child_count(tree.root()) == 2,
            };
            Ok(rooted)
        }
    }
}

fn combine_opt_len(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x + y),
        (Some(x), None) => Some(x),
        (None, Some(y)) => Some(y),
        (None, None) => None,
    }
}

fn compute_node_heights(tree: &Tree, heights: &mut Vec<f64>) {
    heights.clear();
    heights.resize(tree.len(), 0.0);
    let mut depths = vec![0.0f64; tree.len()];
    let mut max_depth = 0.0f64;

    for id in tree.preorder() {
        let idx = id.index();
        let d = if let Some(parent) = tree.parent(id) {
            let bl = tree.branch_length(id).unwrap_or(0.0);
            let bl_nonneg = if bl.is_finite() && bl > 0.0 { bl } else { 0.0 };
            depths[parent.index()] + bl_nonneg
        } else {
            0.0
        };
        depths[idx] = d;
        if tree.is_leaf(id) && d > max_depth {
            max_depth = d;
        }
    }

    for (idx, &d) in depths.iter().enumerate() {
        heights[idx] = (max_depth - d).max(0.0);
    }
}

fn record_node_annots(
    tree: &Tree,
    id: NodeId,
    acc: &mut CladeAnnotationsAcc,
    keys: &mut AnnotationKeyTable,
    store_samples: bool,
) {
    for annot in tree.node_annotations(id) {
        acc.record_annotation(annot, keys, store_samples);
    }
    for annot in tree.branch_annotations(id) {
        acc.record_annotation(annot, keys, store_samples);
    }
}

// ---------- Partial-taxa support: shared-taxon intersection & decisive-tree gCF ----------

fn discover_shared_taxa_in_memory(
    trees: &[Tree],
    target: Option<&Tree>,
) -> Result<(Arc<TaxonTable>, usize)> {
    let mut union = AHashSet::new();
    let mut shared: Option<AHashSet<String>> = None;

    let mut update = |t: &Tree| -> Result<()> {
        let mut cur = AHashSet::new();
        for leaf in t.leaves() {
            let lbl = t
                .label(leaf)
                .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
            union.insert(lbl.to_string());
            cur.insert(lbl.to_string());
        }
        match &mut shared {
            None => shared = Some(cur),
            Some(prev) => prev.retain(|name| cur.contains(name)),
        }
        Ok(())
    };

    if let Some(tgt) = target {
        update(tgt)?;
    }
    for t in trees {
        update(t)?;
    }

    finalize_shared_taxa(shared.unwrap_or_default(), union.len())
}

fn discover_shared_taxa_files(
    plans: &[FilePlan],
    input_fmt: TreeFmt,
    parse_opts: &ParseOptions,
    target: Option<&Tree>,
    progress: Option<&dyn ProgressReporter>,
) -> Result<(Arc<TaxonTable>, usize)> {
    let mut union = AHashSet::new();
    let mut shared: Option<AHashSet<String>> = None;

    let mut observe = |t: &Tree| -> Result<()> {
        let mut cur = AHashSet::new();
        for leaf in t.leaves() {
            let lbl = t
                .label(leaf)
                .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
            union.insert(lbl.to_string());
            cur.insert(lbl.to_string());
        }
        match &mut shared {
            None => shared = Some(cur),
            Some(prev) => prev.retain(|name| cur.contains(name)),
        }
        Ok(())
    };

    if let Some(tgt) = target {
        observe(tgt)?;
    }

    for plan in plans {
        read_trees_chunked_par(
            &plan.path,
            input_fmt,
            parse_opts,
            plan.burnin_skip,
            CHUNK_SIZE,
            progress,
            |_, trees| {
                for t in &trees {
                    observe(t)?;
                }
                Ok(())
            },
        )?;
    }

    finalize_shared_taxa(shared.unwrap_or_default(), union.len())
}

fn finalize_shared_taxa(
    shared: AHashSet<String>,
    union_ntax: usize,
) -> Result<(Arc<TaxonTable>, usize)> {
    if shared.len() < 3 {
        bail!(
            "Shared taxon intersection has only {} taxon/taxa (union has {}); at least 3 shared taxa are required for consensus/MCC/MAP. Use `--method target --target <FILE>` for decisive-tree partial-taxa support.",
            shared.len(),
            union_ntax
        );
    }
    let mut names: Vec<String> = shared.into_iter().collect();
    alphanumeric_sort::sort_str_slice(&mut names);
    let mut table = TaxonTable::with_capacity(names.len());
    for n in names {
        table.intern(&n);
    }
    Ok((Arc::new(table), union_ntax))
}

/// Iteratively prune `tree` down to `keep_taxa`, summing branch lengths across suppressed nodes.
pub(crate) fn prune_tree_to_taxa(tree: &Tree, keep_taxa: &TaxonTable) -> Result<Tree> {
    let n = tree.len();
    let mut keep_node = vec![false; n];
    let mut sub_count = vec![0usize; n];

    for id in tree.postorder() {
        let idx = id.index();
        if tree.is_leaf(id) {
            let lbl = tree
                .label(id)
                .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
            if keep_taxa.get_index_of(lbl).is_some() {
                keep_node[idx] = true;
                sub_count[idx] = 1;
            }
        } else {
            let mut active_children = 0usize;
            for c in tree.children(id) {
                if sub_count[c.index()] > 0 {
                    active_children += 1;
                    sub_count[idx] += sub_count[c.index()];
                }
            }
            if active_children >= 2 {
                keep_node[idx] = true;
            }
        }
    }

    // Locate effective root (highest kept node with sub_count > 0).
    let mut eff_root = None;
    for id in tree.preorder() {
        if keep_node[id.index()] {
            eff_root = Some(id);
            break;
        }
    }
    let eff_root = eff_root.context("No shared taxa found in tree during pruning")?;

    let mut out = Tree::with_capacity(n, 0);
    out.metadata = tree.metadata.clone();
    let out_root = out.add_node();
    out.set_root(out_root);

    // Map kept source nodes -> new output nodes, accumulating branch lengths across unary runs.
    let mut stack: Vec<(NodeId, NodeId, f64, bool)> = Vec::new();
    for c in tree.children(eff_root) {
        if sub_count[c.index()] > 0 {
            let bl = tree.branch_length(c);
            stack.push((c, out_root, bl.unwrap_or(0.0), bl.is_some()));
        }
    }

    let mut parent_last_child: AHashMap<NodeId, NodeId> = AHashMap::new();
    while let Some((src, out_parent, acc_bl, has_bl)) = stack.pop() {
        if !keep_node[src.index()] {
            for c in tree.children(src) {
                if sub_count[c.index()] > 0 {
                    let c_bl = tree.branch_length(c);
                    stack.push((
                        c,
                        out_parent,
                        acc_bl + c_bl.unwrap_or(0.0),
                        has_bl || c_bl.is_some(),
                    ));
                }
            }
            continue;
        }

        let new_id = out.add_node();
        if let Some(lbl) = tree.label(src) {
            out.set_label(new_id, lbl, tree.is_quoted_label(src));
        }
        if has_bl {
            out.set_branch_length(new_id, acc_bl);
        }
        for annot in tree.node_annotations(src) {
            out.add_annotation(new_id, AnnotationTarget::Node, annot.clone());
        }
        for annot in tree.branch_annotations(src) {
            out.add_annotation(new_id, AnnotationTarget::Branch, annot.clone());
        }
        let prev = parent_last_child.get(&out_parent).copied();
        out.attach_child(out_parent, new_id, prev);
        parent_last_child.insert(out_parent, new_id);

        let mut children: Vec<NodeId> = tree
            .children(src)
            .filter(|&c| sub_count[c.index()] > 0)
            .collect();
        children.reverse();
        for c in children {
            let c_bl = tree.branch_length(c);
            stack.push((c, new_id, c_bl.unwrap_or(0.0), c_bl.is_some()));
        }
    }

    Ok(out)
}

fn run_partial_target_in_memory(
    retained: &[Tree],
    target: &Tree,
    rooted: bool,
    options: &TreeSummaryOptions,
    total_read: usize,
    total_burned_in: usize,
) -> Result<SummaryRunData> {
    let store_samples = needs_sample_storage(options);
    let taxa = build_taxon_table_from_tree(target)?;
    let mut table = GlobalCladeTable::new(Arc::clone(&taxa), store_samples);
    let (_, _) = map_tree_nodes_to_clades(target, &mut table, rooted, true)?;
    let target_index = build_taxon_index(&taxa);

    let per_chunk: Vec<PartialTargetChunk> = retained
        .par_chunks(CHUNK_SIZE)
        .map(|chunk| {
            process_partial_target_slice(
                chunk,
                &table,
                &target_index,
                rooted,
                store_samples,
                options.summarize_annotations,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let mut union_names: AHashSet<String> = taxa.iter().map(str::to_string).collect();
    for c in per_chunk {
        union_names.extend(c.extra_taxa);
        for (tid, tip) in c.tip_stats.into_iter().enumerate() {
            table.tip_stats[tid].merge(tip, &c.keys, &mut table.keys);
        }
        table
            .root_stats
            .merge(c.root_stats, &c.keys, &mut table.keys);
        for (cid, stats) in c.clade_stats.into_iter().enumerate() {
            table.clades[cid].merge(stats, &c.keys, &mut table.keys);
        }
    }

    Ok(SummaryRunData {
        table,
        topologies: Vec::new(),
        tree_records: Vec::new(),
        total_read,
        total_burned_in,
        total_summarized: retained.len(),
        rooted,
        union_ntax: union_names.len(),
        is_partial_target: true,
    })
}

fn run_partial_target_files(
    plans: &[FilePlan],
    input_fmt: TreeFmt,
    parse_opts: &ParseOptions,
    target: &Tree,
    rooted: bool,
    options: &TreeSummaryOptions,
    progress: Option<&dyn ProgressReporter>,
) -> Result<(SummaryRunData, Vec<TreeSummaryFileStats>)> {
    let store_samples = needs_sample_storage(options);
    let taxa = build_taxon_table_from_tree(target)?;
    let mut table = GlobalCladeTable::new(Arc::clone(&taxa), store_samples);
    let (_, _) = map_tree_nodes_to_clades(target, &mut table, rooted, true)?;
    let target_index = build_taxon_index(&taxa);
    let mut union_names: AHashSet<String> = taxa.iter().map(str::to_string).collect();

    let mut total_read = 0usize;
    let mut total_burned_in = 0usize;
    let mut total_summarized = 0usize;
    let mut per_file = Vec::with_capacity(plans.len());

    for plan in plans {
        total_read += plan.total_trees;
        total_burned_in += plan.burnin_skip;

        let count = progress::track(progress, &plan.path, || {
            read_trees_chunked_par(
                &plan.path,
                input_fmt,
                parse_opts,
                plan.burnin_skip,
                CHUNK_SIZE,
                progress,
                |_, trees| {
                    let sub_chunks: Vec<PartialTargetChunk> = trees
                        .par_chunks(64)
                        .map(|slice| {
                            process_partial_target_slice(
                                slice,
                                &table,
                                &target_index,
                                rooted,
                                store_samples,
                                options.summarize_annotations,
                            )
                        })
                        .collect::<Result<Vec<_>>>()?;

                    for c in sub_chunks {
                        union_names.extend(c.extra_taxa);
                        for (tid, tip) in c.tip_stats.into_iter().enumerate() {
                            table.tip_stats[tid].merge(tip, &c.keys, &mut table.keys);
                        }
                        table
                            .root_stats
                            .merge(c.root_stats, &c.keys, &mut table.keys);
                        for (cid, stats) in c.clade_stats.into_iter().enumerate() {
                            table.clades[cid].merge(stats, &c.keys, &mut table.keys);
                        }
                    }
                    Ok(())
                },
            )
        })?;

        total_summarized += count;
        per_file.push(TreeSummaryFileStats {
            path: plan.path.display().to_string(),
            trees_read: plan.total_trees,
            trees_burned_in: plan.burnin_skip,
            trees_retained: count,
        });
    }

    Ok((
        SummaryRunData {
            table,
            topologies: Vec::new(),
            tree_records: Vec::new(),
            total_read,
            total_burned_in,
            total_summarized,
            rooted,
            union_ntax: union_names.len(),
            is_partial_target: true,
        },
        per_file,
    ))
}

type SplitMap = AHashMap<(u64, u64), Vec<(Vec<u64>, NodeId)>>;

struct PartialTargetChunk {
    clade_stats: Vec<CladeStats>,
    tip_stats: Vec<CladeStats>,
    root_stats: CladeStats,
    keys: AnnotationKeyTable,
    extra_taxa: AHashSet<String>,
}

fn process_partial_target_slice(
    trees: &[Tree],
    table: &GlobalCladeTable,
    target_index: &AHashMap<&str, usize>,
    rooted: bool,
    store_samples: bool,
    summarize_annotations: bool,
) -> Result<PartialTargetChunk> {
    let n_target_taxa = table.taxa.len();
    let words = table.zobrist.words();
    let mut clade_stats = (0..table.clades.len())
        .map(|_| CladeStats::new(store_samples))
        .collect::<Vec<_>>();
    let mut tip_stats = (0..n_target_taxa)
        .map(|_| CladeStats::new(store_samples))
        .collect::<Vec<_>>();
    let mut root_stats = CladeStats::new(store_samples);
    let mut keys = AnnotationKeyTable::default();
    let mut extra_taxa = AHashSet::new();

    let mut node_bits = Vec::new();
    let mut node_heights = Vec::new();
    let mut leaf_mask = vec![0u64; words];
    let mut restricted = vec![0u64; words];

    for tree in trees {
        if tree.is_empty() {
            bail!("Input contains an empty tree");
        }
        let n_nodes = tree.len();
        node_bits.clear();
        node_bits.resize(n_nodes * words, 0u64);
        leaf_mask.fill(0);
        if rooted {
            compute_node_heights(tree, &mut node_heights);
        }

        // Collect every canonical split present in this tree (projected into target taxon space).
        let mut tree_splits: SplitMap = AHashMap::new();
        for id in tree.postorder() {
            let idx = id.index();
            if tree.is_leaf(id) {
                let lbl = tree
                    .label(id)
                    .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
                if let Some(&tid) = target_index.get(lbl) {
                    set_bit(&mut node_bits[idx * words..(idx + 1) * words], tid);
                    set_bit(&mut leaf_mask, tid);
                    let tip = &mut tip_stats[tid];
                    tip.count += 1;
                    tip.decisive_count += 1;
                    if let Some(len) = tree.branch_length(id) {
                        tip.branch_lengths.push(len);
                    }
                    if rooted {
                        tip.heights.push(node_heights[idx]);
                    }
                    if summarize_annotations {
                        record_node_annots(
                            tree,
                            id,
                            &mut tip.annotations,
                            &mut keys,
                            store_samples,
                        );
                    }
                } else {
                    extra_taxa.insert(lbl.to_string());
                }
            } else {
                let base = idx * words;
                for child in tree.children(id) {
                    let c_base = child.index() * words;
                    for w in 0..words {
                        node_bits[base + w] |= node_bits[c_base + w];
                    }
                }
                if tree.is_root(id) {
                    root_stats.count += 1;
                    root_stats.decisive_count += 1;
                    if rooted {
                        root_stats.heights.push(node_heights[idx]);
                    }
                    if summarize_annotations {
                        record_node_annots(
                            tree,
                            id,
                            &mut root_stats.annotations,
                            &mut keys,
                            store_samples,
                        );
                    }
                }
            }
        }

        let total_shared = popcount(&leaf_mask);
        if total_shared < 3 {
            continue;
        }

        let anchor = first_set_bit(&leaf_mask);
        for id in tree.postorder() {
            if tree.is_leaf(id) || tree.is_root(id) {
                continue;
            }
            let base = id.index() * words;
            let mut bits = node_bits[base..base + words].to_vec();
            let pop = popcount(&bits);
            let nontrivial = if rooted {
                pop >= 2 && pop < total_shared
            } else {
                pop >= 2 && pop <= total_shared.saturating_sub(2)
            };
            if !nontrivial {
                continue;
            }
            if !rooted
                && let Some(anc) = anchor
                && bit_is_set(&bits, anc)
            {
                for (b, &m) in bits.iter_mut().zip(leaf_mask.iter()) {
                    *b = (!*b) & m;
                }
            }
            let hash = table.zobrist.hash_bits(&bits);
            let bucket = tree_splits.entry(hash).or_default();
            if !bucket.iter().any(|(existing, _)| existing == &bits) {
                bucket.push((bits, id));
            }
        }

        // Check every target clade against `leaf_mask` for decisiveness and support.
        for (cid, stats) in clade_stats.iter_mut().enumerate() {
            let target_bits = table.arena.get(CladeId(cid as u32));
            restricted.copy_from_slice(target_bits);
            let mut inside = 0usize;
            for (r, &m) in restricted.iter_mut().zip(leaf_mask.iter()) {
                *r &= m;
                inside += r.count_ones() as usize;
            }
            let outside = total_shared - inside;
            let decisive = if rooted {
                inside >= 2 && outside >= 1
            } else {
                inside >= 2 && outside >= 2
            };
            if !decisive {
                continue;
            }
            stats.decisive_count += 1;
            canonicalize_for_mask(&mut restricted, &leaf_mask, rooted);
            let r_hash = table.zobrist.hash_bits(&restricted);
            if let Some(bucket) = tree_splits.get(&r_hash)
                && let Some((_, matched_node)) = bucket
                    .iter()
                    .find(|(b, _)| b.as_slice() == restricted.as_slice())
            {
                stats.count += 1;
                if let Some(len) = tree.branch_length(*matched_node) {
                    stats.branch_lengths.push(len);
                }
                if rooted {
                    stats.heights.push(node_heights[matched_node.index()]);
                }
                if summarize_annotations {
                    record_node_annots(
                        tree,
                        *matched_node,
                        &mut stats.annotations,
                        &mut keys,
                        store_samples,
                    );
                }
            }
        }
    }

    Ok(PartialTargetChunk {
        clade_stats,
        tip_stats,
        root_stats,
        keys,
        extra_taxa,
    })
}

// ---------- Pass 2: O(N log N) Euler-Tour + Sparse-Table RMQ Common Ancestor Heights ----------

struct SummaryGuide<'a> {
    tree: &'a Tree,
    postorder: Vec<NodeId>,
    leaf_taxon: Vec<Option<usize>>,
    node_clade_ids: Vec<Option<CladeId>>,
    root: NodeId,
    n_nodes: usize,
    n_clades: usize,
    n_taxa: usize,
}

impl<'a> SummaryGuide<'a> {
    fn new(
        summary: &'a Tree,
        node_clade_ids: &[Option<CladeId>],
        taxa: &TaxonTable,
        n_clades: usize,
    ) -> Result<Self> {
        let n_nodes = summary.len();
        let mut leaf_taxon = vec![None; n_nodes];
        for id in summary.leaves() {
            let lbl = summary
                .label(id)
                .ok_or_else(|| anyhow::anyhow!("Summary tree has unlabeled tip"))?;
            let tid = taxa
                .get_index_of(lbl)
                .map(|t| t.index())
                .ok_or_else(|| anyhow::anyhow!("Summary tip '{}' missing from TaxonTable", lbl))?;
            leaf_taxon[id.index()] = Some(tid);
        }
        Ok(Self {
            tree: summary,
            postorder: summary.postorder().collect(),
            leaf_taxon,
            node_clade_ids: node_clade_ids.to_vec(),
            root: summary.root(),
            n_nodes,
            n_clades,
            n_taxa: taxa.len(),
        })
    }
}

struct CaChunkAccum {
    clade_ca: Vec<MetricAccumulator>,
    tip_ca: Vec<MetricAccumulator>,
    root_ca: MetricAccumulator,
}

impl CaChunkAccum {
    fn new(n_clades: usize, n_taxa: usize, store_samples: bool) -> Self {
        Self {
            clade_ca: (0..n_clades)
                .map(|_| MetricAccumulator::new(store_samples))
                .collect(),
            tip_ca: (0..n_taxa)
                .map(|_| MetricAccumulator::new(store_samples))
                .collect(),
            root_ca: MetricAccumulator::new(store_samples),
        }
    }
}

/// Populate `ca_heights` on `table` in $O(N \log N)$ per input tree (in-memory slice).
pub(crate) fn compute_ca_heights_in_memory(
    trees: &[Tree],
    burnin_skip: usize,
    summary: &Tree,
    node_clade_ids: &[Option<CladeId>],
    table: &mut GlobalCladeTable,
    store_samples: bool,
) -> Result<()> {
    let retained = &trees[burnin_skip..];
    let guide = SummaryGuide::new(summary, node_clade_ids, &table.taxa, table.clades.len())?;
    let taxon_index = build_taxon_index(&table.taxa);

    let chunks: Vec<CaChunkAccum> = retained
        .par_chunks(CHUNK_SIZE)
        .map(|slice| compute_ca_slice(slice, &guide, &taxon_index, store_samples))
        .collect::<Result<Vec<_>>>()?;

    merge_ca_chunks(chunks, table);
    Ok(())
}

/// Populate `ca_heights` on `table` in $O(N \log N)$ per input tree across files.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_ca_heights_files(
    plans: &[FilePlan],
    input_fmt: TreeFmt,
    parse_opts: &ParseOptions,
    summary: &Tree,
    node_clade_ids: &[Option<CladeId>],
    table: &mut GlobalCladeTable,
    store_samples: bool,
    prune_to_shared: bool,
    progress: Option<&dyn ProgressReporter>,
) -> Result<()> {
    let taxa_arc = Arc::clone(&table.taxa);
    let guide = SummaryGuide::new(summary, node_clade_ids, &taxa_arc, table.clades.len())?;
    let taxon_index = build_taxon_index(&taxa_arc);

    for plan in plans {
        read_trees_chunked_par(
            &plan.path,
            input_fmt,
            parse_opts,
            plan.burnin_skip,
            CHUNK_SIZE,
            progress,
            |_, mut trees| {
                if prune_to_shared {
                    trees = trees
                        .into_par_iter()
                        .map(|t| prune_tree_to_taxa(&t, &taxa_arc))
                        .collect::<Result<Vec<_>>>()?;
                }
                let chunks: Vec<CaChunkAccum> = trees
                    .par_chunks(64)
                    .map(|slice| compute_ca_slice(slice, &guide, &taxon_index, store_samples))
                    .collect::<Result<Vec<_>>>()?;
                merge_ca_chunks(chunks, table);
                Ok(())
            },
        )?;
    }

    Ok(())
}

fn merge_ca_chunks(chunks: Vec<CaChunkAccum>, table: &mut GlobalCladeTable) {
    for c in chunks {
        table.root_stats.ca_heights.merge(c.root_ca);
        for (i, acc) in c.tip_ca.into_iter().enumerate() {
            table.tip_stats[i].ca_heights.merge(acc);
        }
        for (i, acc) in c.clade_ca.into_iter().enumerate() {
            table.clades[i].ca_heights.merge(acc);
        }
    }
}

fn compute_ca_slice(
    trees: &[Tree],
    guide: &SummaryGuide,
    taxon_index: &AHashMap<&str, usize>,
    store_samples: bool,
) -> Result<CaChunkAccum> {
    let mut accum = CaChunkAccum::new(guide.n_clades, guide.n_taxa, store_samples);
    let mut cached_ptr: Option<*const TaxonTable> = None;
    let mut cached_remap: Vec<Option<usize>> = Vec::new();

    let mut euler_nodes = Vec::new();
    let mut euler_depth = Vec::new();
    let mut tin = vec![u32::MAX; guide.n_taxa];
    let mut heights = Vec::new();
    let mut sparse = Vec::new();
    let mut log2_tab = Vec::new();
    let mut lo = vec![u32::MAX; guide.n_nodes];
    let mut hi = vec![0u32; guide.n_nodes];

    for tree in trees {
        let remap = if let Some(arc_taxa) = tree.taxa() {
            let ptr = Arc::as_ptr(arc_taxa);
            if cached_ptr != Some(ptr) {
                cached_remap.clear();
                cached_remap.reserve(arc_taxa.len());
                for name in arc_taxa.iter() {
                    cached_remap.push(taxon_index.get(name).copied());
                }
                cached_ptr = Some(ptr);
            }
            Some(cached_remap.as_slice())
        } else {
            None
        };

        compute_node_heights(tree, &mut heights);
        build_euler_rmq(
            tree,
            remap,
            taxon_index,
            &mut euler_nodes,
            &mut euler_depth,
            &mut tin,
            &mut sparse,
            &mut log2_tab,
        )?;

        let m = euler_nodes.len();
        if m == 0 {
            continue;
        }

        // Single O(|S|) bottom-up pass computing (min_tin, max_tin) and O(1) RMQ MRCA.
        for &sid in &guide.postorder {
            let s_idx = sid.index();
            if let Some(tid) = guide.leaf_taxon[s_idx] {
                let pos = tin[tid];
                lo[s_idx] = pos;
                hi[s_idx] = pos;
                if pos != u32::MAX {
                    let leaf_node = euler_nodes[pos as usize];
                    accum.tip_ca[tid].push(heights[leaf_node.index()]);
                }
                continue;
            }

            let mut min_t = u32::MAX;
            let mut max_t = 0u32;
            let mut any = false;
            for child in guide.tree.children(sid) {
                let c_lo = lo[child.index()];
                if c_lo == u32::MAX {
                    continue;
                }
                let c_hi = hi[child.index()];
                if !any || c_lo < min_t {
                    min_t = c_lo;
                }
                if !any || c_hi > max_t {
                    max_t = c_hi;
                }
                any = true;
            }
            lo[s_idx] = min_t;
            hi[s_idx] = if any { max_t } else { 0 };
            if !any {
                continue;
            }
            let lca_node = query_rmq_lca(
                min_t as usize,
                max_t as usize,
                m,
                &euler_nodes,
                &euler_depth,
                &sparse,
                &log2_tab,
            );
            let h = heights[lca_node.index()];
            if sid == guide.root {
                accum.root_ca.push(h);
            } else if let Some(cid) = guide.node_clade_ids[s_idx] {
                accum.clade_ca[cid.index()].push(h);
            }
        }
    }

    Ok(accum)
}

#[allow(clippy::too_many_arguments)]
fn build_euler_rmq(
    tree: &Tree,
    remap: Option<&[Option<usize>]>,
    taxon_index: &AHashMap<&str, usize>,
    euler_nodes: &mut Vec<NodeId>,
    euler_depth: &mut Vec<u32>,
    tin: &mut [u32],
    sparse: &mut Vec<u32>,
    log2_tab: &mut Vec<usize>,
) -> Result<()> {
    struct Frame {
        node: NodeId,
        next_child: Option<NodeId>,
        depth: u32,
        visited: bool,
    }

    euler_nodes.clear();
    euler_depth.clear();
    tin.fill(u32::MAX);

    if tree.is_empty() {
        return Ok(());
    }

    let mut stack = Vec::with_capacity(64);
    let root = tree.root();
    stack.push(Frame {
        node: root,
        next_child: tree.first_child(root),
        depth: 0,
        visited: false,
    });

    while let Some(top) = stack.last_mut() {
        let cur_idx = euler_nodes.len();
        euler_nodes.push(top.node);
        euler_depth.push(top.depth);

        if !top.visited {
            top.visited = true;
            if tree.is_leaf(top.node) {
                if let Some(canon_tid) =
                    resolve_optional_leaf_taxon_id(tree, top.node, remap, taxon_index)
                {
                    tin[canon_tid] = cur_idx as u32;
                }
                stack.pop();
                continue;
            }
        }

        if let Some(child) = top.next_child {
            let child_depth = top.depth + 1;
            top.next_child = tree.next_sibling(child);
            stack.push(Frame {
                node: child,
                next_child: tree.first_child(child),
                depth: child_depth,
                visited: false,
            });
        } else {
            stack.pop();
        }
    }

    let m = euler_nodes.len();
    if log2_tab.len() <= m {
        log2_tab.resize(m + 1, 0);
        for i in 2..=m {
            log2_tab[i] = log2_tab[i / 2] + 1;
        }
    }
    let levels = log2_tab[m] + 1;
    sparse.clear();
    sparse.resize(levels * m, 0);
    for (i, slot) in sparse[..m].iter_mut().enumerate() {
        *slot = i as u32;
    }

    for k in 1..levels {
        let span = 1usize << (k - 1);
        let prev_row = (k - 1) * m;
        let cur_row = k * m;
        let limit = m + 1 - (1usize << k);
        for i in 0..limit {
            let a = sparse[prev_row + i];
            let b = sparse[prev_row + i + span];
            sparse[cur_row + i] = if euler_depth[a as usize] <= euler_depth[b as usize] {
                a
            } else {
                b
            };
        }
    }

    Ok(())
}

#[inline]
fn query_rmq_lca(
    l: usize,
    r: usize,
    m: usize,
    euler_nodes: &[NodeId],
    euler_depth: &[u32],
    sparse: &[u32],
    log2_tab: &[usize],
) -> NodeId {
    let (left, right) = if l <= r { (l, r) } else { (r, l) };
    let len = right - left + 1;
    let k = log2_tab[len];
    let row = k * m;
    let a = sparse[row + left];
    let b = sparse[row + right + 1 - (1usize << k)];
    let best = if euler_depth[a as usize] <= euler_depth[b as usize] {
        a
    } else {
        b
    };
    euler_nodes[best as usize]
}

#[inline]
fn resolve_optional_leaf_taxon_id(
    tree: &Tree,
    id: NodeId,
    remap: Option<&[Option<usize>]>,
    taxon_index: &AHashMap<&str, usize>,
) -> Option<usize> {
    if let (Some(map), Some(raw_tid)) = (remap, tree.taxon_id(id))
        && let Some(opt) = map.get(raw_tid.index())
    {
        return *opt;
    }
    let lbl = tree.label(id)?;
    taxon_index.get(lbl).copied()
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use crate::parser::tree::parse_tree_str;

    use super::*;

    #[test]
    fn test_euler_rmq_ca_heights_match_bruteforce_mrca() {
        let t1 = parse_tree_str("(((A:1.0,B:1.0):2.0,C:3.0):1.0,D:4.0);").unwrap();
        let t2 = parse_tree_str("((A:1.5,C:1.5):1.5,(B:0.5,D:0.5):2.5);").unwrap();
        // Summary tree has clade {A, B, C}:
        // In t1, MRCA(A,B,C) has height 3.0.
        // In t2, MRCA(A,B,C) is the root with height 3.0.
        // Clade {A, B}:
        // In t1, MRCA(A,B) has height 1.0.
        // In t2, MRCA(A,B) is the root with height 3.0 -> CA mean = 2.0.
        let summary = parse_tree_str("(((A:1.0,B:1.0):1.0,C:2.0):1.0,D:3.0);").unwrap();
        let taxa = build_taxon_table_from_tree(&summary).unwrap();
        let mut table = GlobalCladeTable::new(Arc::clone(&taxa), true);
        let (node_clades, _) = map_tree_nodes_to_clades(&summary, &mut table, true, true).unwrap();

        let trees = [t1, t2];
        compute_ca_heights_in_memory(&trees, 0, &summary, &node_clades, &mut table, true).unwrap();

        let mut ca_by_size: AHashMap<usize, f64> = AHashMap::new();
        for (idx, stats) in table.clades.iter().enumerate() {
            let pop = popcount(table.arena.get(CladeId(idx as u32)));
            ca_by_size.insert(pop, stats.ca_heights.mean().unwrap());
        }
        assert_approx_eq!(ca_by_size[&2], 2.0, 1e-9);
        assert_approx_eq!(ca_by_size[&3], 3.0, 1e-9);
        assert_approx_eq!(table.root_stats.ca_heights.mean().unwrap(), 3.5, 1e-9);
    }
}
