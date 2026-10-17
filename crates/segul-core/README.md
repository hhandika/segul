# segul-core

Library crate behind [SEGUL](https://www.segul.app/), an ultrafast and memory-efficient tool for phylogenomics. It provides the parsers, writers, summary statistics and core operations used by the `segul` command-line app.

If you only need the command-line app, install [`segul`](https://crates.io/crates/segul) instead.

## Usage

All features are enabled by default. To pull in only what you need:

```toml
[dependencies]
segul-core = { version = "0.24", default-features = false, features = ["align", "sequence"] }
```

| Feature    | Provides                                                           |
|------------|--------------------------------------------------------------------|
| `align`    | Alignment concat, convert, distance, filter, split, summary, trim  |
| `sequence` | Sequence extract, filter, rename, remove, translate, transpose, ID |
| `read`     | Raw read (FASTQ) summary statistics                                |
| `contig`   | Contiguous sequence summary statistics                             |
| `maf`      | Multiple alignment format (MAF) conversion                         |

The base modules (`helper`, the FASTA/NEXUS/PHYLIP parsers, and sequence and partition writers) are always available.

## Migrating from `segul` ≤ 0.23

Replace `segul::` with `segul_core::` in your imports. `segul` 0.24 still re-exports these modules for backward compatibility, but the re-exports will be removed in a future release.

## License

MIT
