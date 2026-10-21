"""Seeded synthetic datasets and real dataset lookup.

Synthetic data is deterministic: the same generator version, dataset id and
parameters always produce the same files, so results from different runs and
machines are comparable. Each dataset is generated once into
`.cache/data/<scale>/<id>/` and reused while its manifest matches.
"""

from __future__ import annotations

import json
import random
import shutil
import zlib
from pathlib import Path

from .config import CACHE_DIR, Dataset

# Bump when the output of a generator changes, to regenerate cached data.
GENERATOR_VERSION = 1
BASES = b"ACGT"
EXTENSIONS = {"nexus": "nex", "phylip": "phy", "fasta": "fas"}


def ensure(dataset: Dataset, scale: str) -> Path | None:
    """Return the directory holding the dataset files, generating synthetic
    data if needed. Returns None for a real dataset whose path is missing."""
    if dataset.kind == "real":
        if dataset.path is None or not dataset.path.exists():
            return None
        return dataset.path

    params = dataset.params.get(scale)
    if params is None:
        raise SystemExit(f"Dataset {dataset.id} has no parameters for scale '{scale}'")
    root = CACHE_DIR / "data" / scale / dataset.id
    files_dir = root / "files"
    manifest_path = root / "manifest.json"
    manifest = {
        "generator": dataset.generator,
        "generator_version": GENERATOR_VERSION,
        "seed": seed_for(dataset.id),
        "params": params,
    }
    if manifest_path.exists() and files_dir.exists():
        cached = json.loads(manifest_path.read_text())
        if {k: cached.get(k) for k in manifest} == manifest:
            return files_dir

    print(f"Generating dataset {dataset.id} ({scale})")
    if root.exists():
        shutil.rmtree(root)
    files_dir.mkdir(parents=True)
    rng = random.Random(manifest["seed"])
    generator = GENERATORS.get(dataset.generator or "")
    if generator is None:
        raise SystemExit(f"Unknown generator {dataset.generator} for dataset {dataset.id}")
    generator(files_dir, rng, params)
    manifest["size_bytes"] = dir_size(files_dir)
    manifest["file_count"] = sum(1 for p in files_dir.iterdir() if p.is_file())
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return files_dir


def describe(dataset: Dataset, path: Path | None) -> dict:
    info = {
        "id": dataset.id,
        "name": dataset.name,
        "kind": dataset.kind,
        "datatype": dataset.datatype,
        "description": dataset.description,
    }
    if dataset.bases:
        info["bases"] = dataset.bases
    if path is not None:
        # Top-level files only: real dataset folders may hold copies in subfolders.
        files = [p for p in path.iterdir() if p.is_file() and not p.name.startswith(".")]
        info["size_mb"] = round(sum(p.stat().st_size for p in files) / 1e6, 1)
        info["files"] = len(files)
    return info


def seed_for(dataset_id: str) -> int:
    return zlib.crc32(dataset_id.encode())


def dir_size(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def make_pool(rng: random.Random, size: int) -> bytes:
    """Random bases to slice sequences from, which is much faster
    than drawing every base."""
    return bytes(rng.choices(BASES, k=size))


def mutate(rng: random.Random, seq: bytes, rate: float, gap_rate: float = 0.0) -> bytes:
    """A copy of `seq` with a fraction of sites changed, so taxa differ."""
    out = bytearray(seq)
    n = len(out)
    for pos in rng.sample(range(n), int(n * rate)):
        out[pos] = rng.choice(BASES)
    for pos in rng.sample(range(n), int(n * gap_rate)):
        out[pos] = ord("-")
    return bytes(out)


def taxon_names(count: int) -> list[bytes]:
    return [f"taxon_{i:04d}".encode() for i in range(count)]


# Writers ----------------------------------------------------------------


def write_alignment(path: Path, fmt: str, names: list[bytes], seqs: list[bytes], layout: str, width: int) -> None:
    nchar = max(len(s) for s in seqs)
    with open(path, "wb") as f:
        if fmt == "nexus":
            f.write(b"#NEXUS\nbegin data;\n")
            f.write(f"dimensions ntax={len(names)} nchar={nchar};\n".encode())
            interleave = b" interleave" if layout == "interleaved" else b""
            f.write(b"format datatype=dna missing=? gap=-" + interleave + b";\nmatrix\n")
            write_rows(f, names, seqs, layout, width, ids_every_block=True)
            f.write(b";\nend;\n")
        elif fmt == "phylip":
            f.write(f"{len(names)} {nchar}\n".encode())
            write_rows(f, names, seqs, layout, width, ids_every_block=False)
        elif fmt == "fasta":
            for name, seq in zip(names, seqs):
                f.write(b">" + name + b"\n")
                if layout == "wrapped":
                    for start in range(0, len(seq), width):
                        f.write(seq[start : start + width] + b"\n")
                else:
                    f.write(seq + b"\n")
        else:
            raise SystemExit(f"Unknown format {fmt}")


def write_rows(f, names, seqs, layout, width, ids_every_block):
    if layout != "interleaved":
        for name, seq in zip(names, seqs):
            f.write(name + b" " + seq + b"\n")
        return
    nchar = max(len(s) for s in seqs)
    for block, start in enumerate(range(0, nchar, width)):
        for name, seq in zip(names, seqs):
            prefix = name + b" " if block == 0 or ids_every_block else b""
            f.write(prefix + seq[start : start + width] + b"\n")
        f.write(b"\n")


# Generators -------------------------------------------------------------


def gen_loci(out: Path, rng: random.Random, p: dict) -> None:
    """Many locus alignments with random taxon subsets, as in phylogenomic datasets."""
    fmt = p.get("format", "nexus")
    names = taxon_names(p["taxa"])
    pool = make_pool(rng, p["max_len"] * 4)
    for locus in range(p["loci"]):
        length = rng.randint(p["min_len"], p["max_len"])
        start = rng.randrange(len(pool) - length)
        base = pool[start : start + length]
        kept = [n for n in names if rng.random() >= p.get("missing", 0.0)] or names[:1]
        seqs = [mutate(rng, base, p.get("variation", 0.01), p.get("gaps", 0.0)) for _ in kept]
        path = out / f"locus_{locus:05d}.{EXTENSIONS[fmt]}"
        write_alignment(path, fmt, kept, seqs, p.get("layout", "sequential"), p.get("width", 80))


def gen_alignment(out: Path, rng: random.Random, p: dict) -> None:
    """A single large alignment."""
    fmt = p.get("format", "nexus")
    names = taxon_names(p["taxa"])
    base = make_pool(rng, p["nchar"])
    seqs = [mutate(rng, base, p.get("variation", 0.01), p.get("gaps", 0.0)) for _ in names]
    path = out / f"alignment.{EXTENSIONS[fmt]}"
    write_alignment(path, fmt, names, seqs, p.get("layout", "sequential"), p.get("width", 1000))


def gen_contigs(out: Path, rng: random.Random, p: dict) -> None:
    """Assembly contigs of varied length in wrapped FASTA, with some Ns."""
    pool = make_pool(rng, p["max_len"] * 4)
    width = p.get("width", 80)
    n_rate = p.get("n_rate", 0.001)
    for file_idx in range(p.get("files", 1)):
        with open(out / f"contigs_{file_idx:03d}.fa", "wb") as f:
            for contig in range(p["contigs"]):
                length = rng.randint(p["min_len"], p["max_len"])
                start = rng.randrange(len(pool) - length)
                seq = bytearray(pool[start : start + length])
                for pos in rng.sample(range(length), int(length * n_rate)):
                    seq[pos] = ord("N")
                f.write(f">contig_{contig}\n".encode())
                for s in range(0, length, width):
                    f.write(seq[s : s + width] + b"\n")


GENERATORS = {"loci": gen_loci, "alignment": gen_alignment, "contigs": gen_contigs}
