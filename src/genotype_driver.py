"""
genotype_driver.py - Random DNA (genotype) data generator
=========================================================
Generates reproducible synthetic data for testing and benchmarking:

    random_genome        random DNA with a chosen GC content
    plant_repeats        copy a motif to several places (so repeats exist)
    mutate               introduce substitutions (SNP-like differences)
    simulate_reads       short reads with sequencing errors, either strand
    pattern_queries      query patterns with a controlled share of hits
    generate_test_suite  one call that returns everything as a dict (JSON-able)

Command line:
    python -m src.genotype_driver --length 5000 --reads 50 --out results/generated
"""

import argparse
import json
import os
import random
from typing import Dict, List, Optional

from .applications import reverse_complement

BASES = "ACGT"


def random_genome(length: int, gc_content: float = 0.5, seed: Optional[int] = None) -> str:
    """Random DNA of `length` bases where P(G) = P(C) = gc/2 and P(A) = P(T) = (1-gc)/2."""
    if length < 0:
        raise ValueError("genome length must be zero or greater")
    if not 0 <= gc_content <= 1:
        raise ValueError("GC content must be between 0 and 1")
    rnd = random.Random(seed)
    weights = [(1 - gc_content) / 2, gc_content / 2, gc_content / 2, (1 - gc_content) / 2]
    return "".join(rnd.choices(BASES, weights=weights, k=length))


def plant_repeats(genome: str, motif_len: int = 200, copies: int = 5,
                  seed: Optional[int] = None) -> str:
    """Overwrite `copies` random places with the same random motif (creates exact repeats)."""
    rnd = random.Random(seed)
    if motif_len * 2 > len(genome):
        return genome
    motif = "".join(rnd.choice(BASES) for _ in range(motif_len))
    g = list(genome)
    for _ in range(copies):
        pos = rnd.randrange(0, len(g) - motif_len)
        g[pos:pos + motif_len] = motif
    return "".join(g)


def mutate(seq: str, rate: float = 0.01, seed: Optional[int] = None) -> str:
    """Substitute each base with probability `rate` by a different base (SNP simulation)."""
    if not 0 <= rate <= 1:
        raise ValueError("mutation rate must be between 0 and 1")
    rnd = random.Random(seed)
    out = []
    for c in seq:
        if rnd.random() < rate:
            out.append(rnd.choice([b for b in BASES if b != c]))
        else:
            out.append(c)
    return "".join(out)


def simulate_reads(genome: str, n_reads: int = 20, read_len: int = 50,
                   error_rate: float = 0.0, both_strands: bool = True,
                   seed: Optional[int] = None) -> List[Dict]:
    """
    Simulate short sequencing reads.
    Each read: {"id", "seq", "true_pos", "strand", "mismatches"} where true_pos is the
    start in the forward genome and strand '-' means the read is the reverse complement.
    """
    if n_reads < 0:
        raise ValueError("number of reads must be zero or greater")
    if read_len <= 0:
        raise ValueError("read length must be positive")
    if n_reads > 0 and read_len > len(genome):
        raise ValueError("read length cannot exceed genome length")
    if not 0 <= error_rate <= 1:
        raise ValueError("error rate must be between 0 and 1")

    rnd = random.Random(seed)
    reads = []
    for i in range(n_reads):
        pos = rnd.randrange(0, len(genome) - read_len + 1)
        frag = genome[pos:pos + read_len]
        noisy = mutate(frag, error_rate, seed=rnd.randrange(1 << 30))
        strand = "-" if both_strands and rnd.random() < 0.5 else "+"
        seq = reverse_complement(noisy) if strand == "-" else noisy
        reads.append({
            "id": f"read_{i:04d}",
            "seq": seq,
            "true_pos": pos,
            "strand": strand,
            "mismatches": sum(x != y for x, y in zip(frag, noisy)),
        })
    return reads


def pattern_queries(genome: str, n: int = 50, length: int = 20, hit_ratio: float = 0.8,
                    seed: Optional[int] = None) -> List[str]:
    """
    `n` query patterns. About hit_ratio of them are real substrings of the genome
    (guaranteed hits); the rest are random strings (almost surely misses).
    """
    if n < 0:
        raise ValueError("number of queries must be zero or greater")
    if length <= 0:
        raise ValueError("pattern length must be positive")
    if not 0 <= hit_ratio <= 1:
        raise ValueError("hit ratio must be between 0 and 1")
    if n > 0 and hit_ratio > 0 and length > len(genome):
        raise ValueError("pattern length cannot exceed genome length when generating hits")
    rnd = random.Random(seed)
    queries = []
    for _ in range(n):
        if rnd.random() < hit_ratio:
            p = rnd.randrange(0, len(genome) - length + 1)
            queries.append(genome[p:p + length])
        else:
            queries.append("".join(rnd.choice(BASES) for _ in range(length)))
    return queries


def write_fasta(path: str, name: str, seq: str, width: int = 70) -> None:
    """Write one sequence as a FASTA file."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def read_fasta(path: str) -> str:
    """Read one non-empty ACGT sequence; reject unsupported symbols or records.

    Lowercase bases are converted to uppercase. A single FASTA header is
    optional. Ambiguous bases such as N are rejected rather than removed,
    preserving sequence positions.
    """
    parts = []
    headers = 0
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                headers += 1
                if headers > 1 or parts:
                    raise ValueError("Please provide one FASTA sequence per file")
                continue
            sequence = line.upper()
            invalid = set(sequence) - set(BASES)
            if invalid:
                raise ValueError(
                    "Unsupported DNA symbols: " + ", ".join(sorted(invalid))
                )
            parts.append(sequence)
    genome = "".join(parts)
    if not genome:
        raise ValueError("The file contains no DNA sequence")
    return genome


def generate_test_suite(length: int = 5000, gc_content: float = 0.5, n_reads: int = 30,
                        read_len: int = 50, error_rate: float = 0.02, n_queries: int = 50,
                        pattern_len: int = 15, seed: int = 7) -> Dict:
    """Build a complete random test suite (genome, a mutated strain, reads, queries)."""
    if length <= 0:
        raise ValueError("genome length must be positive")
    if not 0 <= gc_content <= 1:
        raise ValueError("GC content must be between 0 and 1")
    genome = plant_repeats(random_genome(length, gc_content, seed), seed=seed + 1)
    strain = mutate(genome, rate=0.01, seed=seed + 2)
    return {
        "seed": seed,
        "genome_length": len(genome),
        "gc_content": round((genome.count("G") + genome.count("C")) / len(genome), 4),
        "genome": genome,
        "strain": strain,
        "reads": simulate_reads(genome, n_reads, read_len, error_rate, seed=seed + 3),
        "queries": pattern_queries(genome, n_queries, pattern_len, seed=seed + 4),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate random DNA test data")
    ap.add_argument("--length", type=int, default=5000, help="genome length in bases")
    ap.add_argument("--gc", type=float, default=0.5, help="GC content between 0 and 1")
    ap.add_argument("--reads", type=int, default=30, help="number of simulated reads")
    ap.add_argument("--read-len", type=int, default=50)
    ap.add_argument("--error", type=float, default=0.02, help="per-base read error rate")
    ap.add_argument("--queries", type=int, default=50)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="results/generated", help="output path prefix")
    a = ap.parse_args()

    try:
        suite = generate_test_suite(
            a.length, a.gc, a.reads, a.read_len, a.error, a.queries, seed=a.seed
        )
    except ValueError as exc:
        ap.error(str(exc))
    write_fasta(a.out + "_genome.fasta", "synthetic_genome", suite["genome"])
    write_fasta(a.out + "_strain.fasta", "synthetic_strain", suite["strain"])
    meta = {k: v for k, v in suite.items() if k not in ("genome", "strain")}
    with open(a.out + "_reads_queries.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote {a.out}_genome.fasta, {a.out}_strain.fasta, {a.out}_reads_queries.json")
    print(f"Genome length {suite['genome_length']}, GC {suite['gc_content']}, "
          f"{len(suite['reads'])} reads, {len(suite['queries'])} queries")


if __name__ == "__main__":
    main()
