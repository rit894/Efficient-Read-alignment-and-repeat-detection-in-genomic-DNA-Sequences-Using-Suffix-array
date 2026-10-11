"""
main.py - Entry point of the DNA Suffix Array case study
========================================================
Usage (run from the project folder):

    python main.py                      run the demo (same as 'demo')
    python main.py demo                 step-by-step walkthrough on 'banana' and on random DNA
    python main.py test                 run all unit tests
    python main.py benchmark            time Naive / KMP / Suffix Array / FM-Index
    python main.py generate             write random genome, strain, reads, queries to files
    python main.py search  --fasta G.fasta --pattern ACGTACGT [--method sa|fm|kmp|naive]
    python main.py repeats --fasta G.fasta
    python main.py compare --fasta A.fasta --fasta2 B.fasta
    python main.py align   --fasta G.fasta [--reads 50 --read-len 50 --error 0.02 --mismatch 2]
"""

import argparse
import os
import sys
import time
import unittest

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from src.data_structures import ALL_MATCHERS, BWT, FMIndex, FMIndexMatcher, KMPMatcher, LCPArray, NaiveMatcher, SuffixArray, SuffixArrayMatcher
from src.applications import align_reads, gc_content, longest_common_substring, longest_repeated_substring
from src.genotype_driver import generate_test_suite, plant_repeats, random_genome, read_fasta, simulate_reads, write_fasta
from src.benchmark import run_benchmark, save_csv, save_plot

METHODS = {"naive": NaiveMatcher, "kmp": KMPMatcher, "sa": SuffixArrayMatcher, "fm": FMIndexMatcher}


# ---------------------------------------------------------------------------
# demo
# ---------------------------------------------------------------------------
def demo() -> None:
    print("=" * 66)
    print(" PART 1 - The three structures on the text 'banana'")
    print("=" * 66)
    text = "banana"
    sa = SuffixArray(text)
    lcp = LCPArray(text, sa.sa)
    bwt = BWT(text, sa.sa)
    print(f"{'row':>3} {'SA':>3} {'LCP':>4} {'BWT':>4}   suffix")
    for i, pos in enumerate(sa.sa):
        print(f"{i:>3} {pos:>3} {lcp.get(i):>4} {bwt.get(i):>4}   {sa.text[pos:]}")
    print(f"\nInverse BWT recovers the text: {bwt.inverse()!r}")
    print(f"Binary search for 'ana'      : positions {sa.search('ana')}")

    fm = FMIndex(text, sa.sa, occ_step=2, sa_sample=2)
    print("\nFM-index backward search for 'ana' (reads the pattern right to left):")
    lo, hi = 0, fm.n
    print(f"  start            rows [{lo}, {hi})")
    for c in reversed("ana"):
        lo, hi = fm.C[c] + fm.occ(c, lo), fm.C[c] + fm.occ(c, hi)
        print(f"  after adding '{c}'  rows [{lo}, {hi})")
    print(f"  -> {hi - lo} occurrences, located at {fm.locate('ana')}")
    print(f"  stored SA entries: {fm.stored_sa_entries()} of {fm.n} (rest recovered by LF-walking)")

    print("\n" + "=" * 66)
    print(" PART 2 - Random DNA (genotype driver)")
    print("=" * 66)
    suite = generate_test_suite(length=3000, n_reads=20, read_len=40, error_rate=0.02, seed=7)
    g = suite["genome"]
    print(f"Genome length {len(g)}, GC content {gc_content(g):.2%}")
    sub, length, pos = longest_repeated_substring(g)
    print(f"Longest repeated substring: length {length}, occurs at {pos[:6]}{'...' if len(pos) > 6 else ''}")
    s, l, pa, pb = longest_common_substring(g, suite["strain"])
    print(f"Longest common substring genome vs mutated strain: length {l} (at {pa} / {pb})")

    matcher = SuffixArrayMatcher()
    matcher.preprocess(g)
    aln = align_reads(suite["reads"], matcher, g, max_mismatch=2)
    mapped = correct = 0
    for r in suite["reads"]:
        hits = aln[r["id"]]
        mapped += bool(hits)
        correct += any(p == r["true_pos"] and st == r["strand"] for p, st, _ in hits)
    print(f"Read alignment (<=2 mismatches, both strands): {mapped}/{len(suite['reads'])} mapped, "
          f"{correct} at the true position")


# ---------------------------------------------------------------------------
# sub-commands
# ---------------------------------------------------------------------------
def cmd_test() -> int:
    for tests_dir in (os.path.join(ROOT, "tests"), os.path.join(os.path.dirname(ROOT), "tests")):
        if os.path.isdir(tests_dir):
            suite = unittest.defaultTestLoader.discover(tests_dir)
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            return 0 if result.wasSuccessful() else 1
    print("tests/ folder not found next to this file.")
    return 1


def cmd_benchmark(args) -> None:
    sizes = [int(s) for s in args.sizes.split(",")] if args.sizes else None
    rows = run_benchmark(sizes, n_queries=args.queries, pattern_len=args.pattern_len)
    out = os.path.join(ROOT, "results")
    save_csv(rows, os.path.join(out, "benchmark.csv"))
    print("Saved results/benchmark.csv", "and results/benchmark.png"
          if save_plot(rows, os.path.join(out, "benchmark.png")) else "(matplotlib missing: no plot)")


def cmd_generate(args) -> None:
    suite = generate_test_suite(args.length, 0.5, args.reads, args.read_len, args.error, 50, seed=args.seed)
    prefix = os.path.join(ROOT, "results", "generated")
    write_fasta(prefix + "_genome.fasta", "synthetic_genome", suite["genome"])
    write_fasta(prefix + "_strain.fasta", "synthetic_strain", suite["strain"])
    print(f"Wrote {prefix}_genome.fasta and {prefix}_strain.fasta ({suite['genome_length']} bases each)")


def cmd_search(args) -> None:
    genome = read_fasta(args.fasta)
    matcher = METHODS[args.method]()
    t = time.perf_counter()
    matcher.preprocess(genome)
    t_pre = time.perf_counter() - t
    t = time.perf_counter()
    hits = matcher.find_all(args.pattern.upper())
    t_q = time.perf_counter() - t
    print(f"{matcher.name}: genome {len(genome):,} bases, preprocess {t_pre:.3f}s, query {t_q * 1000:.3f} ms")
    print(f"{len(hits)} occurrence(s) of {args.pattern.upper()}: {hits[:20]}{' ...' if len(hits) > 20 else ''}")


def cmd_repeats(args) -> None:
    genome = read_fasta(args.fasta)
    sub, length, pos = longest_repeated_substring(genome)
    print(f"Genome length {len(genome):,}. Longest repeated substring: {length} bases")
    print(f"Occurs {len(pos)} times at positions {pos[:10]}")
    print(f"Sequence: {sub[:80]}{'...' if length > 80 else ''}")


def cmd_compare(args) -> None:
    a, b = read_fasta(args.fasta), read_fasta(args.fasta2)
    sub, length, pa, pb = longest_common_substring(a, b)
    print(f"Longest common substring: {length} bases (position {pa} in A, {pb} in B)")
    print(f"Sequence: {sub[:80]}{'...' if length > 80 else ''}")


def cmd_align(args) -> None:
    genome = read_fasta(args.fasta)
    reads = simulate_reads(genome, args.reads, args.read_len, args.error, seed=args.seed)
    matcher = SuffixArrayMatcher()
    matcher.preprocess(genome)
    t = time.perf_counter()
    aln = align_reads(reads, matcher, genome, max_mismatch=args.mismatch)
    elapsed = time.perf_counter() - t
    mapped = sum(bool(aln[r["id"]]) for r in reads)
    correct = sum(any(p == r["true_pos"] and s == r["strand"] for p, s, _ in aln[r["id"]]) for r in reads)
    print(f"{len(reads)} simulated reads (length {args.read_len}, error {args.error:.0%}) "
          f"aligned in {elapsed:.3f}s with <= {args.mismatch} mismatches")
    print(f"mapped: {mapped}, mapped to true position: {correct}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="DNA sequence search with suffix arrays, LCP and FM-index")
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("demo", help="step-by-step walkthrough")
    sub.add_parser("test", help="run unit tests")

    b = sub.add_parser("benchmark", help="compare all search methods")
    b.add_argument("--sizes", help="comma separated genome sizes, e.g. 1000,10000,100000")
    b.add_argument("--queries", type=int, default=100)
    b.add_argument("--pattern-len", type=int, default=20)

    g = sub.add_parser("generate", help="generate random DNA test data")
    g.add_argument("--length", type=int, default=5000)
    g.add_argument("--reads", type=int, default=30)
    g.add_argument("--read-len", type=int, default=50)
    g.add_argument("--error", type=float, default=0.02)
    g.add_argument("--seed", type=int, default=7)

    s = sub.add_parser("search", help="find a pattern in a FASTA genome")
    s.add_argument("--fasta", required=True)
    s.add_argument("--pattern", required=True)
    s.add_argument("--method", choices=METHODS, default="sa")

    r = sub.add_parser("repeats", help="longest repeated substring of a genome")
    r.add_argument("--fasta", required=True)

    c = sub.add_parser("compare", help="longest common substring of two genomes")
    c.add_argument("--fasta", required=True)
    c.add_argument("--fasta2", required=True)

    a = sub.add_parser("align", help="simulate reads and align them to a genome")
    a.add_argument("--fasta", required=True)
    a.add_argument("--reads", type=int, default=50)
    a.add_argument("--read-len", type=int, default=50)
    a.add_argument("--error", type=float, default=0.02)
    a.add_argument("--mismatch", type=int, default=2)
    a.add_argument("--seed", type=int, default=7)

    args = ap.parse_args(argv)
    if args.cmd in (None, "demo"):
        demo()
    elif args.cmd == "test":
        return cmd_test()
    elif args.cmd == "benchmark":
        cmd_benchmark(args)
    elif args.cmd == "generate":
        cmd_generate(args)
    elif args.cmd == "search":
        cmd_search(args)
    elif args.cmd == "repeats":
        cmd_repeats(args)
    elif args.cmd == "compare":
        cmd_compare(args)
    elif args.cmd == "align":
        cmd_align(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
