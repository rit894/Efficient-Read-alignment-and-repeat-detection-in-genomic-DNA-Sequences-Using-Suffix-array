"""
main.py - Entry point of the DNA Suffix Array case study

Usage (run from the project folder):

    python -m main.main                      run the demo and ask for a DNA search in Part 2
    python -m main.main demo                 walkthrough, then ask what to search in the 3,000-base genome
    python -m main.main demo --fasta data/phiX174_NC_001422.1.fasta
    python -m main.main test                 run all unit tests
    python -m main.main benchmark            time Naive / KMP / Suffix Array / FM-Index
    python -m main.main generate             write random genome, strain, reads, queries to files
    python -m main.main search              ask for a FASTA file and DNA pattern
    python -m main.main search  --fasta G.fasta --pattern ACGTACGT [--method sa|fm|kmp|naive]
    python -m main.main repeats --fasta G.fasta
    python -m main.main compare --fasta A.fasta --fasta2 B.fasta
    python -m main.main align   --fasta G.fasta [--reads 50 --read-len 50 --error 0.02 --mismatch 2]
"""

# Built-in modules used throughout the project
import argparse
import os
import sys
import time
import unittest

# Get the absolute path of the project root directory
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add project root to Python's module search path
sys.path.insert(0, ROOT)

# Import data structures and search algorithms
from src.data_structures import (
    ALL_MATCHERS,
    BWT,
    FMIndex,
    FMIndexMatcher,
    KMPMatcher,
    LCPArray,
    NaiveMatcher,
    SuffixArray,
    SuffixArrayMatcher
)

# Import genome analysis applications
from src.applications import (
    align_reads,
    gc_content,
    longest_common_substring,
    longest_repeated_substring
)

# Import data generation and FASTA utilities
from src.genotype_driver import (
    generate_test_suite,
    mutate,
    plant_repeats,
    random_genome,
    read_fasta,
    simulate_reads,
    write_fasta
)

# Import benchmarking utilities
from src.benchmark import run_benchmark, save_csv, save_plot

# Dictionary mapping method names to their matcher classes
METHODS = {
    "naive": NaiveMatcher,
    "kmp": KMPMatcher,
    "sa": SuffixArrayMatcher,
    "fm": FMIndexMatcher
}


# ---------------------------------------------------------------------------
# demo
# ---------------------------------------------------------------------------
def demo(query_pattern=None, fasta_path=None) -> None:
    """
    Demonstrates how suffix arrays, LCP arrays,
    BWT and FM-index work using 'banana' and a
    reference loaded from FASTA or a generated DNA genome.
    """

    print("=" * 66)
    print(" PART 1 - The three structures on the text 'banana'")
    print("=" * 66)

    # Example string
    text = "banana"

    # Build suffix array
    sa = SuffixArray(text)

    # Build longest common prefix array
    lcp = LCPArray(text, sa.sa)

    # Build Burrows Wheeler Transform
    bwt = BWT(text, sa.sa)

    # Display table header
    print(f"{'row':>3} {'SA':>3} {'LCP':>4} {'BWT':>4}   suffix")

    # Display each suffix array row
    for i, pos in enumerate(sa.sa):
        print(f"{i:>3} {pos:>3} {lcp.get(i):>4} {bwt.get(i):>4}   {sa.text[pos:]}")

    # Verify inverse BWT reconstructs original text
    print(f"\nInverse BWT recovers the text: {bwt.inverse()!r}")

    # Demonstrate suffix array substring search
    print(f"Binary search for 'ana'      : positions {sa.search('ana')}")

    # Construct FM-index
    fm = FMIndex(text, sa.sa, occ_step=2, sa_sample=2)

    print("\nFM-index backward search for 'ana' (reads the pattern right to left):")

    # Initial row range
    lo, hi = 0, fm.n
    print(f"  start            rows [{lo}, {hi})")

    # Backward FM-index search
    for c in reversed("ana"):
        lo, hi = fm.C[c] + fm.occ(c, lo), fm.C[c] + fm.occ(c, hi)
        print(f"  after adding '{c}'  rows [{lo}, {hi})")

    # Display occurrences found
    print(f"  -> {hi - lo} occurrences, located at {fm.locate('ana')}")

    # Show FM-index compression statistics
    print(f"  stored SA entries: {fm.stored_sa_entries()} of {fm.n} (rest recovered by LF-walking)")

    print("\n" + "=" * 66)
    if fasta_path:
        print(" PART 2 - DNA from a reference FASTA file")
    else:
        print(" PART 2 - Random DNA (genotype driver)")
    print("=" * 66)

    if fasta_path:
        g = read_fasta(fasta_path)
        read_length = min(40, len(g))
        suite = {
            "genome": g,
            "strain": mutate(g, rate=0.01, seed=9),
            "reads": simulate_reads(
                g, n_reads=20, read_len=read_length,
                error_rate=0.02, seed=10
            ),
        }
        print(f"Reference file: {fasta_path}")
        print("The reference sequence is loaded unchanged from the FASTA file.")
        print(
            f"Test data: 20 simulated reads of length {read_length}, "
            "with a 2% substitution probability per base."
        )
        print("Comparison copy: simulated substitutions with a 1% probability per base.")
        print("The sequence is analyzed linearly; end-to-start matches are not included.")
    else:
        suite = generate_test_suite(
            length=3000,
            n_reads=20,
            read_len=40,
            error_rate=0.02,
            seed=7
        )
        g = suite["genome"]

    # GC content is a composition statistic, separate from a pattern search.
    print(f"Genome length {len(g)}, GC content (G + C bases) {gc_content(g):.2%}")

    if query_pattern is None:
        while True:
            try:
                query_pattern = input(
                    "Enter a DNA base or sequence to search in this genome "
                    "(for example A, GC or ACGT): "
                ).strip().upper()
            except EOFError as exc:
                raise ValueError(
                    "no search input received; run demo --pattern ACGT "
                    "when using non-interactive input"
                ) from exc
            if query_pattern and set(query_pattern) <= set("ACGT"):
                break
            print("Please enter a non-empty sequence using only A, C, G and T.")
    else:
        query_pattern = query_pattern.strip().upper()
        if not query_pattern or not set(query_pattern) <= set("ACGT"):
            raise ValueError(
                "search pattern must be non-empty and contain only A, C, G and T"
            )

    # Reuse one reference index for the user search and read alignment.
    matcher = SuffixArrayMatcher()
    matcher.preprocess(g)
    user_hits = matcher.find_all(query_pattern)
    print(f"Search sequence: {query_pattern}")
    print(f"Occurrences in this {len(g):,}-base genome: {len(user_hits)}")
    if user_hits:
        print(
            f"Positions (zero-based): {user_hits[:20]}"
            f"{' ... (showing first 20 positions)' if len(user_hits) > 20 else ''}"
        )
    else:
        print("This sequence was not found in the generated genome.")

    # Longest repeated substring
    sub, length, pos = longest_repeated_substring(g)

    print(
        f"Longest repeated substring: length {length}, "
        f"occurs at {pos[:6]}{'...' if len(pos) > 6 else ''}"
    )

    # Longest common substring between genome and mutated strain
    s, l, pa, pb = longest_common_substring(g, suite["strain"])

    print(
        f"Longest common substring reference vs simulated mutated copy: "
        f"length {l} (at {pa} / {pb})"
    )

    # Align reads allowing up to 2 mismatches
    aln = align_reads(suite["reads"], matcher, g, max_mismatch=2)

    mapped = correct = 0

    # Evaluate alignment quality
    for r in suite["reads"]:
        hits = aln[r["id"]]

        mapped += bool(hits)

        correct += any(
            p == r["true_pos"] and st == r["strand"]
            for p, st, _ in hits
        )

    print(
        f"Simulated read alignment (<=2 mismatches, both strands): "
        f"{mapped}/{len(suite['reads'])} mapped, "
        f"{correct} at the true position"
    )


# ---------------------------------------------------------------------------
# sub-commands
# ---------------------------------------------------------------------------
def cmd_test() -> int:
    """
    Discover and run all unit tests.
    """

    for tests_dir in (
        os.path.join(ROOT, "tests"),
        os.path.join(os.path.dirname(ROOT), "tests")
    ):
        if os.path.isdir(tests_dir):
            suite = unittest.defaultTestLoader.discover(tests_dir)
            result = unittest.TextTestRunner(verbosity=2).run(suite)

            return 0 if result.wasSuccessful() else 1

    print("tests/ folder not found next to this file.")
    return 1


def cmd_benchmark(args) -> None:
    """
    Run benchmarking and save results.
    """

    sizes = [int(s) for s in args.sizes.split(",")] if args.sizes else None

    rows = run_benchmark(
        sizes,
        n_queries=args.queries,
        pattern_len=args.pattern_len
    )

    out = os.path.join(ROOT, "results")

    # Save benchmark data
    save_csv(rows, os.path.join(out, "benchmark.csv"))

    print(
        "Saved results/benchmark.csv",
        "and results/benchmark.png"
        if save_plot(rows, os.path.join(out, "benchmark.png"))
        else "(matplotlib missing: no plot)"
    )


def cmd_generate(args) -> None:
    """
    Generate synthetic genome and strain files.
    """

    suite = generate_test_suite(
        args.length,
        0.5,
        args.reads,
        args.read_len,
        args.error,
        50,
        seed=args.seed
    )

    prefix = os.path.join(ROOT, "results", "generated")

    # Write generated sequences in FASTA format
    write_fasta(prefix + "_genome.fasta", "synthetic_genome", suite["genome"])
    write_fasta(prefix + "_strain.fasta", "synthetic_strain", suite["strain"])

    print(
        f"Wrote {prefix}_genome.fasta and "
        f"{prefix}_strain.fasta "
        f"({suite['genome_length']} bases each)"
    )


def cmd_search(args) -> None:
    """Search a FASTA genome, prompting for any missing search inputs."""
    try:
        if args.fasta is None:
            default_fasta = os.path.join(ROOT, "data", "sample_genome.fasta")
            entered = input(
                "FASTA file path [Enter for data/sample_genome.fasta]: "
            ).strip()
            fasta_path = entered or default_fasta
        else:
            fasta_path = args.fasta.strip()

        # Accept a pasted path enclosed in matching quotes.
        if (len(fasta_path) >= 2
                and fasta_path[0] == fasta_path[-1]
                and fasta_path[0] in ("'", '"')):
            fasta_path = fasta_path[1:-1]

        genome = read_fasta(fasta_path)
        print(f"Loaded {len(genome):,} DNA bases from {fasta_path}")

        if args.pattern is None:
            while True:
                pattern = input(
                    "DNA sequence to find (A, C, G, T; for example ACGT): "
                ).strip().upper()
                if pattern and set(pattern) <= set("ACGT"):
                    break
                print("Please enter a non-empty sequence using only A, C, G and T.")
        else:
            pattern = args.pattern.strip().upper()
            if not pattern or not set(pattern) <= set("ACGT"):
                raise ValueError(
                    "search pattern must be non-empty and contain only A, C, G and T"
                )
    except EOFError as exc:
        raise ValueError(
            "search input ended early; provide --fasta and --pattern "
            "when running without interactive input"
        ) from exc

    matcher = METHODS[args.method]()

    t = time.perf_counter()
    matcher.preprocess(genome)
    t_pre = time.perf_counter() - t

    t = time.perf_counter()
    hits = matcher.find_all(pattern)
    t_q = time.perf_counter() - t

    print(
        f"{matcher.name}: genome {len(genome):,} bases, "
        f"preprocess {t_pre:.3f}s, query {t_q * 1000:.3f} ms"
    )
    print(
        f"{len(hits)} occurrence(s) of {pattern}: "
        f"{hits[:20]}{' ... (showing first 20 positions)' if len(hits) > 20 else ''}"
    )
    if hits:
        print("Positions are zero-based: the first DNA base is at position 0.")
    else:
        print("The DNA sequence was not found in this FASTA genome.")


def cmd_repeats(args) -> None:
    """
    Find longest repeated substring in a genome.
    """

    genome = read_fasta(args.fasta)

    sub, length, pos = longest_repeated_substring(genome)

    print(
        f"Genome length {len(genome):,}. "
        f"Longest repeated substring: {length} bases"
    )

    print(f"Occurs {len(pos)} times at positions {pos[:10]}")

    print(
        f"Sequence: "
        f"{sub[:80]}{'...' if length > 80 else ''}"
    )


def cmd_compare(args) -> None:
    """
    Find longest common substring between two genomes.
    """

    a = read_fasta(args.fasta)
    b = read_fasta(args.fasta2)

    sub, length, pa, pb = longest_common_substring(a, b)

    print(
        f"Longest common substring: "
        f"{length} bases "
        f"(position {pa} in A, {pb} in B)"
    )

    print(
        f"Sequence: "
        f"{sub[:80]}{'...' if length > 80 else ''}"
    )


def cmd_align(args) -> None:
    """
    Simulate reads and align them against a genome.
    """

    genome = read_fasta(args.fasta)

    reads = simulate_reads(
        genome,
        args.reads,
        args.read_len,
        args.error,
        seed=args.seed
    )

    matcher = SuffixArrayMatcher()

    matcher.preprocess(genome)

    # Measure alignment runtime
    t = time.perf_counter()

    aln = align_reads(
        reads,
        matcher,
        genome,
        max_mismatch=args.mismatch
    )

    elapsed = time.perf_counter() - t

    # Count mapped reads
    mapped = sum(bool(aln[r["id"]]) for r in reads)

    # Count correctly mapped reads
    correct = sum(
        any(
            p == r["true_pos"] and s == r["strand"]
            for p, s, _ in aln[r["id"]]
        )
        for r in reads
    )

    print(
        f"{len(reads)} simulated reads "
        f"(length {args.read_len}, "
        f"error {args.error:.0%}) "
        f"aligned in {elapsed:.3f}s "
        f"with <= {args.mismatch} mismatches"
    )

    print(f"mapped: {mapped}, mapped to true position: {correct}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    """
    Main program entry point.
    Parses command-line arguments and dispatches commands.
    """

    # Create argument parser
    ap = argparse.ArgumentParser(
        description="DNA sequence search with suffix arrays, LCP and FM-index"
    )

    # Create subcommand container
    sub = ap.add_subparsers(dest="cmd")

    # Demo command
    d = sub.add_parser("demo", help="walkthrough with a DNA search in Part 2")
    d.add_argument(
        "--pattern", help="DNA sequence to search in Part 2; asks you if omitted"
    )
    d.add_argument(
        "--fasta", help="reference FASTA for Part 2; otherwise generate synthetic DNA"
    )

    # Unit test command
    sub.add_parser("test", help="run unit tests")

    # Benchmark command
    b = sub.add_parser("benchmark", help="compare all search methods")
    b.add_argument("--sizes")
    b.add_argument("--queries", type=int, default=100)
    b.add_argument("--pattern-len", type=int, default=20)

    # Data generator command
    g = sub.add_parser("generate", help="generate random DNA test data")
    g.add_argument("--length", type=int, default=5000)
    g.add_argument("--reads", type=int, default=30)
    g.add_argument("--read-len", type=int, default=50)
    g.add_argument("--error", type=float, default=0.02)
    g.add_argument("--seed", type=int, default=7)

    # Search command
    s = sub.add_parser("search", help="find a pattern in a FASTA genome")
    s.add_argument("--fasta", help="FASTA file path; asks you if omitted")
    s.add_argument("--pattern", help="DNA sequence to find; asks you if omitted")
    s.add_argument("--method", choices=METHODS, default="sa")

    # Repeats command
    r = sub.add_parser("repeats", help="longest repeated substring of a genome")
    r.add_argument("--fasta", required=True)

    # Compare command
    c = sub.add_parser("compare", help="longest common substring of two genomes")
    c.add_argument("--fasta", required=True)
    c.add_argument("--fasta2", required=True)

    # Alignment command
    a = sub.add_parser("align", help="simulate reads and align them to a genome")
    a.add_argument("--fasta", required=True)
    a.add_argument("--reads", type=int, default=50)
    a.add_argument("--read-len", type=int, default=50)
    a.add_argument("--error", type=float, default=0.02)
    a.add_argument("--mismatch", type=int, default=2)
    a.add_argument("--seed", type=int, default=7)

    # Parse command-line arguments
    args = ap.parse_args(argv)

    # Execute command
    if args.cmd in (None, "demo"):
        demo(
            query_pattern=getattr(args, "pattern", None),
            fasta_path=getattr(args, "fasta", None)
        )
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


# Run program when executed directly
if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(2)
