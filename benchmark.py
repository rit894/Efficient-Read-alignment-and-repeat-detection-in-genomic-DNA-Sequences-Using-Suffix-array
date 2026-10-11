"""
benchmark.py - Time and memory comparison
=========================================
Compares Naive, KMP, Suffix Array and FM-Index search on synthetic genomes of
growing size, checks that all four give identical answers, and writes
    results/benchmark.csv      the raw numbers
    results/benchmark.png      query-time and build-time plots (if matplotlib exists)
"""

import csv
import os
import time
import tracemalloc
from typing import Dict, List

from .data_structures import ALL_MATCHERS, FMIndex, LCPArray, SuffixArray
from .genotype_driver import pattern_queries, plant_repeats, random_genome

DEFAULT_SIZES = [1_000, 5_000, 10_000, 50_000, 100_000, 200_000]


def _peak_memory_mb(fn, *args) -> float:
    """Peak extra memory (MB) used while running fn(*args)."""
    tracemalloc.start()
    fn(*args)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak / 1e6


def run_benchmark(sizes: List[int] = None, n_queries: int = 100, pattern_len: int = 20,
                  seed: int = 1, verbose: bool = True) -> List[Dict]:
    """Return one result row per genome size."""
    sizes = sizes or DEFAULT_SIZES
    rows = []
    if verbose:
        print(f"{'n':>8} | {'SA build':>9} {'LCP':>7} {'FM build':>9} | "
              + " ".join(f"{c.name + ' ms/q':>14}" for c in ALL_MATCHERS))
    for n in sizes:
        genome = plant_repeats(random_genome(n, seed=seed), motif_len=min(200, n // 10),
                               copies=5, seed=seed + 1)
        queries = pattern_queries(genome, n_queries, pattern_len, hit_ratio=0.8, seed=seed + 2)

        t = time.perf_counter()
        sa = SuffixArray(genome)
        t_sa = time.perf_counter() - t
        t = time.perf_counter()
        LCPArray(genome, sa.sa)
        t_lcp = time.perf_counter() - t
        t = time.perf_counter()
        fm = FMIndex(genome, sa.sa)
        t_fm = time.perf_counter() - t

        row = {
            "n": n,
            "sa_build_s": round(t_sa, 4),
            "lcp_build_s": round(t_lcp, 4),
            "fm_build_s": round(t_fm, 4),
            "sa_peak_mem_mb": round(_peak_memory_mb(SuffixArray, genome), 2),
            "fm_stored_sa_entries": fm.stored_sa_entries(),
            "sa_entries": len(sa),
        }

        reference = None
        for cls in ALL_MATCHERS:
            m = cls()
            m.preprocess(genome)
            t = time.perf_counter()
            results = [m.find_all(q) for q in queries]
            per_query_ms = (time.perf_counter() - t) / len(queries) * 1000
            if reference is None:
                reference = results
            assert results == reference, f"{m.name} disagrees with Naive at n={n}"
            row[f"{m.name}_ms_per_query".replace(" ", "_").replace("-", "_")] = round(per_query_ms, 5)
        rows.append(row)

        if verbose:
            print(f"{n:>8,} | {t_sa:>8.3f}s {t_lcp:>6.3f}s {t_fm:>8.3f}s | "
                  + " ".join(f"{row[c.name.replace(' ', '_').replace('-', '_') + '_ms_per_query']:>14.5f}"
                             for c in ALL_MATCHERS))
    if verbose:
        print("All four methods returned identical results at every size.")
    return rows


def save_csv(rows: List[Dict], path: str = "results/benchmark.csv") -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def save_plot(rows: List[Dict], path: str = "results/benchmark.png") -> bool:
    """Save query-time and build-time plots. Returns False if matplotlib is missing."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    ns = [r["n"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

    for key, label in [("Naive", "Naive"), ("KMP", "KMP"),
                       ("Suffix_Array", "Suffix Array"), ("FM_Index", "FM-Index")]:
        ax[0].plot(ns, [r[f"{key}_ms_per_query"] for r in rows], marker="o", label=label)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("Genome length n (bases)")
    ax[0].set_ylabel("Time per query (ms)")
    ax[0].set_title("Query time (pattern length 20)")
    ax[0].grid(True, which="both", alpha=0.3)
    ax[0].legend()

    ax[1].plot(ns, [r["sa_build_s"] for r in rows], marker="o", label="Suffix array")
    ax[1].plot(ns, [r["lcp_build_s"] for r in rows], marker="o", label="LCP array")
    ax[1].plot(ns, [r["fm_build_s"] for r in rows], marker="o", label="FM-index (given SA)")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Genome length n (bases)")
    ax[1].set_ylabel("Build time (s)")
    ax[1].set_title("Index construction time")
    ax[1].grid(True, which="both", alpha=0.3)
    ax[1].legend()

    fig.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return True


if __name__ == "__main__":
    data = run_benchmark()
    save_csv(data)
    print("Plot saved." if save_plot(data) else "matplotlib not installed - plot skipped.")
