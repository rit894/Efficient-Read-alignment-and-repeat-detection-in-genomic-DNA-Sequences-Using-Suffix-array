# DNA Sequence Search and Pattern Matching with Suffix Arrays

**DSA Case Study:** Efficient read alignment and repeat detection in genomic DNA
using suffix arrays, LCP arrays and the FM-index (Python).

Team: *(fill in names / roll numbers / section)*

---

## What this project does

| Data structure | Used for |
|---|---|
| **Suffix Array** | exact pattern search in O(m log n) after one-time indexing |
| **LCP Array** | longest repeated substring, comparing two genomes |
| **BWT + FM-Index** | memory-efficient search (backward search, sampled suffix array) |
| Naive, KMP | baselines for the comparison |

Applications: exact search, approximate search (<= k mismatches), read alignment on
both DNA strands, longest repeated substring, longest common substring of two genomes.

## Project layout

```
dna_suffix_array_project/
|-- README.md                  <- this file
|-- main.py                    <- entry point with main() and command-line interface
|-- requirements.txt           <- optional: matplotlib (only for plots)
|-- src/
|   |-- adt.py                 <- ADT definitions (abstract classes) for every structure
|   |-- data_structures.py     <- SuffixArray, LCPArray, BWT, FMIndex, Naive, KMP
|   |-- applications.py        <- repeats, genome comparison, approximate search, read alignment
|   |-- genotype_driver.py     <- random DNA / reads / queries generator
|   `-- benchmark.py           <- timing and memory comparison, CSV + plot
|-- tests/
|   |-- test_cases.json        <- constant, hand-checked test cases
|   `-- test_all.py            <- 17 unit tests (JSON cases + randomized cross-checks)
|-- docs/
|   |-- PROJECT_INFO.md        <- problem statement, objectives, scope, limitations
|   |-- ADT_SPECIFICATION.md   <- ADT of every data structure (operations, invariants)
|   |-- COMPLEXITY_ANALYSIS.md <- time and space complexity + measured results
|   |-- PRESENTATION_OUTLINE.md<- slide-by-slide outline and demo script
|   |-- GLOSSARY.md            <- terms explained simply
|   `-- REFERENCES.md
|-- data/
|   |-- sample_genome.fasta    <- 20,000-base synthetic genome
|   `-- sample_strain.fasta    <- the same genome with 1 % mutations
|-- results/
|   |-- benchmark.csv / benchmark.png / benchmark_output.txt
|   `-- sample_output.txt      <- output of `python main.py demo`
`-- standalone/
    `-- dna_suffix_array_all_in_one.py   <- everything merged into one file with main()
```
## Workflow in one line

FASTA genome → SuffixArray → LCP / BWT → FMIndex → applications (search, align, repeats, compare) → tests and benchmark
## Quick start

Requires Python 3.10 or newer. No packages are needed (matplotlib is optional, for the plot).

```bash
python -m main.main demo
python -m main.main test
python -m main.main search --fasta data/sample_genome.fasta --pattern GAATTC --method sa
python -m main.main repeats --fasta data/sample_genome.fasta
python -m main.main align --fasta data/sample_genome.fasta --reads 50 --mismatch 2
```

Search `--method` options: `naive`, `kmp`, `sa` (default), `fm`.

## Using the code in your own script

```python
from src.data_structures import SuffixArray, LCPArray, BWT, FMIndex

sa = SuffixArray("ACGTACGTAC")
print(sa.search("ACG"))            # [0, 4]
lcp = LCPArray("ACGTACGTAC", sa.sa)
print(lcp.max_lcp())               # longest repeat length and row
fm = FMIndex("ACGTACGTAC")
print(fm.count("ACG"), fm.locate("ACG"))
```

## Headline results (200,000-base genome, pattern length 20)

| Method | Time per query |
|---|---:|
| Naive | 18.9 ms |
| KMP | 19.7 ms |
| FM-Index | 0.029 ms |
| Suffix Array | 0.013 ms |

Index build (suffix array) takes about 1.5 s once. The sampled FM-index keeps only 6.25 % of
the suffix-array entries. Details in `docs/COMPLEXITY_ANALYSIS.md`.

## Reproducibility

All random data uses fixed seeds, so the same command gives the same genome, reads and
test results every time.
