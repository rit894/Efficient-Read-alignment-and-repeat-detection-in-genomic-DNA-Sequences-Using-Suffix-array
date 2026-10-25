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
|-- main/
|   |-- __init__.py
|   `-- main.py                <- entry point with main() and command-line interface
|-- requirements.txt           <- optional: matplotlib (only for plots)
|-- src/
|   |-- adt.py                 <- ADT definitions (abstract classes) for every structure
|   |-- data_structures.py     <- SuffixArray, LCPArray, BWT, FMIndex, Naive, KMP
|   |-- applications.py        <- repeats, genome comparison, approximate search, read alignment
|   |-- genotype_driver.py     <- random DNA / reads / queries generator
|   `-- benchmark.py           <- timing and memory comparison, CSV + plot
|-- tests/
|   |-- test_cases.json        <- constant, hand-checked test cases
|   |-- test_all.py            <- 17 algorithm tests
|   `-- test_reference_cli.py  <- 10 reference-data and command-line tests
|-- docs/
|   |-- PROJECT_INFO.md        <- problem statement, objectives, scope, limitations
|   |-- ADT_SPECIFICATION.md   <- ADT of every data structure (operations, invariants)
|   |-- COMPLEXITY_ANALYSIS.md <- time and space complexity + measured results
|   |-- PRESENTATION_OUTLINE.md<- slide-by-slide outline and demo script
|   |-- GLOSSARY.md            <- terms explained simply
|   |-- REAL_REFERENCE_DATA.md <- NCBI provenance and interpretation
|   `-- REFERENCES.md
|-- data/
|   |-- phiX174_NC_001422.1.fasta <- 5,386-base real NCBI reference
|   |-- sample_genome.fasta    <- 20,000-base synthetic genome
|   `-- sample_strain.fasta    <- the same genome with 1 % mutations
|-- results/
|   |-- benchmark.csv / benchmark.png / benchmark_output.txt
|   `-- sample_output.txt      <- recorded output of the original synthetic demo
`-- standalone/
    `-- dna_suffix_array_all_in_one.py   <- everything merged into one file with main()
```
## Workflow in one line

FASTA genome → SuffixArray → LCP / BWT → FMIndex → applications (search, align, repeats, compare) → tests and benchmark
## Quick start

Requires Python 3.10 or newer. No packages are needed (matplotlib is optional, for the plot). Run these commands from the project root.

Start the demo using the bundled real phiX174 reference:

```bash
python -m main.main demo --fasta data/phiX174_NC_001422.1.fasta
```

Part 2 asks which DNA base or sequence to search for. Enter A, GC or ACGT, for example. The reference contains 5,386 bases. Searching for ACGT reports 19 occurrences. Positions start at zero.

To provide a query without a prompt:

```bash
python -m main.main demo --fasta data/phiX174_NC_001422.1.fasta --pattern ACGT
```

The reference is real DNA downloaded from NCBI. The comparison copy and reads are simulated with substitution errors. The reference stays unchanged. The program analyzes the stored sequence linearly and does not include end-to-start matches. See [dataset provenance and limitations](docs/REAL_REFERENCE_DATA.md).

Other commands:

```bash
python -m main.main test
python -m main.main search
python -m main.main search --fasta data/phiX174_NC_001422.1.fasta --pattern ACGT --method sa
python -m main.main repeats --fasta data/phiX174_NC_001422.1.fasta
python -m main.main align --fasta data/phiX174_NC_001422.1.fasta --reads 50 --mismatch 2
python -m main.main demo --pattern GC
```

The search command asks for a FASTA path and pattern when you omit them. Without --fasta, demo uses the original synthetic 3,000-base genome. Without --pattern, demo asks what sequence to search for.

GC content measures the fraction of individual bases that are G or C. It is separate from searching for the consecutive sequence GC.

Search --method options: naive, kmp, sa (default), fm.

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
