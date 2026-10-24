# Project Information

## Title
**Efficient Read Alignment and Repeat Detection in Genomic DNA Using Suffix Arrays,
LCP Arrays and the FM-Index**

| Field | Details |
|---|---|
| Course | Data Structures and Algorithms (DSA) case study |
| Language | Python 3.10 or newer (uses `bisect` with `key=`) |
| Team members | *(fill in names and roll numbers)* |
| Faculty / Section | *(fill in)* |
| Primary data structure | Suffix Array |
| Supporting structures | LCP Array, Burrows-Wheeler Transform (BWT), FM-Index |
| Baselines for comparison | Naive search, Knuth-Morris-Pratt (KMP) |

---

## 1. Problem statement

Modern DNA sequencers produce millions of short fragments ("reads", typically
50 to 300 bases) that must be located inside a reference genome that is millions
to billions of bases long. Scientists also need to find repeated regions
(transposons, tandem repeats, duplicated genes), compare genomes of related
strains, and tolerate sequencing errors and mutations.

Scanning the genome from start to finish for every read (Naive or KMP) costs
O(n) per query, which is far too slow when there are millions of reads and
n is in the billions. Storing the genome in an *index* lets each query cost
only O(m log n) or O(m), where `m` is the short pattern length.

**Task:** build and analyse text-indexing data structures for DNA, show how they
answer exact and approximate pattern-matching queries, solve repeat and
genome-comparison problems, and measure the speed and memory trade-offs against
classic scanning algorithms.

---

## 2. Objectives

1. Implement a **suffix array** and use it for exact pattern search.
2. Implement the **LCP array** (Kasai) and use it to find the longest repeated substring.
3. Implement the **BWT** and the **FM-index** (backward search, sampled suffix array).
4. Define an **ADT** for every structure and implement it behind a common interface.
5. Solve real bioinformatics tasks: approximate matching, read alignment on both
   strands, repeat detection, longest common substring of two genomes.
6. Generate reproducible random DNA data (genotype driver) and run constant,
   hand-checked test cases.
7. Analyse time and space complexity and confirm it with measurements.

---

## 3. Real-world context

| Real tool / use | Connection to this project |
|---|---|
| **Bowtie, BWA** (read aligners) | Use the BWT and FM-index exactly as implemented here |
| **MUMmer** (genome comparison) | Uses suffix arrays/trees to find maximal unique matches between genomes |
| **Repeat finders** (e.g. for transposons) | Longest repeated substrings via LCP |
| **Variant / SNP detection** | Reads with a few mismatches must still be placed on the genome |
| **bzip2** | Uses the BWT for compression |

---

## 4. Scope and features

| Feature | Module |
|---|---|
| Suffix array (prefix doubling), binary-search matching | `src/data_structures.py` |
| LCP array (Kasai) | `src/data_structures.py` |
| BWT with inverse transform | `src/data_structures.py` |
| FM-index with checkpointed Occ and sampled SA | `src/data_structures.py` |
| Naive and KMP baselines | `src/data_structures.py` |
| Longest repeated substring | `src/applications.py` |
| Longest common substring of two genomes | `src/applications.py` |
| Approximate matching (<= k mismatches, pigeonhole) | `src/applications.py` |
| Read alignment, forward and reverse-complement strand | `src/applications.py` |
| Random genome, repeats, mutations, reads, queries | `src/genotype_driver.py` |
| Benchmark with CSV and plots | `src/benchmark.py` |
| 17 unit tests + constant JSON cases | `tests/` |

---

## 5. Methodology

1. **Design** the ADTs (operations, invariants) before coding.
2. **Implement** each structure and verify it on the small text `banana$`
   whose answers can be checked by hand.
3. **Cross-validate**: all four search methods must return identical positions
   on constant cases and on randomly generated genomes.
4. **Apply** the structures to the bioinformatics tasks.
5. **Measure** build time, query time and stored entries for genome sizes from
   1,000 to 200,000 bases, and compare with theory.

---

## 6. Key results (summary)

* Index-based search stays almost constant (about 0.01 to 0.03 ms per query)
  while Naive and KMP grow linearly with the genome (about 19 ms per query at 200,000 bases).
* The sampled FM-index stores only 6.25 % of the suffix-array entries.
* Reads simulated with 2 % sequencing error map back to their true position
  (both strands) when up to 2 mismatches are allowed, provided the read has no more
  mismatches than the limit.
* Full numbers: `COMPLEXITY_ANALYSIS.md` and `results/benchmark.csv`.

---

## 7. Limitations and future work

* Python is slow for very large genomes: building the suffix array of the 4.6-million-base
  *E. coli* genome with this code will take many seconds to minutes and a lot of memory.
  A linear-time SA-IS construction in C/C++ or a library would scale to the human genome.
* Approximate matching handles **substitutions only** (Hamming distance), not insertions or deletions.
* The FM-index uses a plain string for the BWT; production tools use wavelet trees
  or bit-vectors with rank support for compression.
* Possible extensions: edit-distance (indel) alignment on candidate regions,
  bidirectional FM-index, generalized suffix array for many genomes, a
  Bloom-filter k-mer baseline, SA-IS construction.

---

## 8. Using real data

Download a genome in FASTA format (for example from NCBI) and run:

```bash
python main.py search  --fasta ecoli.fasta --pattern GAATTC
python main.py repeats --fasta ecoli.fasta
python main.py align   --fasta ecoli.fasta --reads 200 --mismatch 2
```

For a first run on real data, use a slice of a few hundred thousand bases.
