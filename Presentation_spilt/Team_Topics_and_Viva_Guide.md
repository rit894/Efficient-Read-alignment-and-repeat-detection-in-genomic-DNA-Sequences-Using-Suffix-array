# Team topics and viva guide

Project: Efficient Read Alignment and Repeat Detection in Genomic DNA Using Suffix Arrays

Use Presenter 1–4 until you add your names. The deck has 20 slides, five per person. Aim for about 3–4 minutes each plus a shared demo. This allocates presentation and study work; it does not claim who originally wrote the code.

## The five ADTs actually defined

| Abstract class | Implementation | Operations in the abstract contract |
|---|---|---|
| SuffixArrayADT | SuffixArray | build, get, bounds, search, count, __len__ |
| LCPArrayADT | LCPArray | build, get, max_lcp, lcp_between, __len__ |
| BWTADT | BWT | build, get, inverse, __len__ |
| FMIndexADT | FMIndex | build, occ, backward_search, count, locate |
| PatternMatcherADT | NaiveMatcher, KMPMatcher, SuffixArrayMatcher, FMIndexMatcher | preprocess, find_all, count |

The abstract classes live in src/adt.py. Concrete implementations live in src/data_structures.py. Other concrete methods include KMPMatcher.failure_function(), FMIndex.stored_sa_entries() and FMIndex._row_to_position(). These are additional methods, not additional ADTs.

## Exact topic inventory

### Data structures and representations

- Suffix array: suffix start positions in lexicographic order.
- LCP array: longest common prefixes of adjacent suffixes.
- Rank arrays and inverse-rank arrays: construction and Kasai traversal.
- BWT string: preceding characters in suffix order.
- FM-index: BWT, C table, checkpointed Occ counts and sampled suffix-array rows.
- Python lists: suffix positions, ranks, LCP, checkpoints, failure tables and results.
- Dictionaries: C table, checkpoints, SA samples, read records and alignment results.
- Sets: alphabet extraction and duplicate candidate removal.
- Strings: DNA, patterns, BWT and temporary prefix slices.
- Tuples: bounds and results such as (position, strand, mismatches).

### Algorithms and DSA concepts

- Naive string matching.
- Knuth–Morris–Pratt (KMP) and the failure function.
- Prefix doubling with packed integer rank-pair keys.
- Comparison sorting through Python list.sort.
- Binary search with bisect_left and bisect_right.
- Kasai's linear-time LCP construction.
- Linear minimum scanning for lcp_between. No RMQ preprocessing exists here.
- Maximum LCP for longest repeated substring.
- Joined suffix array with a unique separator for longest common substring.
- Burrows–Wheeler Transform and inverse BWT.
- LF mapping and FM-index backward search.
- Occ checkpointing and suffix-array sampling.
- Pigeonhole principle and exact seeds for approximate matching.
- Hamming-distance verification of candidate matches.
- Reverse complement and both-strand read alignment.
- Time/space analysis, amortized reasoning for Kasai, and output-sensitive reporting costs.
- ADTs, abstract classes, inheritance and a common matcher interface.

### Supporting project topics

- Fixed random seeds, GC-controlled DNA generation, planted repeats and mutations.
- Known-origin simulated reads and generated query patterns.
- FASTA reading/writing and input validation.
- Unit tests, fixed cases, invariants and randomized cross-checks.
- Timing with perf_counter and peak Python memory with tracemalloc.
- CSV output, JSON metadata, optional matplotlib plots and argparse CLI.
- GC-content calculation.

KMP and Kasai are algorithms, not extra ADTs. BWT is a transform wrapped in an ADT in this code. The project does not implement suffix trees, tries, linked lists, stacks, queues, heaps, segment trees, SA-IS or radix-sort suffix construction. Discussing an alternative does not mean implementing it.

## Four-person allocation

| Person | Slides | Main topics | Extra preparation to balance the work |
|---|---|---|---|
| Presenter 1 | 1–5 | Problem, five ADTs, architecture, Naive and KMP | Genotype driver, FASTA validation, reproducibility and test harness |
| Presenter 2 | 6–10 | Suffix-array representation, prefix doubling, binary search, justification and complexity | Trace banana and explain array versus linked-list access |
| Presenter 3 | 11–15 | LCP, Kasai, repeats, common substring, approximate matching and read alignment | Trace seed offsets, Hamming distance and reverse complement |
| Presenter 4 | 16–20 | BWT, inverse transform, FM-index, checkpoints, sampling, results and live demo | Explain benchmark limits and coordinate the shared demonstration |

Presenter 1 has lighter opening slides, so preparation of generated inputs and the tests balances their workload. Presenter 4 starts the live demonstration; each person explains the output relevant to their section. Everyone should learn the overall pipeline because faculty may question each person individually.

## Presenter 1: study and viva questions

Study src/adt.py, NaiveMatcher, KMPMatcher, src/genotype_driver.py and tests/test_all.py.

1. What is an ADT? A contract for operations and behavior independent of representation.
2. Why one matcher interface? The same tests and benchmark can call each algorithm consistently.
3. Why compare against KMP? It is a linear-scan baseline with no genome index build.
4. Why fixed seeds? Data is reproducible, although timing can vary.
5. What does the corrected FASTA reader support? One nonempty ACGT sequence, lowercase normalization, and clear rejection of ambiguous symbols or multiple records.
6. What do the tests check? Fixed outputs, suffix order, brute-force LCP, inverse BWT, randomized matcher agreement and read mapping.

Handoff: The baselines rescan the reference. Next we explain a suffix-array index that we reuse across queries.

## Presenter 2: study and viva questions

Study terminate(), _prefix_doubling(), SuffixArray and SuffixArrayMatcher.

1. What does SA store? All suffix starting positions, ordered by suffix text.
2. Why store positions? Explicitly storing every full suffix would require quadratic total characters.
3. Why a list rather than a linked list? Binary search needs efficient indexed row access.
4. Why O(n log² n) construction? At most O(log n) rounds, each comparison-sorting n keys in O(n log n).
5. Why two searches? They locate the first and last boundary of the matching suffix interval.
6. Why occ log occ? The implementation sorts returned positions into text order.

Example: banana$ has SA [6,5,3,1,0,4,2]. For ana, bounds are [2,4). The stored positions are [3,1], and search returns [1,3].

Handoff: Adjacent suffixes reveal common prefixes. The LCP array uses these for repeats and sequence comparison.

## Presenter 3: study and viva questions

Study LCPArray and src/applications.py.

1. Why maximum LCP for repeats? Suffixes sharing a prefix form a consecutive group, so a longest shared prefix appears between adjacent suffixes.
2. Why is Kasai linear? It reuses a prefix length that drops at most once per text step, bounding total extensions.
3. Is lcp_between constant time? No. It scans a slice in O(j-i).
4. Why the # separator? It keeps the two genomes distinct during comparison.
5. Why k+1 seeds? At most k substitutions cannot affect every one of k+1 nonempty pieces.
6. Why use a set? Several seed hits can imply the same full-pattern start.
7. What is reverse complement? Complement A/T and C/G, then reverse. AACG becomes CGTT.
8. Are insertion/deletion errors supported? No. Hamming matching checks substitutions without gaps.

Example: ACGTAC with k=1 gives ACG and TAC. A TAC hit at position 13 implies full-pattern start 10 because its seed offset is 3.

Handoff: These applications reuse an indexed matcher. Next we explain the FM-index alternative and its storage trade-offs.

## Presenter 4: study and viva questions

Study BWT, FMIndex, FMIndexMatcher, src/benchmark.py and main/main.py.

1. Does BWT alone compress data? It rearranges characters. This code stores the result as a string.
2. What is C[c]? The number of text characters smaller than c.
3. What is Occ(c,i)? The number of c symbols in BWT[0:i].
4. Why backward search? C and Occ extend a matching interval one character to the left.
5. Why sampled SA? Fewer stored positions, with LF walking to recover others.
6. What happens when B or s increases? Larger B reduces checkpoints but increases local counting work. Larger s reduces SA samples but increases LF walking.
7. Does 6.25% of SA entries mean 6.25% total memory? No. BWT, dictionaries, checkpoints and construction arrays also occupy memory.
8. What do query timings exclude? Preprocessing. The recorded FM build column also assumes an existing SA.

Example: banana$ has BWT annb$aa. Backward search for ana uses [0,7), [1,4), [5,7), [2,4), leaving two matches.

## Complexity reference

n = indexed text length, m = pattern length, occ = reported matches, B = Occ checkpoint spacing, s = SA sample interval.

| Operation | Time |
|---|---|
| SA construction | O(n log² n) |
| LCP given SA | O(n) |
| BWT given SA | O(n) |
| Inverse BWT in this code | O(n log n) |
| SA bounds/count | O(m log n) |
| SA search with sorted output | O(m log n + occ log occ) |
| FM construction given SA | O(sigma n), linear for the fixed DNA alphabet |
| FM count | O(mB), or O(m) for fixed B |
| FM sorted locate | O(mB + occ s B + occ log occ) |
| KMP | O(n+m) |
| Naive | O(nm) worst case |
| Longest repeat | O(n log² n) total |
| Longest common substring | O((a+b) log²(a+b)) |

SA construction uses O(n) auxiliary words. SA search uses O(m+occ) including temporary slices and output. KMP uses O(m+occ) including its failure table and output. Naive also uses temporary slices and output. FM locate uses O(occ) result/sort storage beyond its index. Theoretical packed bit counts do not describe Python object overhead.

Approximate matching: if m <= k, the implementation verifies all starts in O((n-m+1)m). Otherwise, with suffix-array seeds, let H_j be the hit count and l_j the length of seed j, and let C be the number of unique valid candidate starts. Seed searches cost O(sum_j(l_j log n + H_j log H_j)). Candidate sorting costs O(C log C), and verification costs O(Cm). The sum of seed lengths is m. Set insertion uses expected constant time. Each read performs two such searches, reusing one reference index. Repetitive text can produce many seed hits and candidates.

## Live demo

Run from the project root with Python 3.10 or newer:

```bash
python -m main.main test
python -m main.main demo
python -m main.main search --fasta data/sample_genome.fasta --pattern GAATTC --method sa
python -m main.main repeats --fasta data/sample_genome.fasta
python -m main.main align --fasta data/sample_genome.fasta --reads 50 --mismatch 2
```

The reviewed copy passed the 17 supplied tests, including after the genotype_driver update. Rehearse the current team copy after all edits. Slide 19 uses the supplied recorded 200,000-base benchmark, not a new full run.

A simulated read with more errors than the mismatch budget may legitimately remain unmapped. Overlapping repeat occurrences count. Reported positions start from zero.

## Individual GitHub proof

Each member should show their own account and genuine commits, explain what they changed, and understand the overall pipeline. This presentation split does not establish historical contribution. Rehearse by having each member explain their section, trace one example, run one command and answer questions from a teammate.
