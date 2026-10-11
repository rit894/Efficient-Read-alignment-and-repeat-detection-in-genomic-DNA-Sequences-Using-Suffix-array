# Time and Space Complexity Analysis

**Notation:** `n` = text length (bases), `m` = pattern length, `occ` = number of
occurrences, `sigma` = alphabet size (4 for DNA, 5 with the sentinel),
`s` = SA sampling rate, `B` = Occ checkpoint step, `k` = allowed mismatches.

---

## 1. Summary tables

### Construction (done once per genome)

| Structure | Algorithm used | Time | Extra space |
|---|---|---|---|
| Suffix array | Prefix doubling + comparison sort | O(n log^2 n) | O(n) |
| Suffix array (theory) | Prefix doubling + radix sort | O(n log n) | O(n) |
| Suffix array (theory) | SA-IS / DC3 | O(n) | O(n) |
| LCP array | Kasai | O(n) | O(n) |
| BWT | From the suffix array | O(n) | O(n) |
| Inverse BWT | LF-mapping with stable sort | O(n log n) (O(n) with counting sort) | O(n) |
| FM-index (given SA) | C table + Occ checkpoints + sampled SA | O(sigma * n) = O(n) for DNA | see section 3 |

### Query (per pattern)

| Method | Preprocessing | Search / count | Locate all `occ` | Query space |
|---|---|---|---|---|
| Naive | none | O(n * m) worst, ~O(n) typical | same pass | O(1) |
| KMP | O(m) per query | O(n + m) | same pass | O(m) |
| Suffix array | O(n log^2 n) once | O(m log n) | + O(occ log occ) to sort positions | O(1) |
| FM-index | O(n log^2 n) once (SA) | O(m * B) = O(m) | + O(occ * s * B) | O(1) |

### Space of the stored index

| Structure | Theoretical space | What dominates |
|---|---|---|
| Suffix array | n * ceil(log2 n) bits | the n positions |
| LCP array | n * ceil(log2 n) bits | the n values |
| BWT | n * ceil(log2 sigma) bits (about 2-3 bits per base) | the string |
| FM-index | n log sigma + sigma * (n/B) log n + (n/s) log n bits | BWT + small tables |

The FM-index stores only `n/s` of the `n` suffix-array entries (with `s = 16`
that is 6.25 %). This is the memory saving that lets real aligners index the
3-billion-base human genome in a few gigabytes.

---

## 2. Why these bounds hold

### Suffix array construction: O(n log^2 n)
Prefix doubling sorts suffixes by their first `2^j` characters in round `j`.
After `ceil(log2 L)` rounds every suffix is distinguished, where `L` is the
longest repeated substring (worst case `L ~ n`, e.g. `AAAA...A`, giving `log n`
rounds). Each round sorts `n` integer keys with a comparison sort: O(n log n).
Total: **O(n log n) rounds-cost x O(log n) rounds = O(n log^2 n)**.
On random DNA the repeats are short, so only about 4-8 rounds are needed.

### Search on the suffix array: O(m log n)
All suffixes starting with `P` form one contiguous block. Two binary searches
(lower and upper bound) each take `log2 n` steps; each step compares at most `m`
characters. So finding the block is **O(m log n)**. Reading the block and sorting
the positions adds O(occ log occ).

### LCP array (Kasai): O(n)
Process suffixes in text order. If suffix `i` has LCP `h` with its predecessor
in SA order, then suffix `i+1` has LCP at least `h - 1`. The counter `h`
therefore starts each step at the previous value minus 1. It can increase at
most `2n` times in total (it never exceeds `n` and decreases at most `n` times),
so the character comparisons total **O(n)**.

### BWT: O(n)
One pass over the suffix array: `L[i] = T[SA[i] - 1]`.

### FM-index backward search: O(m)
Each of the `m` pattern characters needs two `Occ` lookups. With checkpoints
every `B` positions, one `Occ` is a table read plus counting at most `B`
characters: O(B). Total **O(m * B)**, i.e. O(m) for a constant `B`
(here `B = 64`), independent of `n`.

### FM-index locate: O(occ * s * B)
For each matching row, walk backwards with the LF-mapping until a sampled row is
reached. At most `s - 1` steps are needed (a sampled text position appears
every `s` positions), each costing one `Occ` = O(B).
Trade-off: larger `s` = less memory, slower locate.

### Naive and KMP
Naive tries `n - m + 1` alignments of up to `m` comparisons: O(nm) worst case
(e.g. `T = AAAA...A`, `P = AAA...AC`). KMP uses its failure function to never
move backwards in the text: O(n + m), with O(m) space for the table.

### Approximate matching (<= k substitutions, pigeonhole)
Cut `P` into `k + 1` pieces; any match with at most `k` mismatches contains one
piece exactly. Cost: `(k + 1)` exact searches of length about `m / (k + 1)`
plus O(m) Hamming verification per candidate. On random DNA, short pieces have
few false candidates, so the total is close to `(k + 1)` index queries.

### Applications
| Task | Cost |
|---|---|
| Longest repeated substring | build SA + LCP, then one scan: O(n log^2 n) |
| Longest common substring of two genomes (lengths a, b) | build SA + LCP of the joined text, one scan: O((a + b) log^2 (a + b)) |
| Align R reads of length m (both strands, up to k mismatches) | O(R * 2 * ((k+1) * m log n + candidates * m)) after indexing once |

---

## 3. Measured results (this implementation, Python 3.12)

Synthetic genomes with planted repeats, 100 queries of length 20 (80 % real
substrings, 20 % random). All four methods returned identical results at every
size. Raw data: `results/benchmark.csv`, plot: `results/benchmark.png`.

| n (bases) | SA build | LCP build | FM build | Naive ms/query | KMP ms/query | SA ms/query | FM ms/query |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.002 s | 0.000 s | 0.000 s | 0.086 | 0.082 | 0.0058 | 0.0218 |
| 5,000 | 0.013 s | 0.003 s | 0.001 s | 0.469 | 0.407 | 0.0045 | 0.0200 |
| 10,000 | 0.023 s | 0.004 s | 0.002 s | 0.973 | 0.956 | 0.0055 | 0.0214 |
| 50,000 | 0.226 s | 0.033 s | 0.018 s | 5.566 | 5.363 | 0.0116 | 0.0211 |
| 100,000 | 0.742 s | 0.108 s | 0.055 s | 9.277 | 8.118 | 0.0091 | 0.0265 |
| 200,000 | 1.542 s | 0.194 s | 0.073 s | 18.924 | 19.657 | 0.0128 | 0.0291 |

Stored suffix-array entries at n = 200,000: **200,001** (full suffix array)
versus **12,501** in the sampled FM-index (6.25 %).

### What the numbers show
* **Naive and KMP grow linearly with n** (about 20x slower going from 10k to 200k bases).
  Every query rescans the whole genome.
* **Suffix array and FM-index barely change with n**, matching O(m log n) and O(m):
  after the one-time build, the genome length hardly matters.
* **Index search is far faster than scanning** at 200k bases: about 1,500x for the suffix
  array and about 650x for the FM-index (18.9 ms versus 0.013 ms and 0.029 ms per query).
* The one-time build cost pays for itself after a handful of queries
  (SA build 1.5 s versus about 19 ms saved per query at 200k bases, so around 80 queries).

### Honest notes on interpreting the results
* Python timings reflect interpreter constants, not just the asymptotics. For example,
  `str.count` and string slicing run at C speed, so the naive scan is faster than a
  pure-Python loop would be.
* In this test the plain suffix array is faster per query than the FM-index. The FM-index
  pays extra work per query (checkpoint counting and LF-walking) in exchange for storing
  only about 6 % of the suffix array and a compact BWT. Its advantage is **memory at
  genome scale**, not speed on a 200 kb text.
* Suffix array peak memory in Python grows linearly (about 34 MB at 200k bases) because
  Python integers and lists are large; the theoretical figure of n * log2 n bits is far
  smaller.
* Suffix array construction was measured with Python's built-in sort; an SA-IS
  implementation in a compiled language would make the full human genome feasible.

---

## 4. Choosing a method (decision guide)

| Situation | Best choice |
|---|---|
| One pattern, one-off search | KMP (no index to build) |
| Many patterns against the same genome | Suffix array |
| Repeats, genome comparison, longest common substring | Suffix array + LCP |
| Huge genome where memory is the limit | FM-index (sampled SA) |
| Many patterns at once on a small text | Suffix array or Aho-Corasick |
