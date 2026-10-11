# Abstract Data Types (ADTs) Used in the Project

An **ADT** describes *what* a data structure stores and *which operations* it
supports, without saying *how* it is implemented. The Python interfaces are in
`src/adt.py`; the implementations are in `src/data_structures.py`.

**Notation:** `T` = text (DNA string) ending with the unique sentinel `$`,
`n = |T|` (including `$`), `P` = pattern, `m = |P|`, `occ` = number of matches.

Running example: `T = banana$`

---

## 1. Suffix Array (`SuffixArrayADT` -> `SuffixArray`)

**What it is:** the start positions of all suffixes of `T`, listed in the
alphabetical order of those suffixes.

| Row | SA | Suffix |
|---|---|---|
| 0 | 6 | `$` |
| 1 | 5 | `a$` |
| 2 | 3 | `ana$` |
| 3 | 1 | `anana$` |
| 4 | 0 | `banana$` |
| 5 | 4 | `na$` |
| 6 | 2 | `nana$` |

**Invariant:** `SA` is a permutation of `0..n-1` and `T[SA[0]:] < T[SA[1]:] < ...`.

| Operation | Contract |
|---|---|
| `build(text)` | Pre: text has no `$` inside. Post: invariant holds. |
| `get(i)` | Returns `SA[i]`. Pre: `0 <= i < n`. |
| `bounds(P)` | Returns row interval `[lo, hi)` of suffixes starting with `P`. Pre: `m >= 1`. |
| `search(P)` | Returns sorted start positions of every occurrence of `P`. |
| `count(P)` | Returns `hi - lo`. |
| `len()` | Returns `n`. |

**Why it works for search:** all suffixes starting with `P` are *consecutive*
in sorted order, so two binary searches find the whole block.

**Implementation:** Python list of ints, built by prefix doubling.

---

## 2. LCP Array (`LCPArrayADT` -> `LCPArray`)

**What it is:** `LCP[i]` = length of the longest common prefix of the suffixes
in SA rows `i-1` and `i` (`LCP[0] = 0`).

For `banana$`: `LCP = [0, 0, 1, 3, 0, 0, 2]`
(`a$` & `ana$` share `a`; `ana$` & `anana$` share `ana`; `na$` & `nana$` share `na`.)

**Invariant:** `0 <= LCP[i] <= n - max(SA[i-1], SA[i])`.

| Operation | Contract |
|---|---|
| `build(text, sa)` | Pre: `sa` is the suffix array of `text`. |
| `get(i)` | Returns `LCP[i]`. |
| `max_lcp()` | Returns `(value, index)` of the largest entry = longest repeated substring. |
| `lcp_between(i, j)` | LCP of rows `i < j` = `min(LCP[i+1..j])`. |
| `len()` | Returns `n`. |

**Implementation:** Python list, built with Kasai's algorithm.

---

## 3. Burrows-Wheeler Transform (`BWTADT` -> `BWT`)

**What it is:** a string `L` where `L[i]` is the character *just before* the
i-th smallest suffix (`L[i] = T[SA[i]-1]`, wrapping to the last char when `SA[i] = 0`).

For `banana$`: `L = annb$aa`

| Row | Suffix | Preceding char |
|---|---|---|
| 0 | `$` | a |
| 1 | `a$` | n |
| 2 | `ana$` | n |
| 3 | `anana$` | b |
| 4 | `banana$` | $ |
| 5 | `na$` | a |
| 6 | `nana$` | a |

**Invariant:** `L` is a permutation of the characters of `T`; the transform is
reversible.

| Operation | Contract |
|---|---|
| `build(text, sa)` | Builds `L` from the suffix array. |
| `get(i)` | Returns `L[i]`. |
| `inverse()` | Returns the original text (with `$`) using the LF-mapping. |
| `len()` | Returns `n`. |

**Why it matters:** equal characters cluster together (good for compression),
and `L` is the backbone of the FM-index.

---

## 4. FM-Index (`FMIndexADT` -> `FMIndex`)

**What it is:** a compressed index that answers pattern queries using the BWT
instead of the full suffix array.

**Stored state**

| Part | Meaning |
|---|---|
| `L` | the BWT string |
| `C[c]` | number of characters in `T` that are smaller than `c` |
| `Occ(c, i)` | number of `c` in `L[0:i]`; stored as checkpoints every `occ_step` positions, rest counted on the fly |
| sampled SA | only rows whose text position is a multiple of `sa_sample`; other positions are recovered by walking back with LF |

**Invariant:** after processing a suffix of `P`, the interval `[lo, hi)` is
exactly the set of SA rows whose suffixes start with that suffix of `P`.

**Backward search update** (for each character `c` of `P`, right to left):

```
lo = C[c] + Occ(c, lo)
hi = C[c] + Occ(c, hi)       # stop with "no match" if lo >= hi
```

| Operation | Contract |
|---|---|
| `build(text, sa=None)` | Builds `L`, `C`, checkpoints and the sampled SA. |
| `occ(c, i)` | Count of `c` in `L[0:i]`. |
| `backward_search(P)` | Returns `[lo, hi)`; empty if `P` does not occur. |
| `count(P)` | `hi - lo`. |
| `locate(P)` | Sorted start positions (uses sampled SA + LF-walk). |

---

## 5. Pattern Matcher (`PatternMatcherADT`)

A common interface so that four very different methods can be tested and timed
in exactly the same way.

| Operation | Contract |
|---|---|
| `preprocess(text)` | Prepare the text (build an index, or just store it). |
| `find_all(P)` | Sorted start positions of every occurrence (overlaps allowed). |
| `count(P)` | Number of occurrences. |

| Implementation | Idea |
|---|---|
| `NaiveMatcher` | Compare `P` at every alignment. |
| `KMPMatcher` | Failure function avoids re-comparing characters. |
| `SuffixArrayMatcher` | Binary search on the suffix array. |
| `FMIndexMatcher` | Backward search on the BWT. |

---

## How the structures depend on each other

```
Text T --build--> Suffix Array --Kasai--> LCP Array
                       |
                       +--(preceding char)--> BWT --+C, Occ, sampled SA--> FM-Index
```
