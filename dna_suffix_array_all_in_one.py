"""
dna_suffix_array_all_in_one.py
==============================
Single-file version of the whole case study: ADTs, suffix array, LCP array,
BWT, FM-index, Naive/KMP baselines, bioinformatics applications, random DNA
(genotype) driver, benchmark and a command-line main().

Run:   python dna_suffix_array_all_in_one.py            (demo)
       python dna_suffix_array_all_in_one.py benchmark
       python dna_suffix_array_all_in_one.py --help

(The modular version in src/ plus main.py, tests/ and docs/ is the primary
deliverable; this file is the same code merged for easy submission/demo.)
"""

import argparse
import csv
import json
import os
import random
import sys
import time
import tracemalloc
import unittest
from abc import ABC, abstractmethod
from bisect import bisect_left, bisect_right
from typing import Dict, List, Optional, Tuple

ROOT = os.path.dirname(os.path.abspath(__file__))

# ======================================================================
# SECTION 1 - ADTs
# ======================================================================
from abc import ABC, abstractmethod
from typing import List, Tuple


# ---------------------------------------------------------------------------
# 1. Suffix Array
# ---------------------------------------------------------------------------
class SuffixArrayADT(ABC):
    """
    Suffix Array ADT.

    Abstract state : a permutation SA[0..n-1] of the start positions of all
                     suffixes of T such that T[SA[0]:] < T[SA[1]:] < ... (sorted).
    Invariant      : SA is a permutation of 0..n-1 and suffixes are strictly
                     increasing in lexicographic order.
    """

    @abstractmethod
    def build(self, text: str) -> None:
        """Pre: text is a str without '$' inside. Post: SA satisfies the invariant."""

    @abstractmethod
    def get(self, i: int) -> int:
        """Return SA[i], the start position of the i-th smallest suffix. Pre: 0 <= i < n."""

    @abstractmethod
    def bounds(self, pattern: str) -> Tuple[int, int]:
        """
        Return the half-open row interval [lo, hi) of SA whose suffixes start
        with `pattern`. Empty interval (lo == hi) means no match. Pre: m >= 1.
        """

    @abstractmethod
    def search(self, pattern: str) -> List[int]:
        """Return the sorted list of all start positions of `pattern` in T."""

    @abstractmethod
    def count(self, pattern: str) -> int:
        """Return the number of occurrences of `pattern` in T."""

    @abstractmethod
    def __len__(self) -> int:
        """Return n (number of suffixes, including the sentinel suffix)."""


# ---------------------------------------------------------------------------
# 2. LCP Array
# ---------------------------------------------------------------------------
class LCPArrayADT(ABC):
    """
    LCP Array ADT.

    Abstract state : array LCP[0..n-1] where LCP[i] is the length of the longest
                     common prefix of the suffixes SA[i-1] and SA[i]; LCP[0] = 0.
    Invariant      : 0 <= LCP[i] <= n - max(SA[i-1], SA[i]).
    """

    @abstractmethod
    def build(self, text: str, sa: List[int]) -> None:
        """Pre: sa is the suffix array of text. Post: LCP array is built."""

    @abstractmethod
    def get(self, i: int) -> int:
        """Return LCP[i]. Pre: 0 <= i < n."""

    @abstractmethod
    def max_lcp(self) -> Tuple[int, int]:
        """Return (value, index) of the largest LCP entry (first one if tied)."""

    @abstractmethod
    def lcp_between(self, i: int, j: int) -> int:
        """
        LCP of the suffixes at SA rows i and j (i < j): min(LCP[i+1..j]).
        Pre: 0 <= i < j < n.
        """

    @abstractmethod
    def __len__(self) -> int:
        """Return n."""


# ---------------------------------------------------------------------------
# 3. Burrows-Wheeler Transform
# ---------------------------------------------------------------------------
class BWTADT(ABC):
    """
    BWT ADT.

    Abstract state : string L of length n where L[i] = T[SA[i] - 1]
                     (the character that precedes the i-th smallest suffix,
                     wrapping around to the last character for SA[i] = 0).
    Invariant      : L is a permutation of the characters of T, and the
                     transform is reversible.
    """

    @abstractmethod
    def build(self, text: str, sa: List[int]) -> None:
        """Pre: sa is the suffix array of text. Post: L is built."""

    @abstractmethod
    def get(self, i: int) -> str:
        """Return L[i]. Pre: 0 <= i < n."""

    @abstractmethod
    def inverse(self) -> str:
        """Reconstruct and return the original (sentinel-terminated) text."""

    @abstractmethod
    def __len__(self) -> int:
        """Return n."""


# ---------------------------------------------------------------------------
# 4. FM-Index
# ---------------------------------------------------------------------------
class FMIndexADT(ABC):
    """
    FM-Index ADT (compressed suffix-array-like index).

    Abstract state : BWT string L, the table C[c] (number of characters in T
                     smaller than c), an Occ(c, i) rank structure (number of
                     c's in L[0:i]) and a sampled suffix array.
    Invariant      : after backward search of P, rows [lo, hi) are exactly the
                     SA rows whose suffixes start with P.
    """

    @abstractmethod
    def build(self, text: str, sa: List[int] = None) -> None:
        """Pre: text is a str without '$' inside. Post: index is ready for queries."""

    @abstractmethod
    def occ(self, c: str, i: int) -> int:
        """Return the number of occurrences of character c in L[0:i]."""

    @abstractmethod
    def backward_search(self, pattern: str) -> Tuple[int, int]:
        """Return the SA row interval [lo, hi) for `pattern` (empty if no match)."""

    @abstractmethod
    def count(self, pattern: str) -> int:
        """Return the number of occurrences of `pattern`."""

    @abstractmethod
    def locate(self, pattern: str) -> List[int]:
        """Return the sorted list of start positions of `pattern`."""


# ---------------------------------------------------------------------------
# 5. Common matcher interface (used for the fair benchmark)
# ---------------------------------------------------------------------------
class PatternMatcherADT(ABC):
    """
    Pattern Matcher ADT - the one interface shared by Naive, KMP, Suffix Array
    and FM-Index search so they can be tested and timed identically.
    """

    name: str = "matcher"

    @abstractmethod
    def preprocess(self, text: str) -> None:
        """Prepare the text (build index). Naive and KMP do (almost) nothing here."""

    @abstractmethod
    def find_all(self, pattern: str) -> List[int]:
        """Return the sorted start positions of every occurrence (overlaps allowed)."""

    def count(self, pattern: str) -> int:
        """Number of occurrences. Default: len(find_all)."""
        return len(self.find_all(pattern))


# ======================================================================
# SECTION 2 - Data structures
# ======================================================================
from bisect import bisect_left, bisect_right
from typing import Dict, List, Optional, Tuple


SENTINEL = "$"


def terminate(text: str) -> str:
    """Append the unique sentinel '$' if it is not already the last character."""
    body = text[:-1] if text.endswith(SENTINEL) else text
    if SENTINEL in body:
        raise ValueError("text must not contain the sentinel '$' except as its last character")
    return body + SENTINEL


# ===========================================================================
# 1. Suffix Array
# ===========================================================================
def _prefix_doubling(s: str) -> List[int]:
    """
    Manber-Myers style prefix doubling.
    Round k sorts suffixes by their first 2k characters using the ranks from the
    previous round as (first half, second half) pairs packed into one integer.
    Rounds: O(log n). Each round sorts n integers: O(n log n)  ->  O(n log^2 n).
    """
    n = len(s)
    sa = list(range(n))
    rank = [ord(c) for c in s]
    k = 1
    while True:
        base = max(rank) + 2
        keys = [rank[i] * base + (rank[i + k] + 1 if i + k < n else 0) for i in range(n)]
        sa.sort(key=keys.__getitem__)
        new_rank = [0] * n
        r = 0
        for idx in range(1, n):
            if keys[sa[idx]] != keys[sa[idx - 1]]:
                r += 1
            new_rank[sa[idx]] = r
        rank = new_rank
        if r == n - 1:  # every suffix has a distinct rank -> sorted
            break
        k <<= 1
    return sa


class SuffixArray(SuffixArrayADT):
    """Suffix array stored as a plain Python list of ints."""

    def __init__(self, text: Optional[str] = None):
        self.text: str = ""
        self.sa: List[int] = []
        if text is not None:
            self.build(text)

    def build(self, text: str) -> None:
        self.text = terminate(text)
        self.sa = _prefix_doubling(self.text)

    def get(self, i: int) -> int:
        return self.sa[i]

    def bounds(self, pattern: str) -> Tuple[int, int]:
        if not pattern:
            raise ValueError("pattern must be non-empty")
        m, text = len(pattern), self.text
        key = lambda i: text[i:i + m]  # first m chars of the suffix at SA row
        lo = bisect_left(self.sa, pattern, key=key)
        hi = bisect_right(self.sa, pattern, key=key)
        return lo, hi

    def search(self, pattern: str) -> List[int]:
        lo, hi = self.bounds(pattern)
        return sorted(self.sa[lo:hi])

    def count(self, pattern: str) -> int:
        lo, hi = self.bounds(pattern)
        return hi - lo

    def __len__(self) -> int:
        return len(self.sa)


# ===========================================================================
# 2. LCP Array (Kasai)
# ===========================================================================
class LCPArray(LCPArrayADT):
    """LCP array built with Kasai's linear-time algorithm."""

    def __init__(self, text: Optional[str] = None, sa: Optional[List[int]] = None):
        self.lcp: List[int] = []
        if text is not None and sa is not None:
            self.build(text, sa)

    def build(self, text: str, sa: List[int]) -> None:
        text = terminate(text)
        n = len(text)
        rank = [0] * n
        for i, s in enumerate(sa):
            rank[s] = i
        lcp = [0] * n
        h = 0
        for i in range(n):
            if rank[i] > 0:
                j = sa[rank[i] - 1]
                while i + h < n and j + h < n and text[i + h] == text[j + h]:
                    h += 1
                lcp[rank[i]] = h
                if h > 0:
                    h -= 1  # key insight of Kasai: h drops by at most 1 per step
            else:
                h = 0
        self.lcp = lcp

    def get(self, i: int) -> int:
        return self.lcp[i]

    def max_lcp(self) -> Tuple[int, int]:
        idx = max(range(len(self.lcp)), key=self.lcp.__getitem__)
        return self.lcp[idx], idx

    def lcp_between(self, i: int, j: int) -> int:
        if not 0 <= i < j < len(self.lcp):
            raise ValueError("need 0 <= i < j < n")
        return min(self.lcp[i + 1:j + 1])

    def __len__(self) -> int:
        return len(self.lcp)


# ===========================================================================
# 3. Burrows-Wheeler Transform
# ===========================================================================
class BWT(BWTADT):
    """BWT string derived from a suffix array; supports inversion."""

    def __init__(self, text: Optional[str] = None, sa: Optional[List[int]] = None):
        self.bwt: str = ""
        if text is not None and sa is not None:
            self.build(text, sa)

    def build(self, text: str, sa: List[int]) -> None:
        text = terminate(text)
        self.bwt = "".join(text[i - 1] for i in sa)  # text[-1] when i == 0

    def get(self, i: int) -> str:
        return self.bwt[i]

    def inverse(self) -> str:
        """
        Inverse BWT via the LF-mapping. The j-th occurrence of c in the first
        column F corresponds to the j-th occurrence of c in the last column L.
        A stable sort of L's indices by character gives exactly that mapping.
        """
        n = len(self.bwt)
        order = sorted(range(n), key=self.bwt.__getitem__)  # stable
        lf = [0] * n
        for j, row in enumerate(order):
            lf[row] = j
        row, out = 0, []  # row 0 is the suffix "$"
        for _ in range(n - 1):
            out.append(self.bwt[row])
            row = lf[row]
        return "".join(reversed(out)) + SENTINEL

    def __len__(self) -> int:
        return len(self.bwt)


# ===========================================================================
# 4. FM-Index
# ===========================================================================
class FMIndex(FMIndexADT):
    """
    FM-index with
      * C table               : C[c] = number of characters in T smaller than c
      * checkpointed Occ      : counts stored every `occ_step` positions
      * sampled suffix array  : only rows whose text position is a multiple of
                                `sa_sample` are stored; others are recovered by
                                walking backwards with the LF-mapping.
    Memory therefore shrinks from n stored SA entries to about n / sa_sample.
    """

    def __init__(self, text: Optional[str] = None, sa: Optional[List[int]] = None,
                 occ_step: int = 64, sa_sample: int = 16):
        self.occ_step = occ_step
        self.sa_sample = sa_sample
        self.bwt = ""
        self.n = 0
        self.C: Dict[str, int] = {}
        self.checkpoints: Dict[str, List[int]] = {}
        self.sa_samples: Dict[int, int] = {}
        if text is not None:
            self.build(text, sa)

    def build(self, text: str, sa: Optional[List[int]] = None) -> None:
        t = terminate(text)
        if sa is None:
            sa = SuffixArray(text).sa
        self.n = len(t)
        self.bwt = BWT(t, sa).bwt
        alphabet = sorted(set(t))

        # C table
        total = 0
        self.C = {}
        for c in alphabet:
            self.C[c] = total
            total += t.count(c)

        # Occ checkpoints: checkpoints[c][b] = count of c in bwt[0 : b*step]
        self.checkpoints = {}
        for c in alphabet:
            cp, run = [0], 0
            for start in range(0, self.n, self.occ_step):
                run += self.bwt.count(c, start, min(start + self.occ_step, self.n))
                cp.append(run)
            self.checkpoints[c] = cp

        # Sampled suffix array (keep only rows with text position % sample == 0)
        self.sa_samples = {row: pos for row, pos in enumerate(sa) if pos % self.sa_sample == 0}

    def occ(self, c: str, i: int) -> int:
        block = i // self.occ_step
        return self.checkpoints[c][block] + self.bwt.count(c, block * self.occ_step, i)

    def backward_search(self, pattern: str) -> Tuple[int, int]:
        lo, hi = 0, self.n
        for c in reversed(pattern):
            if c not in self.C:
                return 0, 0
            lo = self.C[c] + self.occ(c, lo)
            hi = self.C[c] + self.occ(c, hi)
            if lo >= hi:
                return 0, 0
        return lo, hi

    def count(self, pattern: str) -> int:
        lo, hi = self.backward_search(pattern)
        return hi - lo

    def _row_to_position(self, row: int) -> int:
        steps = 0
        while row not in self.sa_samples:
            c = self.bwt[row]
            row = self.C[c] + self.occ(c, row)  # LF-mapping: one step back in the text
            steps += 1
        return self.sa_samples[row] + steps

    def locate(self, pattern: str) -> List[int]:
        lo, hi = self.backward_search(pattern)
        return sorted(self._row_to_position(r) for r in range(lo, hi))

    def stored_sa_entries(self) -> int:
        """How many SA entries are really stored (vs n for a full suffix array)."""
        return len(self.sa_samples)


# ===========================================================================
# 5. Pattern matchers sharing one interface
# ===========================================================================
class NaiveMatcher(PatternMatcherADT):
    """Brute force: try every alignment. O(n*m) worst case, O(1) extra space."""
    name = "Naive"

    def preprocess(self, text: str) -> None:
        self.text = text

    def find_all(self, pattern: str) -> List[int]:
        n, m, t = len(self.text), len(pattern), self.text
        return [i for i in range(n - m + 1) if t[i:i + m] == pattern]


class KMPMatcher(PatternMatcherADT):
    """Knuth-Morris-Pratt: O(n + m) per query, O(m) extra space."""
    name = "KMP"

    def preprocess(self, text: str) -> None:
        self.text = text

    @staticmethod
    def failure_function(pattern: str) -> List[int]:
        fail, k = [0] * len(pattern), 0
        for i in range(1, len(pattern)):
            while k and pattern[i] != pattern[k]:
                k = fail[k - 1]
            if pattern[i] == pattern[k]:
                k += 1
            fail[i] = k
        return fail

    def find_all(self, pattern: str) -> List[int]:
        m = len(pattern)
        if m == 0:
            return []
        fail, res, k = self.failure_function(pattern), [], 0
        for i, c in enumerate(self.text):
            while k and c != pattern[k]:
                k = fail[k - 1]
            if c == pattern[k]:
                k += 1
            if k == m:
                res.append(i - m + 1)
                k = fail[k - 1]
        return res


class SuffixArrayMatcher(PatternMatcherADT):
    """Index-based search: build once O(n log^2 n), then O(m log n + occ) per query."""
    name = "Suffix Array"

    def preprocess(self, text: str) -> None:
        self.index = SuffixArray(text)

    def find_all(self, pattern: str) -> List[int]:
        return self.index.search(pattern)


class FMIndexMatcher(PatternMatcherADT):
    """Compressed index: backward search O(m) counting, locate adds LF-walks."""
    name = "FM-Index"

    def __init__(self, occ_step: int = 64, sa_sample: int = 16):
        self.occ_step, self.sa_sample = occ_step, sa_sample

    def preprocess(self, text: str) -> None:
        self.index = FMIndex(text, occ_step=self.occ_step, sa_sample=self.sa_sample)

    def find_all(self, pattern: str) -> List[int]:
        return self.index.locate(pattern)

    def count(self, pattern: str) -> int:
        return self.index.count(pattern)


ALL_MATCHERS = [NaiveMatcher, KMPMatcher, SuffixArrayMatcher, FMIndexMatcher]


# ======================================================================
# SECTION 3 - Applications
# ======================================================================
from typing import Dict, List, Tuple


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


# ---------------------------------------------------------------------------
# DNA helpers
# ---------------------------------------------------------------------------
def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA string (the opposite strand, read 5'->3')."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """Fraction of bases that are G or C."""
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


# ---------------------------------------------------------------------------
# Repeats and comparison of genomes (suffix array + LCP)
# ---------------------------------------------------------------------------
def longest_repeated_substring(text: str) -> Tuple[str, int, List[int]]:
    """
    Longest substring that occurs at least twice.
    Idea: two suffixes share a prefix of length L iff that prefix repeats, and
    adjacent suffixes in sorted order share the longest prefixes -> max(LCP).
    Returns (substring, length, sorted start positions).
    """
    sa = SuffixArray(text)
    lcp = LCPArray(text, sa.sa)
    length, idx = lcp.max_lcp()
    if length == 0:
        return "", 0, []
    sub = sa.text[sa.sa[idx]:sa.sa[idx] + length]
    return sub, length, sa.search(sub)


def longest_common_substring(a: str, b: str) -> Tuple[str, int, int, int]:
    """
    Longest substring shared by genomes a and b (e.g. two virus strains).
    Build the suffix array of a + '#' + b. Adjacent suffixes coming from
    different genomes with the largest LCP give the answer. '#' occurs once so
    no common prefix can cross it.
    Returns (substring, length, position_in_a, position_in_b).
    """
    if "#" in a or "#" in b:
        raise ValueError("inputs must not contain '#'")
    joined = a + "#" + b
    sa = SuffixArray(joined)
    lcp = LCPArray(joined, sa.sa)
    la = len(a)
    best, best_p, best_q = 0, -1, -1
    for i in range(1, len(sa)):
        p, q = sa.sa[i - 1], sa.sa[i]
        if p == la or q == la:  # the '#' separator itself
            continue
        if (p < la) != (q < la) and lcp.get(i) > best:
            best, best_p, best_q = lcp.get(i), p, q
    if best == 0:
        return "", 0, -1, -1
    pos_a, pos_b = (best_p, best_q - la - 1) if best_p < la else (best_q, best_p - la - 1)
    return a[pos_a:pos_a + best], best, pos_a, pos_b


# ---------------------------------------------------------------------------
# Approximate matching and read alignment
# ---------------------------------------------------------------------------
def hamming(a: str, b: str) -> int:
    """Number of positions where equal-length strings differ."""
    return sum(x != y for x, y in zip(a, b))


def approximate_search(text: str, pattern: str, k: int,
                       matcher: PatternMatcherADT) -> List[Tuple[int, int]]:
    """
    All positions where `pattern` matches `text` with at most k substitutions.
    Pigeonhole principle: if the pattern is cut into k+1 pieces, any match with
    <= k mismatches contains at least one piece EXACTLY. So:
        1. find every piece exactly with the (already preprocessed) matcher,
        2. turn each hit into a candidate start position,
        3. verify candidates with the Hamming distance.
    `matcher` must already be preprocessed on `text`.
    Returns sorted (position, mismatches) pairs.
    """
    m, n = len(pattern), len(text)
    if m > n:
        return []
    if m <= k:  # pieces would be empty -> plain scan
        return [(i, d) for i in range(n - m + 1)
                if (d := hamming(text[i:i + m], pattern)) <= k]

    size = m // (k + 1)
    cuts = [i * size for i in range(k + 1)] + [m]  # last piece takes the remainder
    candidates = set()
    for a, b in zip(cuts, cuts[1:]):
        for hit in matcher.find_all(pattern[a:b]):
            start = hit - a
            if 0 <= start <= n - m:
                candidates.add(start)
    out = []
    for start in sorted(candidates):
        d = hamming(text[start:start + m], pattern)
        if d <= k:
            out.append((start, d))
    return out


def align_reads(reads: List[Dict], matcher: PatternMatcherADT, genome: str,
                max_mismatch: int = 0) -> Dict[str, List[Tuple[int, str, int]]]:
    """
    Map each read onto the genome on both strands.
    reads : list of {"id": str, "seq": str}
    Returns {read_id: [(position, strand, mismatches), ...]} (empty list = unmapped).
    """
    result = {}
    for r in reads:
        hits = []
        for strand, seq in (("+", r["seq"]), ("-", reverse_complement(r["seq"]))):
            for pos, d in approximate_search(genome, seq, max_mismatch, matcher):
                hits.append((pos, strand, d))
        result[r["id"]] = sorted(hits)
    return result


# ======================================================================
# SECTION 4 - Genotype (random DNA) driver
# ======================================================================
import argparse
import json
import os
import random
from typing import Dict, List, Optional


BASES = "ACGT"


def random_genome(length: int, gc_content: float = 0.5, seed: Optional[int] = None) -> str:
    """Random DNA of `length` bases where P(G) = P(C) = gc/2 and P(A) = P(T) = (1-gc)/2."""
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
    """Read a FASTA file into one uppercase ACGT string (headers and other symbols dropped)."""
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(c for c in "".join(seq) if c in BASES)


def generate_test_suite(length: int = 5000, gc_content: float = 0.5, n_reads: int = 30,
                        read_len: int = 50, error_rate: float = 0.02, n_queries: int = 50,
                        pattern_len: int = 15, seed: int = 7) -> Dict:
    """Build a complete random test suite (genome, a mutated strain, reads, queries)."""
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


def generator_cli() -> None:
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

    suite = generate_test_suite(a.length, a.gc, a.reads, a.read_len, a.error, a.queries, seed=a.seed)
    write_fasta(a.out + "_genome.fasta", "synthetic_genome", suite["genome"])
    write_fasta(a.out + "_strain.fasta", "synthetic_strain", suite["strain"])
    meta = {k: v for k, v in suite.items() if k not in ("genome", "strain")}
    with open(a.out + "_reads_queries.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Wrote {a.out}_genome.fasta, {a.out}_strain.fasta, {a.out}_reads_queries.json")
    print(f"Genome length {suite['genome_length']}, GC {suite['gc_content']}, "
          f"{len(suite['reads'])} reads, {len(suite['queries'])} queries")


# ======================================================================
# SECTION 5 - Benchmark
# ======================================================================
import csv
import os
import time
import tracemalloc
from typing import Dict, List


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


# ======================================================================
# SECTION 6 - main()
# ======================================================================
import argparse
import os
import sys
import time
import unittest



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
