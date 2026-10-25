"""
data_structures.py - Concrete implementations of every ADT in adt.py
====================================================================
    SuffixArray        (prefix doubling)
    LCPArray           (Kasai's algorithm)
    BWT                (from the suffix array, with inverse transform)
    FMIndex            (BWT + C table + checkpointed Occ + sampled SA)
    NaiveMatcher, KMPMatcher, SuffixArrayMatcher, FMIndexMatcher
"""

from bisect import bisect_left, bisect_right
from typing import Dict, List, Optional, Tuple

from .adt import (
    BWTADT,
    FMIndexADT,
    LCPArrayADT,
    PatternMatcherADT,
    SuffixArrayADT,
)

SENTINEL = "$"


# Ensure the text ends with one unique marker used by the indexing algorithms.
def terminate(text: str) -> str:
    """Append the unique sentinel '$' if it is not already the last character."""
    # Remove an existing final marker first so it is not duplicated.
    body = text[:-1] if text.endswith(SENTINEL) else text
    if SENTINEL in body:
        raise ValueError("text must not contain the sentinel '$' except as its last character")
    return body + SENTINEL


# ===========================================================================
# 1. Suffix Array
# ===========================================================================
# Build a suffix array by repeatedly sorting longer prefixes of each suffix.
def _prefix_doubling(s: str) -> List[int]:
    """
    Manber-Myers style prefix doubling.
    Round k sorts suffixes by their first 2k characters using the ranks from the
    previous round as (first half, second half) pairs packed into one integer.
    Rounds: O(log n). Each round sorts n integers: O(n log n)  ->  O(n log^2 n).
    """
    n = len(s)
    # Start with every suffix represented by its starting index.
    sa = list(range(n))
    # Initial ranks are based on the character codes.
    rank = [ord(c) for c in s]
    k = 1
    while True:
        base = max(rank) + 2
        # Combine the current rank and the rank k positions ahead into one sortable key.
        keys = [rank[i] * base + (rank[i + k] + 1 if i + k < n else 0) for i in range(n)]
        # Sorting these keys updates the order of suffixes.
        sa.sort(key=keys.__getitem__)
        # Assign fresh ranks according to the newly sorted keys.
        new_rank = [0] * n
        r = 0
        # Compare neighboring suffix keys to decide whether their ranks differ.
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

    # Prepare the text and construct its suffix order.
    def build(self, text: str) -> None:
        self.text = terminate(text)
        self.sa = _prefix_doubling(self.text)

    def get(self, i: int) -> int:
        return self.sa[i]

    # Find the interval of suffix-array rows that begin with the pattern.
    def bounds(self, pattern: str) -> Tuple[int, int]:
        if not pattern:
            raise ValueError("pattern must be non-empty")
        m, text = len(pattern), self.text
        # Only the pattern-length prefix is needed for binary searching.
        key = lambda i: text[i:i + m]  # first m chars of the suffix at SA row
        lo = bisect_left(self.sa, pattern, key=key)
        hi = bisect_right(self.sa, pattern, key=key)
        return lo, hi

    # Return matching text positions in ascending order.
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
# The LCP array stores common-prefix lengths for adjacent sorted suffixes.
class LCPArray(LCPArrayADT):
    """LCP array built with Kasai's linear-time algorithm."""

    def __init__(self, text: Optional[str] = None, sa: Optional[List[int]] = None):
        self.lcp: List[int] = []
        if text is not None and sa is not None:
            self.build(text, sa)

    def build(self, text: str, sa: List[int]) -> None:
        text = terminate(text)
        n = len(text)
        # Map each text position to its row in the suffix array.
        rank = [0] * n
        for i, s in enumerate(sa):
            rank[s] = i
        lcp = [0] * n
        h = 0
        # Visit suffixes in text order while reusing the previous match length.
        for i in range(n):
            if rank[i] > 0:
                j = sa[rank[i] - 1]
                # Extend the common prefix while both suffixes still match.
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
# The Burrows-Wheeler Transform stores the character preceding each sorted suffix.
class BWT(BWTADT):
    """BWT string derived from a suffix array; supports inversion."""

    def __init__(self, text: Optional[str] = None, sa: Optional[List[int]] = None):
        self.bwt: str = ""
        if text is not None and sa is not None:
            self.build(text, sa)

    def build(self, text: str, sa: List[int]) -> None:
        text = terminate(text)
        # A suffix starting at index zero wraps around to the final character.
        self.bwt = "".join(text[i - 1] for i in sa)  # text[-1] when i == 0

    def get(self, i: int) -> str:
        return self.bwt[i]

    # Reconstruct the original text by following the last-to-first mapping.
    def inverse(self) -> str:
        """
        Inverse BWT via the LF-mapping. The j-th occurrence of c in the first
        column F corresponds to the j-th occurrence of c in the last column L.
        A stable sort of L's indices by character gives exactly that mapping.
        """
        n = len(self.bwt)
        # Stable sorting groups characters while preserving their occurrence order.
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
# The FM-index supports compact pattern searching using the BWT and sampled positions.
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
        # The alphabet contains each distinct character in sorted order.
        alphabet = sorted(set(t))

        # C maps each character to the number of text characters that sort before it.
        # C table
        total = 0
        self.C = {}
        for c in alphabet:
            self.C[c] = total
            total += t.count(c)

        # Checkpoints avoid recounting the entire BWT for every occurrence query.
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

    # Narrow the matching suffix-array interval from right to left through the pattern.
    def backward_search(self, pattern: str) -> Tuple[int, int]:
        if not pattern:
            raise ValueError("pattern must be non-empty")
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

    # Recover a text position by walking until a sampled suffix-array row is reached.
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
# The naive matcher checks every possible starting position directly.
class NaiveMatcher(PatternMatcherADT):
    """Brute force: try every alignment. O(n*m) worst case, O(1) extra space."""
    name = "Naive"

    def preprocess(self, text: str) -> None:
        self.text = text

    def find_all(self, pattern: str) -> List[int]:
        if not pattern:
            raise ValueError("pattern must be non-empty")
        n, m, t = len(self.text), len(pattern), self.text
        return [i for i in range(n - m + 1) if t[i:i + m] == pattern]


class KMPMatcher(PatternMatcherADT):
    """Knuth-Morris-Pratt: O(n + m) per query, O(m) extra space."""
    name = "KMP"

    def preprocess(self, text: str) -> None:
        self.text = text

    @staticmethod
    # Build the KMP table that tells us how far to fall back after a mismatch.
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
            if not pattern:
                raise ValueError("pattern must be non-empty")
        fail, res, k = self.failure_function(pattern), [], 0
        # Scan the text once, reusing matched-prefix information after mismatches.
        for i, c in enumerate(self.text):
            while k and c != pattern[k]:
                k = fail[k - 1]
            if c == pattern[k]:
                k += 1
            if k == m:
                res.append(i - m + 1)
                k = fail[k - 1]
        return res


# This matcher builds a suffix-array index once and reuses it for queries.
class SuffixArrayMatcher(PatternMatcherADT):
    """Index-based search: build once O(n log^2 n), then O(m log n + occ) per query."""
    name = "Suffix Array"

    def preprocess(self, text: str) -> None:
        self.index = SuffixArray(text)

    def find_all(self, pattern: str) -> List[int]:
        return self.index.search(pattern)


# This matcher wraps the FM-index behind the shared pattern-matcher interface.
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
