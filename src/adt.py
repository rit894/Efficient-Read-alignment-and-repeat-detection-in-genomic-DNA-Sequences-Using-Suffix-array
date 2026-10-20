"""
adt.py - Abstract Data Type (ADT) definitions
=============================================
An ADT says WHAT a data structure does (its operations and their contracts),
not HOW it does it. Every structure used in this case study has an ADT here:

    SuffixArrayADT      - sorted array of all suffix start positions
    LCPArrayADT         - longest-common-prefix values of adjacent suffixes
    BWTADT              - Burrows-Wheeler Transform of a text
    FMIndexADT          - compressed full-text index built on the BWT
    PatternMatcherADT   - common interface so every search method (naive, KMP,
                          suffix array, FM-index) can be benchmarked fairly

Concrete implementations live in data_structures.py.
Conventions used in all docstrings:
    T      the text, always terminated by the sentinel '$' (unique, added
           automatically if missing)
    n      len(T) including the sentinel
    m      len(pattern)
    occ    number of occurrences of the pattern in the text
"""

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
