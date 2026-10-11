"""
applications.py - Bioinformatics problems solved with the data structures
=========================================================================
    longest_repeated_substring   suffix array + LCP
    longest_common_substring     suffix array + LCP on two genomes
    approximate_search           seed-and-verify (pigeonhole) over any matcher
    align_reads                  map short reads (both strands) onto a genome
    reverse_complement, gc_content   small DNA helpers
"""

from typing import Dict, List, Tuple

from .adt import PatternMatcherADT
from .data_structures import LCPArray, SuffixArray

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
