"""
test_all.py - Unit tests
========================
Runs
  * the constant cases in tests/test_cases.json,
  * invariants (sorted suffixes, LCP vs brute force, BWT inverse),
  * randomized cross-checks of all four matchers using the genotype driver.

Run from the project root:   python -m unittest discover -s tests -v
                       or:   python main.py test
"""

import json
import os
import sys
import unittest

# Find the project root so that source modules can be imported.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.applications import (  # noqa: E402
    align_reads,
    approximate_search,
    longest_common_substring,
    longest_repeated_substring,
    reverse_complement,
)
from src.data_structures import (  # noqa: E402
    ALL_MATCHERS,
    BWT,
    FMIndex,
    KMPMatcher,
    LCPArray,
    NaiveMatcher,
    SuffixArray,
    SuffixArrayMatcher,
    terminate,
)
from src.genotype_driver import (  # noqa: E402
    generate_test_suite,
    mutate,
    plant_repeats,
    random_genome,
)

# Load predefined inputs and expected outputs from the JSON test file.
with open(os.path.join(ROOT, "tests", "test_cases.json")) as f:
    CASES = json.load(f)


# This helper finds every exact occurrence of a pattern in a string.
# It serves as a simple reference for checking other search algorithms.
def brute_force(text, pattern):
    return [i for i in range(len(text) - len(pattern) + 1) if text[i:i + len(pattern)] == pattern]


# These tests verify outputs using predefined examples from the JSON file.
# They cover suffix arrays, pattern searching, and substring operations.
class TestConstantCases(unittest.TestCase):
    # Check suffix array, LCP array, and Burrows-Wheeler Transform outputs.
    def test_suffix_array_lcp_bwt(self):
        for c in CASES["suffix_array"]:
            with self.subTest(text=c["text"]):
                sa = SuffixArray(c["text"])
                self.assertEqual(sa.sa, c["expected_sa"])
                self.assertEqual(LCPArray(c["text"], sa.sa).lcp, c["expected_lcp"])
                self.assertEqual(BWT(c["text"], sa.sa).bwt, c["expected_bwt"])

    # Confirm that every matcher finds the expected pattern positions.
    def test_pattern_search_all_matchers(self):
        for c in CASES["pattern_search"]:
            for cls in ALL_MATCHERS:
                with self.subTest(case=c["id"], matcher=cls.name):
                    m = cls()
                    m.preprocess(c["text"])
                    self.assertEqual(m.find_all(c["pattern"]), c["expected"])
                    self.assertEqual(m.count(c["pattern"]), len(c["expected"]))

    # Verify the length and positions of the longest repeated substring.
    def test_longest_repeated_substring(self):
        for c in CASES["longest_repeated_substring"]:
            with self.subTest(text=c["text"]):
                sub, length, positions = longest_repeated_substring(c["text"])
                self.assertEqual(length, c["expected_length"])
                self.assertEqual(len(sub), length)
                if length:
                    self.assertGreaterEqual(len(positions), 2)
                    self.assertEqual(positions, brute_force(c["text"], sub))

    # Check that the longest common substring occurs in both input strings.
    def test_longest_common_substring(self):
        for c in CASES["longest_common_substring"]:
            with self.subTest(a=c["a"], b=c["b"]):
                sub, length, pa, pb = longest_common_substring(c["a"], c["b"])
                self.assertEqual(length, c["expected_length"])
                if length:
                    self.assertEqual(c["a"][pa:pa + length], sub)
                    self.assertEqual(c["b"][pb:pb + length], sub)

    # Compare approximate-search results against the expected test cases.
    def test_approximate_search(self):
        for c in CASES["approximate_search"]:
            with self.subTest(pattern=c["pattern"], k=c["k"]):
                m = SuffixArrayMatcher()
                m.preprocess(c["text"])
                got = [list(x) for x in approximate_search(c["text"], c["pattern"], c["k"], m)]
                self.assertEqual(got, c["expected"])


# These tests check important properties that should always hold.
# A fixed random seed makes the generated test genome repeatable.
class TestInvariants(unittest.TestCase):
    def setUp(self):
        self.text = plant_repeats(random_genome(400, seed=11), motif_len=40, copies=3, seed=12)

    # Verify that suffix-array entries form a permutation and sort suffixes.
    def test_suffixes_are_sorted_and_a_permutation(self):
        sa = SuffixArray(self.text)
        t = sa.text
        self.assertEqual(sorted(sa.sa), list(range(len(t))))
        suffixes = [t[i:] for i in sa.sa]
        self.assertEqual(suffixes, sorted(suffixes))

    # Check LCP values by directly comparing neighboring suffixes.
    def test_lcp_matches_brute_force(self):
        sa = SuffixArray(self.text)
        lcp = LCPArray(self.text, sa.sa)
        t = sa.text
        for i in range(1, len(sa)):
            a, b = t[sa.sa[i - 1]:], t[sa.sa[i]:]
            k = 0
            while k < min(len(a), len(b)) and a[k] == b[k]:
                k += 1
            self.assertEqual(lcp.get(i), k)

    # Confirm that applying the inverse BWT restores the terminated text.
    def test_bwt_inverse_returns_original(self):
        sa = SuffixArray(self.text)
        self.assertEqual(BWT(self.text, sa.sa).inverse(), terminate(self.text))

    # Test LCP queries between suffixes at different ranks.
    def test_lcp_between(self):
        sa = SuffixArray("banana")
        lcp = LCPArray("banana", sa.sa)
        self.assertEqual(lcp.lcp_between(2, 3), 3)  # 'ana$' vs 'anana$'
        self.assertEqual(lcp.lcp_between(1, 3), 1)  # 'a$' vs 'anana$'

    # Verify that FM-index sampling settings preserve search results.
    def test_fm_index_sampling_does_not_change_answers(self):
        sa = SuffixArray(self.text)
        for step, sample in [(1, 1), (4, 3), (64, 16), (7, 50)]:
            fm = FMIndex(self.text, sa.sa, occ_step=step, sa_sample=sample)
            for p in ["ACG", "TTTT", self.text[10:30], "GATTACA"]:
                self.assertEqual(fm.locate(p), brute_force(self.text, p))
                self.assertEqual(fm.count(p), len(brute_force(self.text, p)))

    # Check that the FM-index stores fewer suffix-array entries when sampled.
    def test_fm_index_stores_fewer_sa_entries(self):
        fm = FMIndex(self.text, occ_step=64, sa_sample=16)
        self.assertLess(fm.stored_sa_entries(), len(self.text) // 8)

    # Ensure invalid sentinel usage and empty search bounds raise errors.
    def test_sentinel_validation(self):
        with self.assertRaises(ValueError):
            SuffixArray("AC$GT")
        with self.assertRaises(ValueError):
            SuffixArray("ACGT").bounds("")


# These tests use randomly generated genomes to check algorithm correctness.
# Multiple seeds provide different inputs for additional coverage.
class TestRandomized(unittest.TestCase):
    # Compare all matching algorithms against the brute-force reference.
    def test_all_matchers_agree_on_random_genomes(self):
        for seed in range(5):
            suite = generate_test_suite(length=1500, n_queries=40, pattern_len=8, seed=seed)
            g = suite["genome"]
            matchers = [cls() for cls in ALL_MATCHERS]
            for m in matchers:
                m.preprocess(g)
            for q in suite["queries"] + ["A", "AC", "ACGTACGT"]:
                expected = brute_force(g, q)
                for m in matchers:
                    self.assertEqual(m.find_all(q), expected, f"{m.name} seed={seed} q={q}")

    # Check that reads with an acceptable mismatch count map to their true locations.
    def test_reads_map_back_to_true_position(self):
        suite = generate_test_suite(length=3000, n_reads=25, read_len=40, error_rate=0.03, seed=3)
        g = suite["genome"]
        m = SuffixArrayMatcher()
        m.preprocess(g)
        aln = align_reads(suite["reads"], m, g, max_mismatch=3)
        for r in suite["reads"]:
            if r["mismatches"] <= 3:
                hits = aln[r["id"]]
                strand_pos = [(p, s) for p, s, _ in hits]
                self.assertIn((r["true_pos"], r["strand"]), strand_pos, r["id"])

    # Verify that error-free reads map correctly with zero mismatches.
    def test_exact_reads_map_with_zero_mismatches(self):
        suite = generate_test_suite(length=2000, n_reads=15, read_len=30, error_rate=0.0, seed=9)
        g = suite["genome"]
        m = NaiveMatcher()
        m.preprocess(g)
        aln = align_reads(suite["reads"], m, g, max_mismatch=0)
        for r in suite["reads"]:
            self.assertIn((r["true_pos"], r["strand"], 0), aln[r["id"]])

    # Compare approximate search against direct character-by-character matching.
    def test_approximate_search_matches_brute_force(self):
        g = random_genome(600, seed=21)
        m = KMPMatcher()
        m.preprocess(g)
        pattern = mutate(g[100:124], 0.1, seed=5)
        for k in (0, 1, 2, 3):
            expected = [(i, sum(a != b for a, b in zip(g[i:i + 24], pattern)))
                        for i in range(len(g) - 23)
                        if sum(a != b for a, b in zip(g[i:i + 24], pattern)) <= k]
            self.assertEqual(approximate_search(g, pattern, k, m), expected)

    # Check reverse-complement conversion and applying it twice.
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AACG"), "CGTT")
        self.assertEqual(reverse_complement(reverse_complement("GATTACA")), "GATTACA")


# Run all unit tests and display each test result when this file is executed.
if __name__ == "__main__":
    unittest.main(verbosity=2)
 
