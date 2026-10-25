"""CLI regression tests for interactive searches and the real-reference demo."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'data' / 'phiX174_NC_001422.1.fasta'


class TestReferenceCLI(unittest.TestCase):
    def run_cli(self, *args, input_text='', expected_code=0):
        result = subprocess.run(
            [sys.executable, '-m', 'main.main', *map(str, args)],
            cwd=ROOT, input=input_text, text=True, capture_output=True,
            timeout=20,
        )
        self.assertEqual(result.returncode, expected_code, result.stdout + result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        return result.stdout + result.stderr

    def test_reference_identity_and_sequence_integrity(self):
        lines = REFERENCE.read_text(encoding='utf-8').splitlines()
        self.assertTrue(lines[0].startswith('>NC_001422.1 '))
        self.assertEqual(sum(line.startswith('>') for line in lines), 1)
        sequence = ''.join(lines[1:])
        self.assertEqual(len(sequence), 5386)
        self.assertLessEqual(set(sequence), set('ACGT'))
        self.assertEqual(
            hashlib.sha256(sequence.encode('ascii')).hexdigest(),
            '97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0',
        )

    def test_all_search_methods_report_overlapping_matches(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            genome = Path(folder) / 'genome.fasta'
            genome.write_text('>overlaps\nAAAA\n', encoding='utf-8')
            for method in ('naive', 'kmp', 'sa', 'fm'):
                with self.subTest(method=method):
                    output = self.run_cli('search', '--fasta', genome,
                                          '--pattern', 'aa', '--method', method)
                    self.assertIn('3 occurrence(s) of AA: [0, 1, 2]', output)

    def test_interactive_search_retries_invalid_patterns(self):
        output = self.run_cli('search', '--fasta', REFERENCE,
                              input_text='\nACN\nacgt\n')
        self.assertEqual(output.count('Please enter a non-empty'), 2)
        self.assertIn('19 occurrence(s) of ACGT', output)

    def test_real_reference_demo_preserves_source(self):
        before = REFERENCE.read_bytes()
        output = self.run_cli('demo', '--fasta', REFERENCE, '--pattern', 'ACGT')
        self.assertIn('Genome length 5386', output)
        self.assertIn('GC content (G + C bases) 44.76%', output)
        self.assertIn('Occurrences in this 5,386-base genome: 19', output)
        self.assertIn('Longest repeated substring: length 12', output)
        self.assertIn('reference vs simulated mutated copy', output)
        self.assertIn('Simulated read alignment', output)
        self.assertEqual(REFERENCE.read_bytes(), before)

    def test_real_reference_demo_prompts_for_a_query(self):
        output = self.run_cli('demo', '--fasta', REFERENCE, input_text='acgt\n')
        self.assertIn('Enter a DNA base or sequence', output)
        self.assertIn('Occurrences in this 5,386-base genome: 19', output)

    def test_synthetic_demo_still_works(self):
        output = self.run_cli('demo', '--pattern', 'GC')
        self.assertIn('PART 2 - Random DNA', output)
        self.assertIn('Genome length 3000', output)
        self.assertIn('Occurrences in this 3,000-base genome: 182', output)

    def test_missing_reference_has_a_readable_error(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            output = self.run_cli('demo', '--fasta', Path(folder) / 'missing.fasta',
                                  '--pattern', 'ACGT', expected_code=2)
            self.assertIn('Error:', output)

    def test_empty_reference_is_rejected(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            genome = Path(folder) / 'empty.fasta'
            genome.write_text('>empty\n', encoding='utf-8')
            output = self.run_cli('demo', '--fasta', genome,
                                  '--pattern', 'ACGT', expected_code=2)
            self.assertIn('no DNA sequence', output)

    def test_invalid_explicit_pattern_is_rejected(self):
        output = self.run_cli('demo', '--fasta', REFERENCE,
                              '--pattern', 'ACN', expected_code=2)
        self.assertIn('contain only A, C, G and T', output)

    def test_eof_during_prompt_has_a_readable_error(self):
        output = self.run_cli('demo', '--fasta', REFERENCE, expected_code=2)
        self.assertIn('no search input received', output)


if __name__ == '__main__':
    unittest.main()
