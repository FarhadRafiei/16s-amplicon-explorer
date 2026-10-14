import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analyze_16s import jaccard, kmers, parse_fasta


class TestAnalysis(unittest.TestCase):
    def test_kmers_and_jaccard(self):
        self.assertEqual(kmers("ACGT", 2), {"AC", "CG", "GT"})
        self.assertEqual(jaccard({"A", "B"}, {"B", "C"}), 1 / 3)

    def test_demo_reads_are_parseable(self):
        records = parse_fasta(Path(__file__).parents[1] / "data/demo/reads.fasta")
        self.assertEqual(len(records), 9)


if __name__ == "__main__":
    unittest.main()
