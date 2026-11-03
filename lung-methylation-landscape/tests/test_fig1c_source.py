import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from ingest_fig1c_counts import ingest

class TestFig1cCounts(unittest.TestCase):
    def test_excludes_unknown_and_updates_manifest_idempotently(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source"
            source.mkdir()
            with (source / "DATASETS.csv").open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["figure", "panel", "dataset", "file", "rows", "upstream", "sha256"])
            input_csv = root / "counts.csv"
            with input_csv.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["eligibility", "diagnosis_interval", "ldct_category", "n"])
                w.writerows([
                    ["Eligible", "<3 yrs", "Nodule present", 9],
                    ["Ineligible", ">=3 yrs", "Normal", 45],
                    ["Unknown", "<3 yrs", "Nodule present", 1],
                ])
            self.assertEqual(ingest(input_csv, source), (2, 54))
            self.assertEqual(ingest(input_csv, source), (2, 54))
            with (source / "DATASETS.csv").open(newline="") as f:
                self.assertEqual(len(list(csv.DictReader(f))), 1)
            with (source / "Fig1/ldct_screening_timing_counts.csv").open(newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(sum(int(r["n"]) for r in rows), 54)
            self.assertFalse(any(r["eligibility"] == "Unknown" for r in rows))

if __name__ == "__main__":
    unittest.main()
