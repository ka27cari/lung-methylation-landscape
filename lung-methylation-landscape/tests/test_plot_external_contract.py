"""Static checks for the Python/R external-projection hand-off."""
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
class ExternalPlotContract(unittest.TestCase):
    def test_external_plot_reads_aligned_export(self):
        s=(ROOT/'r/figures/external_pbmc_validation.R').read_text(encoding='utf8')
        self.assertIn('EXTERNAL_PROGRAM_SCORES.csv',s)
        self.assertIn('EXTERNAL_VALIDATION_SUMMARY.csv',s)
        self.assertIn('source_route <- S$route',s)
        self.assertIn('ext_control <- external_control',s)
        self.assertNotIn('source_iqr / external_iqr',s)

if __name__=='__main__':unittest.main()
