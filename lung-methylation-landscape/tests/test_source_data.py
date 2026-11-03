"""Check the figure source-data exporter contract without study data."""
from pathlib import Path
import unittest
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from build_source_data import FIGURES, PANEL_GAPS


class ExportContractTests(unittest.TestCase):
    def test_figure_coverage_contract(self):
        expected={'Fig1','Fig2','Fig3','Fig4','Fig5','SuppFig1','SuppFig2','SuppFig3','SuppFig4'}
        self.assertEqual(set(FIGURES),expected)
        self.assertEqual(set(PANEL_GAPS),expected)
        for fig,entries in FIGURES.items():
            labels=[label for panel,label,spec in entries]
            self.assertEqual(len(labels),len(set(labels)),fig)
            for panel,label,spec in entries:
                self.assertIn(spec.split('/')[0],{'tables','model','derived'})

    def test_no_sample_level_source_tables(self):
        forbidden={'SOURCE_STATE_POSTERIOR.csv','LONGITUDINAL_FATE_HEMO_DISPLACEMENT.csv',
                   'LONGITUDINAL_TERMINAL_MASS_PAIRS.csv','EXTERNAL_TM.csv',
                   'SOURCE_PROGRAM_SCORES.csv','EXTERNAL_PROGRAM_SCORES.csv'}
        for entries in FIGURES.values():
            for _,_,spec in entries:
                self.assertNotIn(spec.split('/')[-1], forbidden)


if __name__=='__main__':unittest.main()
