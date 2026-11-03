"""Public entrypoints refer to the current, descriptive file layout."""
from pathlib import Path
import sys
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import run


class RepositoryLayoutTests(unittest.TestCase):
    def test_every_pipeline_stage_resolves(self):
        self.assertEqual(len(run.STAGES), len(set(label for label, _ in run.STAGES)))
        for _, script in run.STAGES:
            self.assertTrue((REPO / script).is_file(), script)
            self.assertEqual((REPO / script).parent, REPO / 'pipeline')

    def test_renamed_r_environment_paths_resolve(self):
        self.assertTrue((REPO / 'r/analysis/analysis_environment.R').is_file())
        self.assertTrue((REPO / 'r/figures/figure_environment.R').is_file())
        self.assertTrue((REPO / 'r/figures/external_pbmc_validation.R').is_file())
        for f in (REPO / 'r').rglob('*'):
            if f.suffix.lower() == '.r':
                source = f.read_text(encoding='utf8')
                self.assertNotIn('"preprocessing", "cohorts", "00_environment.r"', source, f.name)
                self.assertNotIn('"plots", "common.R"', source, f.name)

    def test_original_private_inputs_not_tracked(self):
        self.assertIn('/rebuild/', (REPO/'.gitignore').read_text(encoding='utf8'))


if __name__ == '__main__':
    unittest.main()
