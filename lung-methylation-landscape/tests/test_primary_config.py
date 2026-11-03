"""Verify the configured primary Bayesian adjustment without using patient data."""
import unittest
from pathlib import Path

import yaml


class ConfigContractTests(unittest.TestCase):
    def test_primary_model_keeps_lineage_composition(self):
        path = Path(__file__).resolve().parents[1] / 'config/model.yaml'
        cfg = yaml.safe_load(path.read_text())['bayesian_me']
        primary = next(x for x in cfg['covariate_configs'] if x['id'] == cfg['primary_config'])
        self.assertEqual(primary['blood_pcs'], 0)
        self.assertTrue(all(primary[k] for k in ['age', 'sex', 'smoking_status', 'smoking_amount', 'origin']))
        self.assertEqual(cfg['primary_config'], 'C_ASSA')


if __name__ == '__main__':
    unittest.main()
