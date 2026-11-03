"""Test ME contract checks with synthetic result rows."""
from pathlib import Path
import sys, tempfile, unittest
import pandas as pd
BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE/'tools'))
from verify_me import check_table, EXPECTED


class TestMEContract(unittest.TestCase):
    def test_expected_supported_pairs(self):
        rows = [{'config':'C_ASSA','condition':c,'feature':f,'supported':True} for c,f in EXPECTED]
        rows.append({'config':'C_ASSA','condition':'Pooled','feature':'ME1','supported':False})
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'results.csv'
            pd.DataFrame(rows).to_csv(p,index=False)
            self.assertTrue(check_table(p)['matches_frozen_fig4'])


if __name__ == '__main__':
    unittest.main()
