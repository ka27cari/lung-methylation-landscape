#!/usr/bin/env python3
"""Check terminal-conditioned ME results against the chosen analysis contract."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

EXPECTED = {
    ('T1','ME0'), ('T1','ME3'),
    ('T2','ME2'),
    ('T4','ME0'), ('T4','ME1'),
}


def check_table(path: Path, config: str = 'C_ASSA') -> dict:
    df = pd.read_csv(path)
    selected = df.loc[df['config'].eq(config)]
    if selected.empty:
        raise ValueError(f'Configuration {config} absent from result table')
    supported = selected['supported'].astype(str).str.lower().eq('true')
    found = set(zip(selected.loc[supported, 'condition'], selected.loc[supported, 'feature']))
    pooled_supported = any(c == 'Pooled' for c, _ in found)
    return {'config': config,
            'supported': sorted([f'{c}:{m}' for c,m in found]),
            'matches_frozen_fig4': found == EXPECTED and not pooled_supported}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]/'rebuild')
    args = parser.parse_args()
    out = check_table(args.root/'analysis/tables/TERMINAL_CONDITIONED_ME_HORIZON.csv')
    for k,v in out.items():
        print(k, v)
    return 0 if out['matches_frozen_fig4'] else 1


if __name__=='__main__':
    raise SystemExit(main())
