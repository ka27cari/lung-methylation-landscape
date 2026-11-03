#!/usr/bin/env python3
"""Aggregate-only check of external-control calibration and terminal mass."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'src'))
from pfl.external_projection import frozen_projection


def verify(root: Path):
    src = pd.read_csv(root / 'cache/model/SOURCE_PROGRAM_SCORES.csv')
    post = pd.read_csv(root / 'cache/model/SOURCE_POSTERIOR.csv')
    raw_path = root / 'cache/external/EXTERNAL_PROGRAM_SCORES_RAW.csv'
    ext = pd.read_csv(raw_path if raw_path.is_file() else root / 'cache/external/EXTERNAL_PROGRAM_SCORES.csv')
    runtime = np.load(root / 'cache/model/MODEL_RUNTIME_OUTPUT.npz', allow_pickle=False)
    program = [f'program_{j:02d}' for j in range(1, 10)]
    controls = post.set_index('sentrix_id').loc[src.sentrix_id, 'group'].eq('control').to_numpy()
    ext_controls = ext.group.eq('Control').to_numpy()
    y = ext.group.eq('LUAD').to_numpy()
    xs = src[program].to_numpy(float)
    xe = ext[program].to_numpy(float)
    result = frozen_projection(xe, xs, controls, ext_controls,
                               runtime['centroids'], runtime['state_scales'],
                               float(runtime['tau'][0]), runtime['terminal_states'].astype(int))
    score_raw = result['raw_terminal_mass']
    score_aligned = result['terminal_mass']
    return {'n_controls': int((~y).sum()), 'n_cases': int(y.sum()),
            'auroc_unaligned': float(roc_auc_score(y, score_raw)),
            'auroc_aligned': float(roc_auc_score(y, score_aligned))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=BASE/'rebuild')
    p.add_argument('--expected-aligned', type=float, default=.7114285714285714)
    p.add_argument('--tolerance', type=float, default=1e-8)
    a = p.parse_args()
    r = verify(a.root)
    for key, value in r.items():
        print(f'{key}: {value}')
    return 0 if abs(r['auroc_aligned'] - a.expected_aligned) < a.tolerance else 1


if __name__ == '__main__':
    raise SystemExit(main())
