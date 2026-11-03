#!/usr/bin/env python3
"""Remove generated outputs before a clean rebuild.

The expensive rNMF discovery checkpoint is preserved by default. Pass
--recompute-discovery to remove it as well and rebuild directly from the six
input files.
"""
from pathlib import Path
import argparse, shutil

ROOT = Path(__file__).resolve().parents[1] / "rebuild"
parser = argparse.ArgumentParser()
parser.add_argument('--recompute-discovery', action='store_true')
args = parser.parse_args()

for rel in [
    'cache/model',
    'cache/longitudinal',
    'cache/external',
    'cache/posterior_draws',
    'analysis/tables',
    'analysis/figure_data',
    'analysis/figures',
    'analysis/spec',
]:
    p = ROOT / rel
    if p.exists():
        shutil.rmtree(p)
    p.mkdir(parents=True, exist_ok=True)

if args.recompute_discovery:
    p = ROOT / 'cache/discovery'
    if p.exists():
        shutil.rmtree(p)
    p.mkdir(parents=True, exist_ok=True)

print('Cleaned generated outputs.' + (' Discovery will be recomputed.' if args.recompute_discovery else ' Discovery checkpoint preserved.'))
