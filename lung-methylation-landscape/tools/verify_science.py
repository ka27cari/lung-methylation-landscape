#!/usr/bin/env python3
"""Compare aggregate figure contracts against a locally staged study run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from verify_external import verify as verify_external
from verify_me import check_table

REPO = Path(__file__).resolve().parents[1]
EXPECTED = json.loads((REPO / 'tests/expected_aggregates.json').read_text())


def route_program_concordance(root: Path) -> tuple[int, int]:
    """Count route-level program direction matches across source and external samples."""
    model = root / 'cache/model'
    ext_dir = root / 'cache/external'
    src = pd.read_csv(model / 'SOURCE_PROGRAM_SCORES.csv')
    source_post = pd.read_csv(model / 'SOURCE_POSTERIOR.csv').set_index('sentrix_id').loc[src.sentrix_id]
    route = pd.read_csv(model / 'SOURCE_ROUTE_SCORES.csv').set_index('sentrix_id').loc[src.sentrix_id]
    ext = pd.read_csv(ext_dir / 'EXTERNAL_PROGRAM_SCORES.csv')
    ext_post = pd.read_csv(ext_dir / 'EXTERNAL_POSTERIOR.csv').set_index('sample').loc[ext['sample']]
    absorption = pd.read_csv(root / 'analysis/tables/STATE_TO_TERMINAL_ABSORPTION.csv').sort_values('state')

    terminals = [int(c.removeprefix('route_to_T')) for c in route if c.startswith('route_to_T')]
    programs = [c for c in src if c.startswith('program_')]
    src_route = np.argmax(route[[f'route_to_T{i}' for i in terminals]].to_numpy(), axis=1)
    ext_posterior = ext_post[[f'P_state_{j:02d}' for j in range(len(absorption))]].to_numpy()
    ext_route = np.argmax(
        ext_posterior @ absorption[[f'T{i}' for i in terminals]].to_numpy(), axis=1
    )
    x = src[programs].to_numpy()
    y = ext[programs].to_numpy()
    src_control = source_post.group.eq('control').to_numpy()
    ext_control = ext.group.eq('Control').to_numpy()
    src_baseline = x[src_control].mean(axis=0)
    ext_baseline = y[ext_control].mean(axis=0)

    n_match = n_total = 0
    for j in range(len(terminals)):
        src_cases = source_post.group.eq('prediagnostic').to_numpy() & (src_route == j)
        ext_cases = ext.group.eq('LUAD').to_numpy() & (ext_route == j)
        if not src_cases.any() or not ext_cases.any():
            raise ValueError(f'No case samples in route {j + 1}')
        source_change = x[src_cases].mean(axis=0) - src_baseline
        external_change = y[ext_cases].mean(axis=0) - ext_baseline
        n_match += int(np.sum(np.sign(source_change) == np.sign(external_change)))
        n_total += len(programs)
    return n_match, n_total


def check(root: Path) -> list[str]:
    issues: list[str] = []
    summary = json.loads((root / 'analysis/tables/ANALYSIS_SUMMARY.json').read_text())
    for key, expected_key in (
        ('external_auc', 'external_auc_control_aligned'),
        ('terminal_aware_coupling_rho', 'terminal_aware_longitudinal_spearman'),
        ('pooled_coupling_rho', 'pooled_longitudinal_spearman'),
        ('blood_geometry_rho', 'program_blood_geometry_spearman'),
        ('blood_geometry_permutation_p', 'program_blood_geometry_permutation_p'),
    ):
        if abs(float(summary[key]) - EXPECTED[expected_key]) > 1e-8:
            issues.append(f'Unexpected {key}: {summary[key]}')
    ext = verify_external(root)
    if abs(ext['auroc_aligned'] - EXPECTED['external_auc_control_aligned']) > 1e-8:
        issues.append(f'External AUROC changed: {ext}')
    me = check_table(root / 'analysis/tables/TERMINAL_CONDITIONED_ME_HORIZON.csv')
    if not me['matches_frozen_fig4']:
        issues.append(f'Bayesian terminal ME changed: {me}')
    agree = route_program_concordance(root)
    if f'{agree[0]}/{agree[1]}' != EXPECTED['program_by_route_direction_concordance']:
        issues.append(f'Route/program concordance changed: {agree}')
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=REPO / 'rebuild')
    args = parser.parse_args()
    problems = check(args.root)
    if problems:
        print('\n'.join(problems))
        return 1
    print('Aggregate scientific contracts: PASS')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
