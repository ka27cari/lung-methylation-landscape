#!/usr/bin/env python3
"""Execute the methylation landscape pipeline on locally staged input data."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent / 'rebuild'
REPO = ROOT.parent
# The order here, not numeric filename prefixes, defines the pipeline.
STAGES = (
    ('clean', 'pipeline/clean_generated_outputs.py'),
    ('anchors', 'pipeline/resolve_reference_anchors.py'),
    ('discovery', 'pipeline/discover_recurrent_programs.py'),
    ('programs', 'pipeline/select_stable_programs.py'),
    ('landscape', 'pipeline/fit_microstate_landscape.py'),
    ('posterior', 'pipeline/calibrate_terminal_posterior.py'),
    ('transport', 'pipeline/infer_sinkhorn_transport.py'),
    ('downstream', 'pipeline/analyze_landscape_outcomes.py'),
    ('context', 'pipeline/fit_terminal_blood_associations.py'),
    ('extensions', 'pipeline/summarize_validation_extensions.py'),
    ('oof', 'pipeline/compare_out_of_fold_models.py'),
    ('figdata', 'pipeline/export_figure_tables.py'),
    ('diagnostic-plots', 'pipeline/render_analysis_diagnostics.py'),
    ('validate', 'pipeline/check_analysis_outputs.py'),
)

REQUIRED_INPUTS = (
    'nested.csv', 'longitudinal.csv', 'm_tr.csv',
    'm_te.csv', 'ph-ext.csv', 'm-ext.csv',
)
PRIMARY_MODEL = 'C_ASSA'

STAGE_DESCRIPTIONS = {
    'clean': 'Remove generated output (preserve NMF checkpoint by default)',
    'anchors': 'Resolve stable reference CpGs',
    'discovery': 'Discover recurrent NMF programs (computationally expensive)',
    'programs': 'Select stable program basis',
    'landscape': 'Fit microstate landscape and terminal states',
    'posterior': 'Calibrate terminal membership probabilities',
    'transport': 'Estimate net transport between microstates',
    'downstream': 'Analyze trajectories and external cohort',
    'context': 'Estimate terminal-conditioned methylation associations',
    'extensions': 'Prepare supplementary model summaries',
    'oof': 'Compare models out of fold',
    'figdata': 'Write numerical figure tables',
    'diagnostic-plots': 'Render internal diagnostics (not manuscript artwork)',
    'validate': 'Check generated analysis outputs',
}



def stage_sequence(recompute_discovery: bool):
    checkpoint = ROOT / 'cache/discovery/RECURRENT_FAMILIES_ALL.csv'
    return [
        (label, script)
        for label, script in STAGES
        if label != 'discovery' or recompute_discovery or not checkpoint.is_file()
    ]


def health():
    """Report missing scripts/inputs and confirm the primary model selection."""
    missing_scripts = [
        script for _, script in STAGES if not (REPO / script).is_file()
    ]
    if not (REPO / 'data/derived_evidence/OOF_ANCHOR_SEGMENT_SUMMARY.csv').is_file():
        missing_scripts.append('data/derived_evidence/OOF_ANCHOR_SEGMENT_SUMMARY.csv')
    missing_inputs = [
        filename for filename in REQUIRED_INPUTS
        if not (ROOT / 'data' / filename).is_file()
    ]
    config_file = REPO / 'config/model.yaml'
    text = config_file.read_text() if config_file.is_file() else ''
    match = re.search(r'^\s*primary_config:\s*(\w+)\s*$', text, re.MULTILINE)
    primary = match.group(1) if match else None

    print('Missing scripts:', ', '.join(missing_scripts) or 'none')
    print('Missing private inputs:', ', '.join(missing_inputs) or 'none')
    print('Primary Bayesian configuration:', primary or 'not found')
    print('External projection:', 'source-control median/IQR alignment in Python')
    return missing_scripts, missing_inputs, primary


def run(recompute_discovery: bool):
    environment = os.environ.copy()
    for name in (
        'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS',
    ):
        environment[name] = '1'
    for label, script in stage_sequence(recompute_discovery):
        args = [sys.executable, str(REPO / script)]
        if label == 'clean' and recompute_discovery:
            args.append('--recompute-discovery')
        print(f'RUN {label}: {script}', flush=True)
        subprocess.run(args, cwd=ROOT, env=environment, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['list', 'check', 'run'])
    parser.add_argument('--recompute-discovery', action='store_true')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--allow-clean', action='store_true')
    args = parser.parse_args()

    if args.command == 'list':
        for label, script in stage_sequence(args.recompute_discovery):
            print(f'{label:>16}  {script}\n                  {STAGE_DESCRIPTIONS[label]}')
        return 0

    missing_scripts, missing_inputs, primary = health()
    ready = not missing_scripts and not missing_inputs and primary == PRIMARY_MODEL
    if args.command == 'check':
        return 0 if ready else 2

    if not args.execute or not args.allow_clean:
        print('Run requires both --execute and --allow-clean to overwrite generated results.')
        return 2
    if not ready:
        print('Execution blocked by preflight checks.')
        return 2
    run(args.recompute_discovery)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
