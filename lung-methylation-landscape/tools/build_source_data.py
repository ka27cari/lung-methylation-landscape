#!/usr/bin/env python3
"""Export figure-labelled numerical tables from the completed private analysis run.

Exports contain aggregate results only. Participant-level posterior matrices,
clinical metadata, methylation measurements and paired displacements are excluded.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

FIGURES = {
    'Fig1': [],
    'Fig2': [],
    'Fig3': [
        ('d', 'anchor_selection', 'tables/ANCHOR_RESOLUTION_SOURCE_ONLY_AUDIT.csv'),
        ('e', 'program_coverage', 'tables/PROGRAM_DUAL_RESOLUTION_COVERAGE_CURVE.csv'),
        ('e', 'program_dictionary', 'tables/PROGRAM_DICTIONARY_BASIS_MAP.csv'),
        ('f', 'program_loading_cosine', 'derived/program_loading_cosine.csv'),
    ],
    'Fig4': [
        ('b', 'terminal_mass_group_summaries', 'tables/ROUTE_STRATIFIED_TERMINAL_MASS.csv'),
        ('c-d', 'state_context_aggregates', 'tables/SOURCE_STATE_CONTEXT_SUMMARY.csv'),
        ('e', 'terminal_route_flows', 'tables/TERMINAL_ROUTE_SUMMARY.csv'),
        ('e', 'sinkhorn_net_current', 'tables/SINKHORN_NET_CURRENT_EDGES.csv'),
        ('f', 'program_state_centroids', 'model/STATE_CENTROIDS.csv'),
        ('h', 'matched_pair_oof_aucs', 'tables/INCREMENTAL_INFORMATION_OOF_SUMMARY.csv'),
        ('i', 'terminal_conditioned_ME', 'tables/BAYESIAN_ME_PROXIMITY_BY_TERMINAL.csv'),
    ],
    'Fig5': [
        ('a', 'terminal_mass_group_summaries', 'tables/ROUTE_STRATIFIED_TERMINAL_MASS.csv'),
        ('b', 'paired_route_transitions', 'derived/paired_route_transitions.csv'),
        ('c', 'baseline_state_future_route_AUCs', 'derived/baseline_state_future_route_AUCs.csv'),
        ('d', 'terminal_blood_associations', 'tables/BAYESIAN_TERMINAL_BLOOD_ASSOCIATIONS.csv'),
        ('e', 'terminal_aware_longitudinal_summary', 'tables/LONGITUDINAL_COUPLING_SUMMARY.csv'),
        ('f', 'external_PBMC_calibrated_summary', 'tables/EXTERNAL_VALIDATION_SUMMARY.csv'),
        ('g', 'route_program_direction_contrasts', 'derived/route_program_direction_contrasts.csv'),
        ('h', 'preferred_route_group_counts', 'derived/preferred_route_group_counts.csv'),
    ],
    'SuppFig1': [],
    'SuppFig2': [
        ('a', 'program_loading_cosine', 'derived/program_loading_cosine.csv'),
    ],
    'SuppFig3': [
        ('a', 'posterior_temperature_calibration', 'tables/POSTERIOR_TAU_CALIBRATION.csv'),
        ('b', 'K_selection_stability', 'tables/K_SELECTION_STABILITY_AUDIT.csv'),
        ('c', 'terminal_consensus', 'tables/SELECTED_BASIN_TERMINAL_SUPPORT.csv'),
        ('d', 'sinkhorn_net_current', 'tables/SINKHORN_NET_CURRENT_EDGES.csv'),
        ('e', 'terminal_route_affinity', 'tables/SINKHORN_TERMINAL_ROUTE_AFFINITY.csv'),
        ('f', 'route_score_grid', 'tables/ROUTE_WEIGHTED_SCORE_GRID.csv'),
    ],
    'SuppFig4': [
        ('a', 'longitudinal_terminal_mass_group_summaries', 'tables/ROUTE_STRATIFIED_TERMINAL_MASS.csv'),
        ('b', 'K_selection_stability', 'tables/K_SELECTION_STABILITY_AUDIT.csv'),
        ('c', 'route_program_direction_contrasts', 'derived/route_program_direction_contrasts.csv'),
    ],
}

# These panels lack the underlying plotted observations or require the R-only
# analysis objects and/or consent for participant-level source-data release.
PANEL_GAPS = {
    'Fig1': 'a,b,c,d,e,f',
    'Fig2': 'a,b,c,d,e,f,g',
    'Fig3': 'a,b,c,g,h; panels d-f only partially represented',
    'Fig4': 'a,b,c,d,e,f,g,h,i: aggregate support only, original plotted points incomplete',
    'Fig5': 'a,d,e,f: aggregate support only, plotted individual points withheld; c,g,h derived points/counts provided',
    'SuppFig1': 'flowchart construction data and source cohort files are not included',
    'SuppFig2': 'b,c,d,e,f: full plotted data require R-derived objects',
    'SuppFig3': 'b,c,d,e,f: summary-level support, not all plot-level values',
    'SuppFig4': 'a,b,c: per-specimen or CpG-level underlying values still needed',
}


def build_derived(root: Path, dest: Path) -> None:
    model = root / 'cache/model'
    external = root / 'cache/external'
    long = root / 'cache/longitudinal'
    tables = root / 'analysis/tables'
    dest.mkdir(parents=True, exist_ok=True)

    h = np.load(model / 'SELECTED_PROGRAM_LOADINGS.npz', allow_pickle=False)['loadings']
    h = h / np.maximum(np.linalg.norm(h, axis=1, keepdims=True), 1e-12)
    cos = h @ h.T
    pd.DataFrame([
        {'program_1':f'P{i+1}', 'program_2':f'P{j+1}', 'loading_cosine':float(cos[i,j])}
        for i in range(9) for j in range(9)
    ]).to_csv(dest / 'program_loading_cosine.csv',index=False)

    terms = [2,6,7,8]
    routes = ['T1','T2','T3','T4']
    prog = [f'program_{j:02d}' for j in range(1,10)]
    S = pd.read_csv(model / 'SOURCE_PROGRAM_SCORES.csv')
    P = pd.read_csv(model / 'SOURCE_POSTERIOR.csv').set_index('sentrix_id').loc[S.sentrix_id].reset_index()
    R = pd.read_csv(model / 'SOURCE_ROUTE_SCORES.csv').set_index('sentrix_id').loc[S.sentrix_id].reset_index()
    E = pd.read_csv(external / 'EXTERNAL_PROGRAM_SCORES.csv')
    EP = pd.read_csv(external / 'EXTERNAL_POSTERIOR.csv').set_index('sample').loc[E['sample']].reset_index()
    A = pd.read_csv(tables / 'STATE_TO_TERMINAL_ABSORPTION.csv').sort_values('state')
    SR = np.asarray(routes)[np.argmax(R[[f'route_to_T{t}' for t in terms]].to_numpy(float),axis=1)]
    ER = np.asarray(routes)[np.argmax(EP[[f'P_state_{i:02d}' for i in range(16)]].to_numpy(float) @ A[[f'T{t}' for t in terms]].to_numpy(float),axis=1)]
    Xs = S[prog].to_numpy(float)
    Xe = E[prog].to_numpy(float)
    sc = P.group.eq('control').to_numpy()
    sp = P.group.eq('prediagnostic').to_numpy()
    ec = E.group.eq('Control').to_numpy()
    ep = E.group.eq('LUAD').to_numpy()
    sc_mean = Xs[sc].mean(0)
    ec_mean = Xe[ec].mean(0)
    rows=[]
    for route in routes:
        sm = sp & (SR == route)
        em = ep & (ER == route)
        if not (sm.any() and em.any()):
            raise ValueError(f'Missing specimens for external direction comparison, {route}')
        ds = Xs[sm].mean(0)-sc_mean
        de = Xe[em].mean(0)-ec_mean
        for j in range(9):
            rows.append({'route':route,'program':f'P{j+1}',
                'source_prediagnostic_minus_control':float(ds[j]),
                'external_LUAD_minus_control':float(de[j]),
                'direction_agreement':bool(np.sign(ds[j])==np.sign(de[j]))})
    directions = pd.DataFrame(rows)
    assert int(directions.direction_agreement.sum())==31 and len(directions)==36
    directions.to_csv(dest/'route_program_direction_contrasts.csv',index=False)

    counts=[]
    for stage,groupings,rt in [
        ('Discovery',[('Control',sc),('Prediagnostic',sp)],SR),
        ('External',[('Control',ec),('LUAD',ep)],ER),
    ]:
        for group, mask in groupings:
            n=int(mask.sum())
            for route in routes:
                value=int((mask & (rt==route)).sum())
                counts.append({'cohort':stage,'group':group,'route':route,'n':value,'group_n':n,'proportion':value/n})
    pd.DataFrame(counts).to_csv(dest/'preferred_route_group_counts.csv',index=False)

    lm=pd.read_csv(long/'LONGITUDINAL_POSTERIOR.csv')
    lr=pd.read_csv(long/'LONGITUDINAL_ROUTE_SCORES.csv')
    match=lr.merge(lm[['sentrix_id','matched_index','group']],on='sentrix_id',validate='one_to_one')
    match['T1/T2']=match['route_to_T2']+match['route_to_T6']
    match['T3']=match['route_to_T7']
    match['T4']=match['route_to_T8']
    branches=['T1/T2','T3','T4']
    m=match[match[branches].apply(np.isfinite).all(axis=1)].copy()
    m=m[m[branches].sum(axis=1)>0]
    m['dominant']=np.select([(m['T1/T2']>=m.T3)&(m['T1/T2']>=m.T4),
       (m.T3>=m['T1/T2'])&(m.T3>=m.T4)],branches[:2],default='T4')
    pair=m.pivot(index='matched_index',columns='group',values='dominant').dropna(subset=['prediagnostic','case'])
    assert len(pair)==39 and int((pair.prediagnostic==pair.case).sum())==24
    trans=[]
    for after in branches:
        subset=pair[pair.case==after]
        for before in branches:
            n=int((subset.prediagnostic==before).sum())
            trans.append({'diagnosis_route':after,'prediagnostic_route':before,'n':n,
                'diagnosis_route_total':len(subset),'fraction_within_diagnosis_route':n/len(subset) if len(subset) else 0})
    pd.DataFrame(trans).to_csv(dest/'paired_route_transitions.csv',index=False)

    roles=pd.read_csv(tables/'FINAL_STATE_ROLE_CLASSIFICATION.csv')
    states=roles[(roles.n>=5) | roles.macrostate_role.eq('Terminal')].copy()
    states['role_rank']=states.macrostate_role.map({'Source':0,'Intermediate':1,'Terminal':2})
    states=states.sort_values(['role_rank','state'])
    states['j']=states.groupby('macrostate_role').cumcount()+1
    states['state_label']=states.macrostate_role.map({'Source':'B','Intermediate':'I','Terminal':'T'})+states.j.astype(str)
    labels=states[states.macrostate_role.eq('Terminal')].set_index('state').state_label.to_dict()
    aff=pd.read_csv(tables/'SINKHORN_TERMINAL_ROUTE_AFFINITY.csv')
    aff=aff.merge(pd.DataFrame({'terminal_state':list(labels),'terminal_label':list(labels.values())}),on='terminal_state')
    aff['branch']=aff.terminal_label.replace({'T1':'T1/T2','T2':'T1/T2'})
    pre=lm[lm.group.eq('prediagnostic')]
    case=m[m.group.eq('case')][['matched_index','dominant']]
    pre=pre.merge(case,on='matched_index',validate='one_to_one')
    assert len(pre)==39
    results=[]
    for _,state in states[states.macrostate_role.isin(['Source','Intermediate'])].iterrows():
        xs=aff[aff.state==state.state]
        val=xs.groupby('branch').net95_absorption_affinity.sum()
        weights=np.array([val.get(b,0.) for b in branches])
        if weights.sum()<=0:continue
        weights/=weights.sum()
        values=pre[f'P_state_{int(state.state):02d}'].to_numpy(float)
        for j,br in enumerate(branches):
            mask=pre.dominant.eq(br).to_numpy()
            xa=values[mask];xb=values[~mask]
            auc=(rankdata(values)[mask].sum()-len(xa)*(len(xa)+1)/2)/(len(xa)*len(xb))
            pv=mannwhitneyu(xa,xb,alternative='greater',method='asymptotic').pvalue
            selected=bool(weights[j]>=.60 and abs(weights[j]-weights.max())<1e-12)
            results.append({'state':int(state.state),'state_label':state.state_label,
              'route_branch':br,'n_branch':len(xa),'n_other':len(xb),
              'AUC':float(auc),'p':float(pv),'source_affinity':float(weights[j]),'selected':selected})
    aucdf=pd.DataFrame(results)
    aucdf['q']=multipletests(aucdf.p,method='fdr_bh')[1]
    aucdf['supported']=aucdf.selected & (aucdf.q<.05)&(aucdf.AUC>.5)
    assert int(aucdf.supported.sum())==5
    aucdf.to_csv(dest/'baseline_state_future_route_AUCs.csv',index=False)


def materialize(root: Path, out: Path) -> None:
    t=root/'analysis/tables'
    dest=out/'derived'
    build_derived(root,dest)
    out.mkdir(parents=True,exist_ok=True)
    source={'tables':t,'model':root/'cache/model','derived':dest}
    index=[]
    for figure,entries in FIGURES.items():
        folder=out/figure
        folder.mkdir(parents=True,exist_ok=True)
        for panel,label,spec in entries:
            kind,filename=spec.split('/',1)
            src=source[kind]/filename
            if not src.is_file():
                raise FileNotFoundError(f'{figure}/{label} <- {src}')
            target=folder/(label+'.csv')
            shutil.copyfile(src,target)
            with target.open(newline='') as fh:
                n=max(0,sum(1 for _ in fh)-1)
            index.append({'figure':figure,'panel':panel,'dataset':label,'file':str(target.relative_to(out)),
                          'rows':n,'upstream':spec,
                          'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
    with (out/'DATASETS.csv').open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=list(index[0]));w.writeheader();w.writerows(index)
    with (out/'COVERAGE.csv').open('w',newline='',encoding='utf-8') as fh:
        w=csv.writer(fh);w.writerow(['figure','remaining_requirement'])
        for fig,gap in PANEL_GAPS.items():w.writerow([fig,gap])
    (out/'README.md').write_text('''# Figure numerical data: author-review candidate\n\n'
'Aggregate and derived numerical tables are organized by figure and labelled by panel.\n'
'Figure 1, Figure 2 and Supplementary Figure 1 have no tables because their\n'
'original R-only plotting inputs were not supplied. The included aggregates\n'
'do not replace specimen-level plotted points, and several panels remain\n'
'incomplete (see COVERAGE.csv). This is not the final journal Source Data file.\n'
'No original sample identifiers, clinical dates or methylation matrices are\n'
'included. Data custodians must authorize public release of these derivatives.\n'
'No participant-level samples are redistributed.\n'''.replace("'\n'",''),encoding='utf-8')
    print('Exported',len(index),'aggregate datasets', 'across',len(FIGURES),'figure folders')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True,help='Private completed rebuild workspace')
    p.add_argument('--out',type=Path,required=True,help='Output directory outside Git repository')
    a=p.parse_args()
    materialize(a.root.resolve(),a.out.resolve())

if __name__=='__main__':main()
