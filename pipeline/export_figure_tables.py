#!/usr/bin/env python3
"""Export stable plot-ready tables for figures.pptx / visualization workflows.

Only the current clean model is represented. Stable generic plotting filenames are emitted where their semantics are current; terminal-specific files are generated from the current terminal IDs.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1] / "rebuild"
TAB=ROOT/'analysis/tables'; FIG=ROOT/'analysis/figure_data'
MOD=ROOT/'cache/model'; LONG=ROOT/'cache/longitudinal'; EXT=ROOT/'cache/external'
FIG.mkdir(parents=True,exist_ok=True)

spec=json.load(open(TAB/'FINAL_MODEL_SPEC.json'))
terminal_states=[int(x) for x in spec['terminal_states']]

# ---------- utilities ----------
def hbin(t):
    if t < 1: return '0–1'
    if t < 3: return '1–3'
    if t < 5: return '3–5'
    if t < 8: return '5–8'
    if t < 12: return '8–12'
    return '12+'
HORDER=['0–1','1–3','3–5','5–8','8–12','12+']

def binned(df, value_col):
    q=df[df.group.eq('prediagnostic')].copy()
    q['horizon_bin']=q.time.map(hbin)
    rows=[]
    for b in HORDER:
        g=q[q.horizon_bin.eq(b)]
        if len(g):
            rows.append(dict(horizon_bin=b,x=float(g.time.median()),
                             median=float(g[value_col].median()),
                             q25=float(g[value_col].quantile(.25)),
                             q75=float(g[value_col].quantile(.75)),n=int(len(g))))
    return pd.DataFrame(rows)

# ---------- longitudinal terminal mass ----------
pairs=pd.read_csv(FIG/'LONGITUDINAL_TERMINAL_MASS_PAIRS.csv')
scat=pd.read_csv(FIG/'LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_DATA.csv')
binned(scat,'TM').to_csv(FIG/'LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_BINS.csv',index=False)
compat_tm=pairs[['matched_index','time','pred_TM','case_TM']].copy()
compat_tm['pred_sink_weighted_TM']=compat_tm.pred_TM
compat_tm['case_sink_weighted_TM']=compat_tm.case_TM
compat_tm['delta_TM']=compat_tm.case_TM-compat_tm.pred_TM
compat_tm['delta_sink_weighted_TM']=compat_tm.delta_TM
compat_tm.to_csv(FIG/'V13_LONGITUDINAL_SINK_WEIGHTED_TM.csv',index=False)

# ---------- direct copies into figure_data ----------
copy_map={
 'BAYESIAN_TERMINAL_BLOOD_ASSOCIATIONS.csv':'TERMINAL_BLOOD_ASSOCIATIONS.csv',
 'BAYESIAN_TERMINAL_ME_ASSOCIATIONS.csv':'TERMINAL_ME_ASSOCIATIONS.csv',
 'BAYESIAN_ME_PROXIMITY_BY_TERMINAL.csv':'ME_PROXIMITY_BY_TERMINAL.csv',
 'SOURCE_STATE_CONTEXT_SUMMARY.csv':'STATE_CONTEXT_FINGERPRINTS.csv',
 'TERMINAL_CONTEXT_FINGERPRINTS.csv':'TERMINAL_CONTEXT_FINGERPRINTS.csv',
 'TERMINAL_ROUTE_SUMMARY.csv':'TERMINAL_ROUTE_SUMMARY.csv',
 'LANDSCAPE_COARSE_FLOW.csv':'LANDSCAPE_COARSE_FLOW.csv',
 'ANCHOR_RESOLUTION_SOURCE_ONLY_AUDIT.csv':'ANCHOR_SELECTION.csv',
 'PROGRAM_DUAL_RESOLUTION_COVERAGE_CURVE.csv':'PROGRAM_SELECTION.csv',
 'K_SELECTION_STABILITY_AUDIT.csv':'K_SELECTION.csv',
 'POSTERIOR_TAU_CALIBRATION.csv':'TAU_SELECTION.csv',
 'ROUTE_AWARE_TEMPORAL_COMPARISON.csv':'ROUTE_AWARE_TEMPORAL_COMPARISON.csv',
 'PROGRAM_BLOOD_GEOMETRY_ASSOCIATION.csv':'PROGRAM_BLOOD_GEOMETRY_ASSOCIATION.csv',
 'LONGITUDINAL_COUPLING_SUMMARY.csv':'LONGITUDINAL_COUPLING_SUMMARY.csv',
 'FINAL_STATE_ROLE_CLASSIFICATION.csv':'FINAL_STATE_ROLE_CLASSIFICATION.csv',
 'MACROSTATE_MAPPING_FINAL.csv':'MACROSTATE_MAPPING_FINAL.csv',
 'STATE_VISUALIZATION_INDEX.csv':'STATE_VISUALIZATION_INDEX.csv',
 'STATE_TO_TERMINAL_ABSORPTION.csv':'STATE_TO_TERMINAL_ABSORPTION.csv',
 'MACROSTATE_OT_PLAN_FINAL.csv':'MACROSTATE_OT_PLAN_FINAL.csv',
 'MACROSTATE_OT_ROUTES_FINAL.csv':'MACROSTATE_OT_ROUTES_FINAL.csv',
 'TERMINAL_PARTIAL_OT_SENSITIVITY.csv':'TERMINAL_PARTIAL_OT_SENSITIVITY.csv',
 'TERMINAL_ROUTE_HORIZON_SLOPES.csv':'TERMINAL_ROUTE_HORIZON_SLOPES.csv',
 'ROUTE_STRATIFIED_TERMINAL_MASS.csv':'ROUTE_STRATIFIED_TERMINAL_MASS.csv',
 'TERMINAL_ROUTE_TEMPORAL_PROFILES.csv':'TERMINAL_ROUTE_TEMPORAL_PROFILES.csv',
 'ROUTE_WEIGHTED_SCORE_GRID.csv':'ROUTE_WEIGHTED_SCORE_GRID.csv',
 'SINKHORN_NET_CURRENT_EDGES.csv':'SINKHORN_NET_CURRENT_EDGES.csv',
 'SINKHORN_TERMINAL_ROUTE_AFFINITY.csv':'SINKHORN_TERMINAL_ROUTE_AFFINITY.csv',
 'SINKHORN_ROUTE_SENSITIVITY.csv':'SINKHORN_ROUTE_SENSITIVITY.csv',
}
for src,dst in copy_map.items():
    p=TAB/src
    if p.exists(): pd.read_csv(p).to_csv(FIG/dst,index=False)

# ---------- sample-level tables ----------
ext=pd.read_csv(EXT/'EXTERNAL_POSTERIOR.csv')
ext[['sample','group','TM']].to_csv(FIG/'EXTERNAL_TM.csv',index=False)

sp=pd.read_csv(MOD/'SOURCE_POSTERIOR.csv')
keep=[c for c in ['sentrix_id','group','time','type','hard_state','TM'] if c in sp.columns] + [c for c in sp.columns if c.startswith('P_state_')]
sp[keep].to_csv(FIG/'SOURCE_STATE_POSTERIOR.csv',index=False)

pd.read_csv(MOD/'SELECTED_PROGRAMS.csv').to_csv(FIG/'SELECTED_PROGRAMS.csv',index=False)

# ---------- longitudinal macrostate / terminal posterior ----------
lp=pd.read_csv(LONG/'LONGITUDINAL_POSTERIOR.csv')
roles=pd.read_csv(TAB/'FINAL_STATE_ROLE_CLASSIFICATION.csv')
role_map=roles.set_index('state').macrostate_role.to_dict()
state_cols={int(c.split('_')[-1]):c for c in lp.columns if c.startswith('P_state_')}
for role_name,out_name in [
    ('Source','LONGITUDINAL_SOURCE_MACROSTATE_SCATTER_BINS.csv'),
    ('Intermediate','LONGITUDINAL_INTERMEDIATE_MACROSTATE_SCATTER_BINS.csv'),
    ('Terminal','LONGITUDINAL_TERMINAL_MACROSTATE_SCATTER_BINS.csv')]:
    states=[k for k,v in role_map.items() if v==role_name and k in state_cols]
    vals=lp[[state_cols[k] for k in states]].sum(axis=1) if states else np.zeros(len(lp))
    binned(lp.assign(_score=vals),'_score').to_csv(FIG/out_name,index=False)

for t in terminal_states:
    binned(lp.assign(_score=lp[state_cols[t]]),'_score').to_csv(FIG/f'LONGITUDINAL_T{t}_POSTERIOR_SCATTER_BINS.csv',index=False)

# ---------- route-aware data ----------
route=pd.read_csv(LONG/'LONGITUDINAL_ROUTE_SCORES.csv')
qm=lp[['sentrix_id','matched_index','group','time']].merge(route,on='sentrix_id',how='left')
binned(qm,'role_stage').to_csv(FIG/'LONGITUDINAL_ROUTE_STAGE_SCATTER_BINS.csv',index=False)
# robust-trend filename retained; values are horizon-bin medians/IQR, not a fitted smoother.
rtrend=binned(qm,'role_stage').rename(columns={'median':'route_stage_median','q25':'route_stage_q25','q75':'route_stage_q75'})
rtrend.to_csv(FIG/'LONGITUDINAL_ROUTE_STAGE_ROBUST_TREND.csv',index=False)
# Relay-affinity summary for the current landscape.
binned(qm,'eventual_relay_affinity').to_csv(FIG/'LONGITUDINAL_MEDIATION_AFFINITY_SCATTER_BINS.csv',index=False)
for t in terminal_states:
    c=f'route_to_T{t}'
    if c in qm.columns:
        binned(qm,c).to_csv(FIG/f'LONGITUDINAL_ROUTE_T{t}_AFFINITY_SCATTER_BINS.csv',index=False)

# ---------- combined terminal-affinity forest input ----------
comb=[]
n_pred=int((sp.group=='prediagnostic').sum()) if 'group' in sp else np.nan
for domain,src in [('Blood','BAYESIAN_TERMINAL_BLOOD_ASSOCIATIONS.csv'),('ME','BAYESIAN_TERMINAL_ME_ASSOCIATIONS.csv')]:
    z=pd.read_csv(TAB/src)
    for _,r in z.iterrows():
        comb.append(dict(domain=domain,terminal_state=int(r.terminal_state),predictor=r.feature,
                         posterior_median=float(r.posterior_median),
                         **{'CrI_2.5':float(r.CrI_2_5),'CrI_97.5':float(r.CrI_97_5)},
                         n=n_pred,model='source terminal-affinity regression'))
pd.DataFrame(comb).to_csv(FIG/'V13_BAYESIAN_TERMINAL_AFFINITY_ASSOCIATIONS.csv',index=False)

# ---------- plot-file index ----------

idx=[]
for p in sorted(FIG.iterdir()):
    if p.is_file(): idx.append(dict(filename=p.name,relative_path=str(p.relative_to(ROOT)),bytes=p.stat().st_size))
pd.DataFrame(idx).to_csv(FIG/'FIGURE_DATA_INDEX.csv',index=False)

contract={
  'anchor_selection':'ANCHOR_SELECTION.csv',
  'program_selection':'PROGRAM_SELECTION.csv',
  'K_selection':'K_SELECTION.csv',
  'tau_selection':'TAU_SELECTION.csv',
  'source_posterior':'SOURCE_STATE_POSTERIOR.csv',
  'state_roles':'FINAL_STATE_ROLE_CLASSIFICATION.csv',
  'state_visualization_index':'STATE_VISUALIZATION_INDEX.csv',
  'macrostate_mapping':'MACROSTATE_MAPPING_FINAL.csv',
  'state_to_terminal_absorption':'STATE_TO_TERMINAL_ABSORPTION.csv',
  'longitudinal_scatter':'LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_DATA.csv',
  'longitudinal_bins':'LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_BINS.csv',
  'longitudinal_coupling':'LONGITUDINAL_TERMINAL_AWARE_COUPLING.csv',
  'terminal_blood':'TERMINAL_BLOOD_ASSOCIATIONS.csv',
  'terminal_ME':'TERMINAL_ME_ASSOCIATIONS.csv',
  'terminal_routes':'TERMINAL_ROUTE_SUMMARY.csv',
  'sinkhorn_net_current_edges':'SINKHORN_NET_CURRENT_EDGES.csv',
  'sinkhorn_terminal_route_affinity':'SINKHORN_TERMINAL_ROUTE_AFFINITY.csv',
  'sinkhorn_route_sensitivity':'SINKHORN_ROUTE_SENSITIVITY.csv',
  'external_TM':'EXTERNAL_TM.csv',
  'n_programs':int(spec['n_programs']),'K':int(spec['K']),
  'terminal_states':terminal_states,'tau':float(spec['tau'])
}
json.dump(contract,open(FIG/'FIGURE_DATA_CONTRACT.json','w'),indent=2)
print('figure-data export complete:',len(list(FIG.glob('*'))),'files')
