#!/usr/bin/env python3
"""Final interpretation extensions for the frozen  model.

This script performs no model selection and no frequentist ME significance
screen. Route outputs use the frozen Sinkhorn net-current layer. ME summaries
are copied from the Bayesian Normal-Inverse-Gamma analysis in fit_terminal_blood_associations.py.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1] / "rebuild"
TAB=ROOT/'analysis/tables'; FIG=ROOT/'analysis/figure_data'; MOD=ROOT/'cache/model'; EXT=ROOT/'cache/external'
for p in [TAB,FIG]: p.mkdir(parents=True,exist_ok=True)
roles=pd.read_csv(TAB/'FINAL_STATE_ROLE_CLASSIFICATION.csv')
state=pd.read_csv(TAB/'STATE_TABLE.csv') if (TAB/'STATE_TABLE.csv').exists() else pd.DataFrame()
terms=list(map(int,json.load(open(TAB/'FINAL_MODEL_SPEC.json'))['terminal_states']))

# Source-state structure, preserving legacy filename.
if len(state):
    struc=state.copy()
    if 'role' not in struc.columns:
        struc=struc.merge(roles[['state','macrostate_role']].rename(columns={'macrostate_role':'role'}),on='state',how='left')
    struc.to_csv(TAB/'SOURCE_STATE_STRUCTURE.csv',index=False)
    struc.to_csv(FIG/'SOURCE_STATE_STRUCTURE.csv',index=False)

# Fine route paths from Sinkhorn net95, not the discovery LP95 plan.
J=pd.read_csv(MOD/'SINKHORN_NET_CURRENT_95.csv',index_col=0).to_numpy(float)
role=dict(zip(roles.state.astype(int),roles.macrostate_role.astype(str)))
edges={(i,j):float(J[i,j]) for i in range(J.shape[0]) for j in range(J.shape[1]) if i!=j and J[i,j]>1e-15}
sources=[s for s,r in role.items() if r=='Source']
route_rows=[]
def dfs(path,target,max_edges=5):
    cur=path[-1]
    if len(path)>1 and cur==target:
        vals=[edges[(path[i],path[i+1])] for i in range(len(path)-1)]
        route_rows.append(dict(terminal_state=target,path='→'.join(f'S{x}' for x in path),
                               role_sequence='→'.join(role.get(x,'?')[0] for x in path),n_edges=len(vals),
                               path_capacity=min(vals),flow_product=float(np.prod(vals)),
                               edge_masses=';'.join(f'{v:.8g}' for v in vals),transport_basis='Sinkhorn_net95'))
        return
    if len(path)-1>=max_edges:return
    for (a,b),w in edges.items():
        if a!=cur or b in path:continue
        if b in terms and b!=target:continue
        dfs(path+[b],target,max_edges)
for t in terms:
    for s in sources:dfs([s],t,5)
routes=pd.DataFrame(route_rows)
if len(routes):
    routes=routes.sort_values(['terminal_state','path_capacity','flow_product'],ascending=[True,False,False])
    routes['route_rank']=routes.groupby('terminal_state').cumcount()+1
    routes=routes[routes.route_rank<=20]
routes.to_csv(TAB/'TERMINAL_STATE_ROUTE_PATHS.csv',index=False); routes.to_csv(FIG/'TERMINAL_STATE_ROUTE_PATHS.csv',index=False)

# External terminal direction, descriptive only.
ext_path=EXT/'EXTERNAL_POSTERIOR.csv'
if ext_path.exists():
    ext=pd.read_csv(ext_path)
    rows=[]
    for i,t in enumerate(terms):
        c=f'P_state_{t:02d}'
        if c in ext:
            rows.append(dict(terminal_state=t,terminal_label=f'T{i+1}',control_median=float(ext.loc[ext.group.astype(str).str.lower().eq('control'),c].median()),
                             LUAD_median=float(ext.loc[ext.group.eq('LUAD'),c].median()),
                             median_difference=float(ext.loc[ext.group.eq('LUAD'),c].median()-ext.loc[ext.group.astype(str).str.lower().eq('control'),c].median())))
    pd.DataFrame(rows).to_csv(TAB/'EXTERNAL_TERMINAL_DIRECTION.csv',index=False)
    pd.DataFrame(rows).to_csv(FIG/'EXTERNAL_TERMINAL_DIRECTION.csv',index=False)

# Bayesian-only ME conditioning compatibility summaries. No p-values/FDR.
me_audit=pd.read_csv(TAB/'ME_COVARIATE_CONFIG_AUDIT.csv')
me_all=pd.read_csv(TAB/'TERMINAL_CONDITIONED_ME_HORIZON.csv')
me_audit.to_csv(TAB/'ME_CONDITIONING_MODEL_COMPARISON.csv',index=False)
me_all[me_all.primary_config.astype(bool)].to_csv(TAB/'ME_CONDITIONING_INTERACTION_COEFFICIENTS.csv',index=False)
pd.DataFrame([dict(method='Bayesian Normal-Inverse-Gamma',
                        primary_config=str(me_audit.loc[me_audit.primary_config.astype(bool),'config'].iloc[0]),
                        posterior_support='95% CrI excludes 0 and posterior sign probability >=0.975',
                        FDR_or_Wald_used=False,
                        selection_adjusted=False,
                        note='Primary model is specified in the configuration; alternatives are sensitivity models')]).to_csv(TAB/'ME_CONDITIONING_SUMMARY.csv',index=False)

# Program-basis sensitivity compatibility summary if a current P14 build exists.
p14=ROOT/'analysis/sensitivity/P14/SUMMARY.csv'
comp=[]
prim=pd.read_csv(TAB/'ANALYSIS_SUMMARY.csv').iloc[0] if (TAB/'ANALYSIS_SUMMARY.csv').exists() else None
if prim is not None:
    comp.append(dict(model='Primary P9/K16 terminal-coherent',n_programs=9,K=16,n_terminal=len(terms),
                     source_far_median=float(prim.source_far_median),source_near_median=float(prim.source_near_median),
                     long_far_median=float(prim.long_far_median),long_near_median=float(prim.long_near_median),
                     long_case_median=float(prim.long_case_median),external_auc=float(prim.external_auc)))
if p14.exists():
    s=pd.read_csv(p14).iloc[0]
    comp.append(dict(model='P14/K16 basis sensitivity',n_programs=14,K=16,n_terminal=int(s.n_terminal),
                     source_far_median=float(s.source_far_median),source_near_median=float(s.source_near_median),
                     long_far_median=float(s.long_far_median),long_near_median=float(s.long_near_median),
                     long_case_median=float(s.long_case_median),external_auc=float(s.external_auc)))
pd.DataFrame(comp).to_csv(TAB/'P14_PROGRAM_BASIS_SENSITIVITY.csv',index=False)

# Explicit transport-basis note.
json.dump(dict(discovery_transport='exact partial LP95 used only for source terminal-support selection',
               route_transport='balanced Sinkhorn net current, epsilon=0.10, top 95% net-current mass',
               terminal_states_frozen=terms,
               ME_inference='Bayesian Normal-Inverse-Gamma only; no OLS/Wald/FDR used for the primary ME claim'),
          open(TAB/'FINAL_INTERPRETATION_LAYER_SPEC.json','w'),indent=2)
print('final interpretation extensions written; ME inference is Bayesian-only')
