#!/usr/bin/env python3
"""Structural and numerical regression checks for  release."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1] / "rebuild";TAB=ROOT/'analysis/tables';MOD=ROOT/'cache/model';DRAW=ROOT/'cache/posterior_draws';FIG=ROOT/'analysis/figure_data'
problems=[]
anchor=json.load(open(TAB/'ANCHOR_RESOLUTION_POLICY.json')); prog=json.load(open(TAB/'PROGRAM_SELECTION_POLICY.json')); ksel=json.load(open(TAB/'K_SELECTION_POLICY.json')); spec=json.load(open(TAB/'FINAL_MODEL_SPEC.json')); basin=json.load(open(TAB/'TERMINAL_COHERENT_BASIN_POLICY.json'))
if abs(float(anchor['selected_anchor_fraction'])-.02)>1e-12: problems.append('anchor != 0.02')
anchor_audit=pd.read_csv(TAB/'ANCHOR_RESOLUTION_SOURCE_ONLY_AUDIT.csv')
anchor_min=float(anchor_audit.loc[anchor_audit.counterfactual_PCA30_participation_ratio.idxmin(),'anchor_fraction'])
if abs(anchor_min-float(anchor['selected_anchor_fraction']))>1e-12: problems.append(f'anchor is not global minimum participation-ratio grid point: min={anchor_min}')
if anchor.get('criterion')!='counterfactual_PCA30_participation_ratio' or anchor.get('direction')!='minimize': problems.append('anchor selection policy is not simplified participation-ratio minimization')
if bool(anchor.get('selection_uses_adjacent_representation_correlation',True)): problems.append('adjacent representation correlation unexpectedly used for anchor selection')
if bool(anchor.get('selection_uses_control_reconstruction_MAE',True)): problems.append('control reconstruction MAE unexpectedly used for anchor selection')
if int(prog['n_selected_programs'])!=9: problems.append('selected programs != 9')
if int(prog.get('stable_dictionary_size',prog.get('n_stable_distinct_pool',-1)))!=14: problems.append('stable program dictionary != 14')
if int(prog.get('compact_landscape_basis_size',prog.get('n_selected_programs',-1)))!=9: problems.append('compact landscape basis != 9')
if bool(prog.get('downstream_endpoints_used_for_basis_size_selection',True)): problems.append('downstream endpoint unexpectedly used for P9 basis selection')
curve=pd.read_csv(TAB/'PROGRAM_DUAL_RESOLUTION_COVERAGE_CURVE.csv')
u8=int(curve.loc[curve.n_programs.eq(8),'uncovered_n'].iloc[0]);u9=int(curve.loc[curve.n_programs.eq(9),'uncovered_n'].iloc[0])
if not (u8==27 and u9==7):problems.append(f'P8/P9 source noncoverage changed {u8},{u9}')
thrs=pd.read_csv(TAB/'PROGRAM_DUAL_COVERAGE_LIMIT_10_20.csv')
if len(thrs)!=11 or not thrs.selected_basis_size.eq(9).all():problems.append('10-20 uncovered sample limits do not all select P9')
bmap=pd.read_csv(TAB/'PROGRAM_DICTIONARY_BASIS_MAP.csv')
if len(bmap)!=14 or int(bmap.selected_in_compact_basis.astype(bool).sum())!=9: problems.append('program dictionary/basis map is not 14 -> 9')
if int(ksel['selected_K'])!=16: problems.append('K != 16')
kstab=pd.read_csv(TAB/'K_SELECTION_STABILITY_AUDIT.csv');r16=kstab[kstab.K.eq(16)].iloc[0];r17=kstab[kstab.K.eq(17)].iloc[0]
if not (bool(r16.passes_typical_support) and not bool(r17.passes_typical_support)): problems.append('K16 is not the finest typically-supported resolution')
kscan=pd.read_csv(TAB/'K_START_ADMISSIBILITY_SCAN.csv'); n16=int(((kscan.K==16)&kscan.passes).sum())
if n16!=42: problems.append(f'K16 source-admissible starts !=42: {n16}')
if list(map(int,basin['selected_starts']))!=[34,38,45,56,61,62,74,95]: problems.append('terminal-coherent start basin changed')
if int(basin['selected_cut'])!=4: problems.append('coarsest coherent cut changed')
for key in ['residual_time_used','longitudinal_used','external_used','blood_used','ME_used','histology_used','ARI_used','requested_terminal_count_used']:
    if bool(basin.get(key,False)): problems.append(f'{key} unexpectedly used in basin selection')
if list(map(int,spec['terminal_states']))!=[2,6,7,8]: problems.append('terminal states changed')
if abs(float(spec['tau'])-.20)>1e-12: problems.append('tau != 0.20')
post=pd.read_csv(MOD/'SOURCE_POSTERIOR.csv');pc=[c for c in post if c.startswith('P_state_')];err=float(np.max(np.abs(post[pc].sum(axis=1).to_numpy()-1)))
if err>1e-8: problems.append(f'posterior normalization error {err}')
roles=pd.read_csv(TAB/'FINAL_STATE_ROLE_CLASSIFICATION.csv'); terms=sorted(roles.loc[roles.final_microstate_terminal.astype(bool),'state'].astype(int).tolist())
if terms!=[2,6,7,8]: problems.append('role-table terminal states changed')
# Sinkhorn route is downstream interpretation only.
rs=json.load(open(TAB/'SINKHORN_ROUTE_LAYER_SPEC.json'))
if list(map(int,rs.get('terminal_states',[])))!=[2,6,7,8]: problems.append('Sinkhorn terminal IDs changed')
if abs(float(rs.get('epsilon',np.nan))-.10)>1e-12: problems.append('Sinkhorn epsilon !=0.10')
if int(rs.get('retained_net_current_edges',-1))!=33: problems.append('Sinkhorn net95 retained edge count !=33')
if not bool(rs.get('sinkhorn_converged',False)): problems.append('Sinkhorn did not converge')
# Bayesian ME rigor checks.
me=pd.read_csv(TAB/'TERMINAL_CONDITIONED_ME_HORIZON.csv'); ma=pd.read_csv(TAB/'ME_COVARIATE_CONFIG_AUDIT.csv'); primary=ma.loc[ma.primary_config.astype(bool),'config'].tolist()
if primary!=['C_ASSA']: problems.append(f'primary ME covariate config != C_ASSA: {primary}')
for forbidden in ['FDR','p_value','pvalue','interaction_omnibus_p','Wald']:
    if any(forbidden.lower() in c.lower() for c in me.columns): problems.append(f'frequentist significance column in Bayesian ME table: {forbidden}')
pm=me[me.primary_config.astype(bool)]
if pm[pm.condition.eq('Pooled')].supported.astype(bool).any(): problems.append('primary pooled ME proximity effect unexpectedly supported')
if not pm[~pm.condition.eq('Pooled')].supported.astype(bool).any(): problems.append('no terminal-conditioned ME proximity effect supported')
if not {'Pr_gt0','Pr_lt0','Pr_direction','CrI_2.5','CrI_97.5','supported'}.issubset(me.columns): problems.append('Bayesian posterior evidence fields missing')
# numerical regression targets
summary=json.load(open(TAB/'ANALYSIS_SUMMARY.json'))
expected={'source_control_median':0.0007963560645184901,'source_far_median':0.021300148071155014,'source_near_median':0.5681812025860795,'long_far_median':0.24553821977140952,'long_near_median':0.6646813237916472,'long_case_median':0.9311709877801614,'terminal_aware_coupling_rho':0.667004048582996,'pooled_coupling_rho':0.32246963562753034,'external_auc':0.7114285714285714}
for k,v in expected.items():
    if k not in summary or abs(float(summary[k])-v)>1e-8: problems.append(f'numerical regression changed: {k}={summary.get(k)} expected {v}')
# expensive/core artifacts
for f in ['ADMISSIBLE_COASSIGNMENT.npz','TERMINAL_CONDITIONED_ME_DRAWS.npz','BAYESIAN_ME_PROXIMITY_DRAWS.npz','LONGITUDINAL_COUPLING_DRAWS.npz','BLOOD_GEOMETRY_PERMUTATION_DRAWS.npz','EXTERNAL_SOURCE_ADAPTER.npz']:
    if not (DRAW/f).exists(): problems.append('missing object '+f)
for f in ['SINKHORN_OT_PLAN.csv','SINKHORN_NET_CURRENT.csv','SINKHORN_NET_CURRENT_95.csv','SOURCE_ROUTE_SCORES.csv']:
    if not (MOD/f).exists(): problems.append('missing model output '+f)
for f in ['LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_DATA.csv','LONGITUDINAL_TERMINAL_AWARE_COUPLING.csv','EXTERNAL_TM.csv','FIGURE_DATA_CONTRACT.json','FIGURE_DATA_INDEX.csv']:
    if not (FIG/f).exists(): problems.append('missing figure data '+f)
out={'status':'pass' if not problems else 'fail','anchor_fraction':anchor['selected_anchor_fraction'],'anchor_selection_rule':anchor.get('rule'),'anchor_selected_participation_ratio':anchor.get('selected_participation_ratio'),'stable_program_dictionary_size':int(prog.get('stable_dictionary_size',prog.get('n_stable_distinct_pool'))),'n_programs':prog['n_selected_programs'],'compact_landscape_basis_size':int(prog.get('compact_landscape_basis_size',prog.get('n_selected_programs'))),'P8_dual_uncovered':u8,'P9_dual_uncovered':u9,'K':ksel['selected_K'],'K16_source_admissible_starts':n16,'selected_basin_starts':basin['selected_starts'],'terminal_states':spec['terminal_states'],'tau':spec['tau'],'primary_ME_config':primary[0] if primary else None,'sinkhorn_net95_edges':rs.get('retained_net_current_edges'),'posterior_max_sum_error':err,'problems':problems}
json.dump(out,open(TAB/'BUILD_VALIDATION.json','w'),indent=2)
print(json.dumps(out,indent=2))
if problems: raise SystemExit(1)
