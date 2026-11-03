#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1] / "rebuild"
TAB=ROOT/'analysis/tables'; FIGD=ROOT/'analysis/figure_data'; OUT=ROOT/'analysis/figures';OUT.mkdir(parents=True,exist_ok=True)

def save(fig,name):
    fig.tight_layout();fig.savefig(OUT/name,dpi=220,bbox_inches='tight');fig.savefig(OUT/name.replace('.png','.pdf'),bbox_inches='tight');plt.close(fig)

# 1 Anchor: normal-like CpG rank segments plus compact reference-adjusted PCA30 structure
seg=pd.read_csv(Path(__file__).resolve().parents[1]/'data/derived_evidence/OOF_ANCHOR_SEGMENT_SUMMARY.csv')
a=pd.read_csv(TAB/'ANCHOR_RESOLUTION_SOURCE_ONLY_AUDIT.csv')
fig,axs=plt.subplots(1,3,figsize=(13.6,3.9))
for ax,col,ylab,title in [(axs[0],'heldout_control_MAD_mean','Control CpG MAD','A | Normal control variation'),
                          (axs[1],'heldout_abs_case_control_shift_mean','Absolute matched shift','B | Case–control shift')]:
    ax.bar(seg.segment.astype(str),seg[col]);ax.tick_params(axis='x',labelrotation=30);ax.set(ylabel=ylab,title=title);ax.grid(axis='y',alpha=.2)
axs[2].plot(a.anchor_fraction*100,a.counterfactual_PCA30_participation_ratio,'o-')
q=a[np.isclose(a.anchor_fraction,.02)].iloc[0]
axs[2].scatter([2],[q.counterfactual_PCA30_participation_ratio],s=110,marker='*');axs[2].set(xlabel='Anchor fraction (%)',ylabel='PCA30 effective dimension',title='C | Reference-adjusted structure');axs[2].grid(alpha=.2)
fig.suptitle('Normal-like anchors and compact residual structure')
save(fig,'01_anchor_resolution.png')
# 2 Programs: source samples uncovered in both discovery resolutions, first crossing at P9
x=pd.read_csv(TAB/'PROGRAM_DUAL_RESOLUTION_COVERAGE_CURVE.csv')
fig,ax=plt.subplots(figsize=(7,4.4));ax.plot(x.n_programs,x.uncovered_n,'o-')
ax.axhspan(10,20,alpha=.15);ax.axvline(9,linestyle='--',linewidth=1)
ax.scatter([9],[int(x.loc[x.n_programs.eq(9),'uncovered_n'].iloc[0])],s=100,marker='D')
ax.set(xlabel='Reproducibility-ordered program prefix',ylabel='Doubly uncovered source samples',
       title='P9 by source dual-resolution discovery-context coverage')
for n in [8,9,10,11]:
    r=x[x.n_programs.eq(n)]
    if len(r):ax.annotate(f'P{n}: {int(r.uncovered_n.iloc[0])}',(n,int(r.uncovered_n.iloc[0])),xytext=(3,8),textcoords='offset points')
ax.grid(alpha=.2);save(fig,'02_program_selection.png')

# 3 K resolution: typical terminal support and admissible fraction
x=pd.read_csv(TAB/'K_SELECTION_STABILITY_AUDIT.csv');pol=json.load(open(TAB/'K_SELECTION_POLICY.json'));sel=int(pol['selected_K']);thr=float(pol['typical_terminal_support_min']);fig,ax=plt.subplots(figsize=(6.8,4.3));ax.plot(x.K,x.median_smallest_terminal_pred_support,marker='o',label='Median smallest terminal support');ax.axhline(thr,linestyle='--',linewidth=1,label='Support threshold');ax.axvline(sel,linestyle='--',linewidth=1);ax.set(xlabel='K',ylabel='Prediagnostic samples',title=f'Landscape resolution (selected K={sel})');ax2=ax.twinx();ax2.plot(x.K,x.admissible_fraction,marker='s',linewidth=1,label='Admissible fraction');ax2.set_ylabel('Admissible-start fraction');lines=ax.get_lines()+ax2.get_lines();labels=[l.get_label() for l in lines if not l.get_label().startswith('_')];lines=[l for l in lines if not l.get_label().startswith('_')];ax.legend(lines,labels,frameon=False,loc='best');ax.grid(alpha=.2);save(fig,'03_K_resolution.png')

# 4 tau: source temporal contrast under control fidelity
x=pd.read_csv(TAB/'POSTERIOR_TAU_CALIBRATION.csv');fig,ax=plt.subplots(figsize=(6.7,4.2));ax.plot(x.tau,x.loo_binary_nll,marker='o',label='LOO terminal/non-terminal NLL');sel=x[x.selected].iloc[0];ax.axvline(sel.tau,linestyle='--',linewidth=1);ax.set(xlabel='Posterior softness τ',ylabel='Leave-one-out binary negative log-likelihood',title=f'Source-only posterior calibration (τ={sel.tau:.2f})');ax.legend(frameon=False);ax.grid(alpha=.2);save(fig,'04_tau_calibration.png')

# 5 longitudinal TM scatter
x=pd.read_csv(FIGD/'LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_DATA.csv');fig,ax=plt.subplots(figsize=(7.4,4.8));pre=x[x.group.eq('prediagnostic')];case=x[x.group.eq('case')];ax.scatter(pre.time,pre.TM,s=28,label='Prediagnostic');ax.scatter(case.time,case.TM,s=28,marker='x',label='Case');
for mid,g in x.groupby('matched_index'):
    if len(g)==2:
        g=g.sort_values('time',ascending=False);ax.plot(g.time,g.TM,linewidth=.55,alpha=.25)
ax.set(xlabel='Years to diagnosis',ylabel='Terminal mass',title='Longitudinal terminal occupancy');ax.invert_xaxis();ax.legend(frameon=False);ax.grid(alpha=.2);save(fig,'05_longitudinal_terminal_mass.png')

# 6 coupling
x=pd.read_csv(FIGD/'LONGITUDINAL_TERMINAL_AWARE_COUPLING.csv');fig,ax=plt.subplots(figsize=(5.2,4.7));ax.scatter(x.terminal_aware_fate_displacement,x.terminal_aware_hemo_displacement,s=34);b=np.polyfit(x.terminal_aware_fate_displacement,x.terminal_aware_hemo_displacement,1);xx=np.linspace(x.terminal_aware_fate_displacement.min(),x.terminal_aware_fate_displacement.max(),100);ax.plot(xx,b[0]*xx+b[1],linewidth=1);summ=pd.read_csv(TAB/'LONGITUDINAL_COUPLING_SUMMARY.csv').iloc[0];ax.set(xlabel='Terminal-aware methylation displacement',ylabel='Terminal-aware hematopoietic displacement',title=f'Longitudinal coupling (Spearman ρ={summ.terminal_aware_rho:.2f})');ax.grid(alpha=.2);save(fig,'06_longitudinal_terminal_aware_coupling.png')

# 7 external
x=pd.read_csv(FIGD/'EXTERNAL_TM.csv');groups=['Control','LUAD'];vals=[x.loc[x.group.eq(g),'TM'].to_numpy() for g in groups];fig,ax=plt.subplots(figsize=(5.0,4.6));ax.boxplot(vals,tick_labels=groups,showfliers=False);rng=np.random.default_rng(1)
for i,v in enumerate(vals,1):ax.scatter(i+rng.normal(0,.035,len(v)),v,s=16,alpha=.65)
s=pd.read_csv(TAB/'EXTERNAL_VALIDATION_SUMMARY.csv').iloc[0];ax.set(ylabel='Terminal mass',title=f'External cohort (AUC={s.AUC:.3f})');ax.grid(axis='y',alpha=.2);save(fig,'07_external_terminal_mass.png')

# 8 blood coefficient forest
x=pd.read_csv(TAB/'BAYESIAN_TERMINAL_BLOOD_ASSOCIATIONS.csv');x['label']='T'+x.terminal_state.astype(str)+' · '+x.feature;xx=np.arange(len(x));fig,ax=plt.subplots(figsize=(7.2,max(4.5,len(x)*.26)));ax.errorbar(x.posterior_median,xx,xerr=[x.posterior_median-x.CrI_2_5,x.CrI_97_5-x.posterior_median],fmt='o',capsize=2);ax.axvline(0,linewidth=.8);ax.set_yticks(xx,x.label);ax.invert_yaxis();ax.set(xlabel='Standardized coefficient (95% interval)',title='Terminal-specific blood associations');ax.grid(axis='x',alpha=.2);save(fig,'08_terminal_blood_associations.png')

# 9 ME coefficient forest
x=pd.read_csv(TAB/'BAYESIAN_TERMINAL_ME_ASSOCIATIONS.csv');x['label']='T'+x.terminal_state.astype(str)+' · '+x.feature;xx=np.arange(len(x));fig,ax=plt.subplots(figsize=(7.2,max(4.5,len(x)*.24)));ax.errorbar(x.posterior_median,xx,xerr=[x.posterior_median-x.CrI_2_5,x.CrI_97_5-x.posterior_median],fmt='o',capsize=2);ax.axvline(0,linewidth=.8);ax.set_yticks(xx,x.label);ax.invert_yaxis();ax.set(xlabel='Adjusted standardized coefficient (95% interval)',title='Terminal-specific ME associations');ax.grid(axis='x',alpha=.2);save(fig,'09_terminal_ME_associations.png')

# 10 state context heatmap
x=pd.read_csv(TAB/'SOURCE_STATE_CONTEXT_SUMMARY.csv');features=[c for c in ['median_residual_time','LUAD_fraction','LUSC_fraction','Mye','Lym','CD8tot','Btot','Treg']+[f'ME{i}' for i in range(7)] if c in x.columns];A=x[features].to_numpy(float);mu=np.nanmean(A,axis=0);sd=np.nanstd(A,axis=0);Az=(A-mu)/(sd+1e-12);fig,ax=plt.subplots(figsize=(8.2,5.8));im=ax.imshow(Az,aspect='auto');ax.set_xticks(np.arange(len(features)),features,rotation=60,ha='right');ax.set_yticks(np.arange(len(x)),['S'+str(s) for s in x.state]);ax.set_title('State context fingerprints (feature-wise z scores)');fig.colorbar(im,ax=ax,label='z');save(fig,'10_state_context_fingerprints.png')

# 11 route-aware sensitivity
x=pd.read_csv(TAB/'ROUTE_AWARE_TEMPORAL_COMPARISON.csv');fig,ax=plt.subplots(figsize=(6.5,4.2));q=x.pivot(index='score',columns='dataset',values='spearman_rho');q.plot(kind='bar',ax=ax);ax.set(ylabel='Spearman ρ vs diagnosis proximity',xlabel='',title='Route-aware temporal score sensitivity');ax.axhline(0,linewidth=.8);ax.legend(frameon=False);ax.grid(axis='y',alpha=.2);save(fig,'11_route_aware_temporal_sensitivity.png')
# 12 route-stratified terminal-mass horizon summary
p=TAB/'ROUTE_STRATIFIED_TERMINAL_MASS.csv'
if p.exists():
    x=pd.read_csv(p)
    labs=x.route_group.tolist(); xpos=np.arange(len(x)); w=.24
    fig,ax=plt.subplots(figsize=(7.2,4.5))
    ax.bar(xpos-w,x.long_far_median,width=w,label='Longitudinal far')
    ax.bar(xpos,x.long_near_median,width=w,label='Longitudinal near')
    ax.bar(xpos+w,x.long_case_median,width=w,label='Diagnosis')
    ax.set_xticks(xpos,labs,rotation=12,ha='right');ax.set_ylabel('Median route-stratified terminal mass');ax.set_title('Route-stratified longitudinal terminal occupancy');ax.legend(frameon=False);ax.grid(axis='y',alpha=.2)
    save(fig,'12_route_stratified_terminal_mass.png')

print('rendered',len(list(OUT.glob('*.png'))),'PNG figures')
