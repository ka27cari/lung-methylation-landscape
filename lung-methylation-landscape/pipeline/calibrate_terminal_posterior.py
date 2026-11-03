#!/usr/bin/env python3
"""Source-only posterior-temperature calibration for the frozen clean consensus.

Tau is selected without residual-time, longitudinal, external, blood, or ME
information.  It minimizes leave-one-out binary negative log-likelihood for
hard terminal-vs-nonterminal membership.  For each source sample, only its own
assigned state centroid is recomputed leave-one-out; other centroids and the
state scales remain frozen from the consensus.
"""
import sys, json, yaml
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT=Path(__file__).resolve().parents[1] / "rebuild"
sys.path.insert(0,str(ROOT.parent/'src'))
from pfl.io_source import load_source

CFG=yaml.safe_load((ROOT.parent/'config/model.yaml').read_text())
MOD=ROOT/'cache/model'; TAB=ROOT/'analysis/tables'; EXP=ROOT/'cache/posterior_draws'
score=pd.read_csv(MOD/'SOURCE_PROGRAM_SCORES.csv')
Pcols=[c for c in score if c.startswith('program_')]
Z=score[Pcols].to_numpy(float)
X,ids,cpg,meta=load_source(ROOT/'data')
z=np.load(EXP/'ADMISSIBLE_COASSIGNMENT.npz',allow_pickle=True)
lab=z['labels'].astype(int); cent=z['centroids'].astype(float); terms=z['terminal_states'].astype(int).tolist(); K=int(z['K'][0])
ctrl=meta.group.eq('control').to_numpy(); pred=meta.group.eq('prediagnostic').to_numpy(); tt=pd.to_numeric(meta.time,errors='coerce').to_numpy(float)


def scales(X,l,c,f=.15):
    d=cdist(X,c,'sqeuclidean')
    s=np.array([np.median(d[l==j,j]) for j in range(len(c))])
    g=s[np.isfinite(s)&(s>0)]; floor=max(1e-4,np.median(g)*f)
    return np.where(np.isfinite(s)&(s>floor),s,floor)


def post_from_dist(d,s,tau):
    a=-d/(2*tau*s[None,:]); a-=a.max(1,keepdims=True); e=np.exp(a); return e/e.sum(1,keepdims=True)

s=scales(Z,lab,cent)
dfull=cdist(Z,cent,'sqeuclidean')
# Leave-one-out distance to assigned state centroid; distances to other states remain fixed.
dloo=dfull.copy()
for st in range(K):
    idx=np.where(lab==st)[0]
    if len(idx)<=1: continue
    total=Z[idx].sum(axis=0)
    for i in idx:
        c_loo=(total-Z[i])/(len(idx)-1)
        dloo[i,st]=float(np.sum((Z[i]-c_loo)**2))

hard_terminal=np.isin(lab,terms).astype(float)
rows=[]
for tau in [float(x) for x in CFG['posterior']['tau_grid']]:
    Ploo=post_from_dist(dloo,s,tau); tmloo=Ploo[:,terms].sum(1)
    p=np.clip(tmloo,1e-12,1-1e-12)
    loss=-(hard_terminal*np.log(p)+(1-hard_terminal)*np.log(1-p))
    # Full-data posterior is evaluated only for descriptive calibration columns.
    Pfull=post_from_dist(dfull,s,tau); tm=Pfull[:,terms].sum(1)
    rows.append(dict(tau=tau,loo_binary_nll=float(loss.mean()),loo_binary_nll_se=float(loss.std(ddof=1)/np.sqrt(len(loss))),
                     loo_terminal_nll=float(loss[hard_terminal==1].mean()),loo_nonterminal_nll=float(loss[hard_terminal==0].mean()),
                     loo_assigned_state_posterior_median=float(np.median(Ploo[np.arange(len(Ploo)),lab])),
                     source_control_median=float(np.median(tm[ctrl])),
                     source_far_median=float(np.median(tm[pred&(tt>=5)])),
                     source_near_median=float(np.median(tm[pred&(tt<3)]))))
R=pd.DataFrame(rows)
ix=int(R.loo_binary_nll.idxmin()); tau=float(R.loc[ix,'tau']); R['selected']=False; R.loc[ix,'selected']=True
R.to_csv(TAB/'POSTERIOR_TAU_CALIBRATION.csv',index=False)

P=post_from_dist(dfull,s,tau); tm=P[:,terms].sum(1)
out=pd.DataFrame({'sentrix_id':ids,'group':meta.group.astype(str),'time':pd.to_numeric(meta.time,errors='coerce'),'type':meta.type,'hard_state':lab,'TM':tm})
for k in range(K): out[f'P_state_{k:02d}']=P[:,k]
out.to_csv(MOD/'SOURCE_POSTERIOR.csv',index=False)
ld=np.load(MOD/'SELECTED_PROGRAM_LOADINGS.npz')
np.savez_compressed(MOD/'MODEL_RUNTIME_OUTPUT.npz',centroids=cent,state_centers=cent,state_scales=s,tau=np.array([tau]),
                    terminal_states=np.array(terms),hard_state=lab,posterior=P,TM=tm,
                    selected_loadings=ld['loadings'],program_center=ld['center'],program_scale=ld['scale'],
                    scaler_center=ld['center'],scaler_scale=ld['scale'])

basin=json.load(open(TAB/'TERMINAL_COHERENT_BASIN_POLICY.json'))
prog_policy=json.load(open(TAB/'PROGRAM_SELECTION_POLICY.json'))
spec=dict(K=K,n_programs=len(Pcols),stable_program_dictionary_size=int(prog_policy['stable_dictionary_size']),
          compact_landscape_basis_size=int(prog_policy['compact_landscape_basis_size']),
          program_basis_interpretation='all 14 stable recurrent families are retained as the reproducible program dictionary; the 9-program primary set is the parsimonious landscape basis, not a truth filter',
          terminal_states=terms,tau=tau,
          terminal_rule='source occupancy gain 95% CrI > 0 + destination-only in primary exact 95% partial OT + support/purity',
          consensus_rule='terminal-coherent basin selected without residual time; ordinary co-assignment consensus within the basin',
          selected_starts=basin['selected_starts'],selected_hierarchical_cut=basin['selected_cut'],
          tau_rule='minimize source leave-one-out binary NLL for frozen hard terminal vs non-terminal membership',
          residual_time_used_in_model_selection=False,longitudinal_used_in_model_selection=False,external_used_in_model_selection=False,
          blood_used_in_model_selection=False,ME_used_in_model_selection=False)
json.dump(spec,open(TAB/'FINAL_MODEL_SPEC.json','w'),indent=2)
pd.DataFrame([spec]).to_csv(TAB/'MODEL_SELECTION_SUMMARY.csv',index=False)

print('selected tau',tau)
print(R.loc[max(ix-4,0):min(ix+4,len(R)-1)].to_string(index=False))
print('source medians',float(np.median(tm[ctrl])),float(np.median(tm[pred&(tt>=5)])),float(np.median(tm[pred&(tt<3)])))
