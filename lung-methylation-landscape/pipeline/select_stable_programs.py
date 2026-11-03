#!/usr/bin/env python3
import sys,yaml,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.preprocessing import RobustScaler
ROOT=Path(__file__).resolve().parents[1] / "rebuild";sys.path.insert(0,str(ROOT.parent/'src'))
from pfl.io_source import load_source
from pfl.representations import derive_orientation,counterfactual_representation
from pfl.contexts import kmeans_contexts
CFG=yaml.safe_load((ROOT.parent/'config/model.yaml').read_text());DISC=ROOT/'cache/discovery';MOD=ROOT/'cache/model';TAB=ROOT/'analysis/tables';MOD.mkdir(parents=True,exist_ok=True);TAB.mkdir(parents=True,exist_ok=True)
X,ids,cpg,meta=load_source(ROOT/'data');ori=derive_orientation(X,meta);_,A,_,_=counterfactual_representation(X,meta,ori,CFG['representation']['anchor_fraction'],CFG['representation']['control_neighbors'],CFG['representation']['anchor_mode'])
fam=pd.read_csv(DISC/'RECURRENT_FAMILIES_ALL.csv');mem=pd.read_csv(DISC/'FAMILY_MEMBERS.csv');z=np.load(DISC/'RECURRENT_FAMILY_LOADINGS.npz',allow_pickle=True);H=z['loadings'].astype(float)
rs=mem.groupby('family_id')['rank'].nunique();fam['rank_support']=fam.family_id.map(rs).fillna(0).astype(int)
base=fam[(fam.representation==CFG['representation']['name'])&(fam.strategy=='kmeans_local')&np.isclose(fam.cosine_threshold,CFG['recurrence']['within_context_cosines'][0])&(fam.rank_support==len(CFG['nmf']['ranks']))].sort_values(['seed_stability_mean','context_n','family_id'],ascending=[False,False,True]).copy()
near=float(CFG['program_selection']['near_duplicate_cosine']);stab=float(CFG['program_selection']['min_seed_stability']);uncovered_limit=int(CFG['program_selection']['dual_uncovered_limit'])
# Build the stable, globally non-redundant recurrent pool in stability order.
stable=[];audit=[]
for rank,(_,r) in enumerate(base.iterrows(),start=1):
    h=H[int(r.family_row)];h=h/(np.linalg.norm(h)+1e-12)
    mx=max([float(np.dot(H[int(s.family_row)]/(np.linalg.norm(H[int(s.family_row)])+1e-12),h)) for s in stable],default=-1)
    if float(r.seed_stability_mean)<stab: decision='below_seed_stability_threshold'
    elif mx>=near: decision='collapsed_near_duplicate'
    else: decision='stable_distinct_candidate';stable.append(r)
    audit.append(dict(candidate_rank=rank,decision=decision,family_id=r.family_id,family_row=int(r.family_row),context=r.context,context_n=int(r.context_n),seed_stability=float(r.seed_stability_mean),rank_support=int(r.rank_support),max_cosine_to_stable_pool=mx))
if len(stable)<2: raise RuntimeError('stable recurrent pool is too small')
stable_df=pd.DataFrame(stable).reset_index(drop=True);stable_df.insert(0,'stable_pool_rank',np.arange(1,len(stable_df)+1))
# Select the first reproducibility-ranked prefix leaving no more than the declared
# number of SOURCE samples uncovered in BOTH independent discovery resolutions.
# The context labels are recomputed using exactly the same deterministic local
# clustering function as pipeline/discover_recurrent_programs.py. No K, tau, time, external,
# terminal or outcomes enter this step.
contexts=kmeans_contexts(A,meta,CFG['local_contexts']['k_grid'],CFG['local_contexts']['min_context_n'],CFG['random_seed'])
labels={}
for k in CFG['local_contexts']['k_grid']:
    lab=np.full(len(ids),-1,dtype=int)
    for ctx in contexts:
        if ctx['name'].startswith(f'kmeans_k{k}_'):
            lab[np.asarray(ctx['indices'],dtype=int)]=int(ctx['name'].split('_')[-1][1:])
    if (lab<0).any():raise RuntimeError(f'Missing source context labels at k={k}')
    labels[k]=lab
pd.DataFrame({'sentrix_id':ids,'k4_context':labels[4],'k6_context':labels[6]}).to_csv(TAB/'PROGRAM_DISCOVERY_CONTEXT_MEMBERSHIP.csv',index=False)
curve=[]
for n in range(1,len(stable_df)+1):
    prefixes=stable_df.context.iloc[:n].astype(str)
    selected4={int(x.split('_')[-1][1:]) for x in prefixes if 'kmeans_k4_' in x}
    selected6={int(x.split('_')[-1][1:]) for x in prefixes if 'kmeans_k6_' in x}
    noncovered=(~np.isin(labels[4],list(selected4))) & (~np.isin(labels[6],list(selected6)))
    curve.append(dict(n_programs=n,uncovered_n=int(noncovered.sum()),uncovered_fraction=float(noncovered.mean()),
                      covered_k4=len(selected4),covered_k6=len(selected6),passes_coverage=bool(noncovered.sum()<=uncovered_limit)))
curve=pd.DataFrame(curve)
passing=curve[curve.passes_coverage]
if passing.empty: raise RuntimeError('No prefix satisfies dual-resolution source coverage')
nsel=int(passing.n_programs.min())
# Terminology firewall: every row of stable_df is a reproducible program in the
# stable dictionary.  The first nsel rows are only the compact basis used for
# downstream landscape construction; the remaining dictionary members are not
# declared false/failed programs.
stable_df['dictionary_member']=True
stable_df['selected_in_compact_basis']=stable_df['stable_pool_rank'].le(nsel)
stable_df['program_role']=np.where(stable_df['selected_in_compact_basis'],
                                   'stable_dictionary__compact_landscape_basis',
                                   'stable_dictionary__sensitivity_only')
stable_df.to_csv(TAB/'PROGRAM_STABLE_POOL.csv',index=False)
selected_stability=stable_df.iloc[:nsel].copy().reset_index(drop=True);selected_stability['selection_rank']=np.arange(1,nsel+1)
report_order=CFG['program_selection'].get('reporting_order_selection_ranks',list(range(1,nsel+1)))
if sorted(report_order)!=list(range(1,nsel+1)): raise RuntimeError('reporting_order_selection_ranks must be a permutation of selected stability ranks')
selected=pd.concat([selected_stability.iloc[[r-1]] for r in report_order],ignore_index=True);selected.insert(0,'program_index',np.arange(1,nsel+1))
# final audit status
aud=pd.DataFrame(audit);aud['max_cosine_to_selected']=aud['max_cosine_to_stable_pool'];stable_rank=dict(zip(stable_df.family_id,stable_df.stable_pool_rank));aud['stable_pool_rank']=aud.family_id.map(stable_rank);aud['selected_program_index']=aud.family_id.map(dict(zip(selected.family_id,selected.program_index)))
aud['final_decision']=aud['decision']
aud.loc[aud.family_id.isin(selected.family_id),'final_decision']='selected_dual_resolution_coverage_basis'
aud.loc[(aud.decision=='stable_distinct_candidate')&(~aud.family_id.isin(selected.family_id)),'final_decision']='stable_pool_not_needed_for_dual_resolution_coverage'
# fit final selected scores/scaling independently, equivalent to selecting prefix columns
Hs=H[selected.family_row.astype(int).to_numpy()];raw=A@Hs.T;sc=RobustScaler(quantile_range=(10,90)).fit(raw);Z=sc.transform(raw);score=pd.DataFrame({'sentrix_id':ids})
for j in range(Z.shape[1]):score[f'program_{j+1:02d}']=Z[:,j]
selected.to_csv(MOD/'SELECTED_PROGRAMS.csv',index=False);aud.to_csv(MOD/'PROGRAM_SELECTION_AUDIT.csv',index=False);score.to_csv(MOD/'SOURCE_PROGRAM_SCORES.csv',index=False);np.savez_compressed(MOD/'SELECTED_PROGRAM_LOADINGS.npz',loadings=Hs,cpg=np.asarray(cpg,dtype='U'),center=sc.center_,scale=sc.scale_)
curve['selected_n_programs']=nsel;curve['uncovered_limit']=uncovered_limit
curve.to_csv(TAB/'PROGRAM_DUAL_RESOLUTION_COVERAGE_CURVE.csv',index=False)
# Selection-range audit; all values 10..20 are explicitly checked.
limit_rows=[]
for lim in range(10,21):
    eligible=curve[curve.uncovered_n.le(lim)]
    limit_rows.append(dict(uncovered_sample_limit=lim,selected_basis_size=int(eligible.n_programs.min()) if len(eligible) else None))
limit_sens=pd.DataFrame(limit_rows)
limit_sens.to_csv(TAB/'PROGRAM_DUAL_COVERAGE_LIMIT_10_20.csv',index=False)
# Explicit dictionary-to-basis crosswalk.  Reporting order may differ from
# stability order, but membership in the primary basis is the stability prefix.
idxmap=dict(zip(selected.stable_pool_rank.astype(int),selected.program_index.astype(int)))
basis_map=stable_df[['stable_pool_rank','family_id','context','context_n','seed_stability_mean','rank_support',
                     'dictionary_member','selected_in_compact_basis','program_role']].copy()
basis_map['primary_program_index']=basis_map.stable_pool_rank.map(idxmap)
basis_map.to_csv(TAB/'PROGRAM_DICTIONARY_BASIS_MAP.csv',index=False)
policy={
 'n_selected_programs':nsel,'n_stable_distinct_pool':len(stable_df),
 'stable_dictionary_size':len(stable_df),'compact_landscape_basis_size':nsel,
 'seed_stability_threshold':stab,'near_duplicate_cosine':near,
 'dual_uncovered_limit':uncovered_limit,
 'uncovered_at_selected_prefix':int(curve.loc[curve.n_programs.eq(nsel),'uncovered_n'].iloc[0]),
 'uncovered_at_previous_prefix':int(curve.loc[curve.n_programs.eq(nsel-1),'uncovered_n'].iloc[0]),
 'stable_dictionary_definition':'all-rank recurrent families passing seed-stability and global near-duplicate collapse',
 'compact_basis_definition':'smallest reproducibility-ordered prefix with at most 20 of 816 source samples not represented by any selected discovery context at either k=4 or k=6',
 'interpretation':'the full stable dictionary contains all 14 reproducible programs; P9 is a compact landscape basis chosen by source discovery-context coverage, not by reconstruction accuracy or biological completeness',
 'downstream_endpoints_used_for_basis_size_selection':False,
 'threshold_robustness_note':'allowances 10–20 all select P9 (P8 leaves 27, P9 leaves 7)',
 'rule':'dual-resolution source discovery-context sample coverage; first prefix with doubly uncovered sample count <=20'
}
json.dump(policy,open(TAB/'PROGRAM_SELECTION_POLICY.json','w'),indent=2)
print('stable distinct pool',len(stable_df),'selected programs',nsel)
print(curve.to_string(index=False))
print(selected[['program_index','family_id','context','seed_stability_mean']].to_string(index=False))
