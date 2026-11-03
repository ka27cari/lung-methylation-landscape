#!/usr/bin/env python3
"""Fit the source landscape and form the terminal-coherent consensus.

Selection is source-only. Residual time, longitudinal/external outcomes, blood,
ME, histology, requested terminal count, ARI, and any named historical start do
not participate in basin selection.

Pipeline
--------
1. Fit deterministic K-means starts over the source program scores.
2. Define per-start terminal candidates by prediagnostic occupancy gain plus
   destination-only status in the PRIMARY exact 95% partial-OT plan, support,
   and control-contamination limits. 90/99% are sensitivity only.
3. Select K as the finest resolution whose median smallest terminal support
   across starts remains >= the declared support threshold.
4. At selected K, compare source-admissible starts by label-free terminal-set
   overlap (Hungarian-matched Jaccard), cluster starts by average linkage, and
   scan hierarchical cuts from coarse to fine.
5. A basin is terminal-coherent when its ordinary co-assignment consensus is
   source-admissible and every consensus terminal is matched at Jaccard > 0.5
   in >50% of member starts. Select the first/coarsest cut containing any such
   basin; if more than one exists, choose largest then strongest support.
6. Form the final ordinary co-assignment consensus only from that basin.
"""
import os, sys, json, yaml
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='1'
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.spatial.distance import cdist
from scipy.optimize import linprog, linear_sum_assignment

ROOT=Path(__file__).resolve().parents[1] / "rebuild"
sys.path.insert(0,str(ROOT.parent/'src'))
from pfl.io_source import load_source
from pfl.representations import derive_orientation, counterfactual_representation

CFG=yaml.safe_load((ROOT.parent/'config/model.yaml').read_text())
MOD=ROOT/'cache/model'; TAB=ROOT/'analysis/tables'; EXP=ROOT/'cache/posterior_draws'
for p in [MOD,TAB,EXP]: p.mkdir(exist_ok=True,parents=True)

# Program scores are the canonical output of step 03. Reuse them here rather than
# recomputing the expensive counterfactual representation a second time.
score=pd.read_csv(MOD/'SOURCE_PROGRAM_SCORES.csv')
ids=score.sentrix_id.astype(str).to_numpy()
meta0=pd.read_csv(ROOT/'data/nested.csv'); meta0.sentrix_id=meta0.sentrix_id.astype(str)
meta=meta0.set_index('sentrix_id',drop=False).loc[ids].reset_index(drop=True)
Pcols=[c for c in score.columns if c.startswith('program_')]
Z=score[Pcols].to_numpy(float)
ctrl=meta.group.eq('control').to_numpy(); pred=meta.group.eq('prediagnostic').to_numpy()
Nc=int(ctrl.sum()); Np=int(pred.sum())

TDEF=CFG['terminal_definition']
PRIMARY=float(TDEF.get('primary_partial_ot_mass',0.95))
SENS=[float(x) for x in TDEF.get('partial_ot_masses',[0.90,0.95,0.99])]
JTHR=float(CFG['consensus'].get('terminal_identity_jaccard_threshold',0.5))
MAJ=float(CFG['consensus'].get('terminal_identity_majority_fraction',0.5))


def partial_plan(a,b,C,mass):
    K=len(a); Au=[]; bu=[]
    for i in range(K):
        q=np.zeros(K*K); q[i*K:(i+1)*K]=1; Au.append(q); bu.append(a[i])
    for j in range(K):
        q=np.zeros(K*K); q[j::K]=1; Au.append(q); bu.append(b[j])
    r=linprog(C.ravel(),A_ub=np.asarray(Au),b_ub=np.asarray(bu),
              A_eq=np.ones((1,K*K)),b_eq=[mass],bounds=(0,None),method='highs',
              options={'time_limit':5.0})
    if not r.success: raise RuntimeError(r.message)
    return r.x.reshape(K,K)


def occupancy_posterior(cc,pp,seed,draws):
    rng=np.random.default_rng(int(seed)); qlo=[]; qhi=[]; pg=[]
    alpha=(1-float(TDEF['occupancy_gain_credible_level']))/2
    for s in range(len(cc)):
        d=(rng.beta(pp[s]+.5,Np-pp[s]+.5,int(draws))-
           rng.beta(cc[s]+.5,Nc-cc[s]+.5,int(draws)))
        lo,hi=np.quantile(d,[alpha,1-alpha])
        qlo.append(float(lo)); qhi.append(float(hi)); pg.append(float(np.mean(d>0)))
    return np.asarray(qlo),np.asarray(qhi),np.asarray(pg)


def detect(lab,cent,draw_seed,draws):
    K=len(cent); cc=np.bincount(lab[ctrl],minlength=K); pp=np.bincount(lab[pred],minlength=K); sz=cc+pp
    q025,q975,pg=occupancy_posterior(cc,pp,draw_seed,draws)
    C=cdist(cent,cent,'sqeuclidean'); pos=C[C>0]; C=C/(np.median(pos) if len(pos) else 1.0)
    plans={m:partial_plan(cc/Nc,pp/Np,C,m) for m in sorted(set(SENS+[PRIMARY]))}
    export_by_mass={m:np.array([G[s].sum()-G[s,s] for s in range(K)]) for m,G in plans.items()}
    export_primary=export_by_mass[PRIMARY]
    destination_primary=export_primary<=1e-9
    terms=[s for s in range(K) if q025[s]>0 and destination_primary[s]]
    pn=[int(pp[t]) for t in terms]
    cf=[float(cc[t]/sz[t]) for t in terms]
    ok=(len(terms)>=int(TDEF['min_terminal_count']) and
        min(pn,default=0)>=int(TDEF['min_prediagnostic_n']) and
        max(cf,default=1)<=float(TDEF['max_control_fraction']))
    return dict(terms=terms,pred_n=pn,ctrl_frac=cf,ok=bool(ok),cc=cc,pp=pp,prob=pg,
                q025=q025,q975=q975,plans=plans,export_by_mass=export_by_mass,
                destination_primary=destination_primary,C=C,draw_seed=int(draw_seed),draws=int(draws))


def fit_one(K,st):
    seed=int(CFG['random_seed']+K*1009+st*7919)
    km=KMeans(K,n_init=1,max_iter=500,random_state=seed).fit(Z)
    d=detect(km.labels_,km.cluster_centers_,K*100000+st,int(TDEF['start_occupancy_draws']))
    row=dict(K=K,start=st,seed=seed,occupancy_draw_seed=d['draw_seed'],passes=d['ok'],
             n_terminal=len(d['terms']),terminal_states=';'.join(map(str,d['terms'])),
             min_pred=min(d['pred_n'],default=0),max_ctrl=max(d['ctrl_frac'],default=1),
             primary_partial_ot_mass=PRIMARY)
    return row,km.labels_.astype(np.int16),km.cluster_centers_.astype(np.float32),d


def consensus_from_labels(labels,K):
    S=np.zeros((len(Z),len(Z)),dtype=float)
    for l in labels: S+=(l[:,None]==l[None,:]).astype(float)
    S/=len(labels)
    cl=AgglomerativeClustering(n_clusters=K,metric='precomputed',linkage=CFG['consensus']['linkage']).fit_predict(1-S)
    cent=np.vstack([Z[cl==j].mean(0) for j in range(K)])
    d=detect(cl,cent,int(TDEF['consensus_occupancy_draw_seed']),int(TDEF['consensus_occupancy_draws']))
    return dict(S=S,lab=cl,cent=cent,detect=d)


def jaccard_mask(a,b):
    inter=np.sum(a & b); union=np.sum(a | b)
    return float(inter/union) if union else 1.0


def terminal_set_similarity(lab_a,terms_a,lab_b,terms_b):
    A=[lab_a==t for t in terms_a]; B=[lab_b==t for t in terms_b]
    if not A or not B: return 0.0
    M=np.array([[jaccard_mask(a,b) for b in B] for a in A],float)
    rr,cc=linear_sum_assignment(-M)
    return float(M[rr,cc].sum()/max(len(A),len(B)))


def consensus_terminal_support(member_records,cons):
    cterms=cons['detect']['terms']; cv=[cons['lab']==t for t in cterms]
    A=np.zeros((len(member_records),len(cterms)),float)
    matched_state=np.full((len(member_records),len(cterms)),-1,int)
    for i,r in enumerate(member_records):
        sv=[r['lab']==t for t in r['detect']['terms']]
        if not sv: continue
        M=np.array([[jaccard_mask(a,b) for b in sv] for a in cv],float)
        n=max(len(cv),len(sv)); P=np.zeros((n,n)); P[:len(cv),:len(sv)]=M
        rr,cc=linear_sum_assignment(-P)
        for q,row in enumerate(rr):
            if row>=len(cv): continue
            col=cc[q]; A[i,row]=P[row,col]
            if col<len(r['detect']['terms']): matched_state[i,row]=int(r['detect']['terms'][col])
    frac=(A>JTHR).mean(axis=0) if A.size else np.array([])
    med=np.median(A,axis=0) if A.size else np.array([])
    return dict(matrix=A,matched_state=matched_state,fraction=frac,median=med,
                qualifies=bool(len(frac) and np.all(frac>MAJ)))

# -----------------------------------------------------------------------------
# Fit K grid and choose finest typically-supported resolution.
# -----------------------------------------------------------------------------
start_rows=[]; perK={}
for K in CFG['state_resolution']['K_grid']:
    # Keep HiGHS solves on the main thread. This is faster and more reliable than
    # invoking scipy.optimize.linprog repeatedly inside a Python thread pool.
    res=[fit_one(K,st) for st in range(CFG['state_resolution']['starts_per_K'])]
    res=sorted(res,key=lambda x:x[0]['start'])
    start_rows.extend([x[0] for x in res])
    valid_support=[x[0]['min_pred'] for x in res if x[0]['n_terminal']>=int(TDEF['min_terminal_count'])]
    typical=float(np.median(valid_support)) if valid_support else 0.0
    perK[K]=dict(res=res,n_admissible=sum(bool(x[0]['passes']) for x in res),typical=typical)
    pd.DataFrame(start_rows).to_csv(TAB/'K_START_ADMISSIBILITY_SCAN.csv',index=False)
    print(f'K={K}: admissible {perK[K]["n_admissible"]}/{len(res)}, median smallest terminal support={typical:.1f}',flush=True)

pd.DataFrame(start_rows).to_csv(TAB/'K_START_ADMISSIBILITY_SCAN.csv',index=False)
support_min=float(CFG['state_resolution']['typical_terminal_support_min'])
passing_K=[K for K,v in perK.items() if v['typical']>=support_min and v['n_admissible']>=2]
selected_K=max(passing_K)
json.dump(dict(selected_K=selected_K,
               rule='choose the finest resolution whose median smallest primary-LP95 terminal-candidate support across deterministic starts remains at least the declared source support threshold',
               typical_terminal_support_min=support_min,primary_partial_ot_mass=PRIMARY),
          open(TAB/'K_SELECTION_POLICY.json','w'),indent=2)

# all-admissible K audit is retained as a sensitivity summary; final selection below uses terminal-coherent basin.
crows=[]
for K,v in perK.items():
    adm=[x for x in v['res'] if x[0]['passes']]
    crows.append(dict(K=K,n_admissible=len(adm),admissible_fraction=len(adm)/CFG['state_resolution']['starts_per_K'],
                      median_smallest_terminal_pred_support=v['typical'],passes_typical_support=v['typical']>=support_min,
                      selected=K==selected_K))
Cdf=pd.DataFrame(crows); Cdf.to_csv(TAB/'K_SELECTION_STABILITY_AUDIT.csv',index=False); Cdf.to_csv(TAB/'K_CONSENSUS_AUDIT.csv',index=False)

# -----------------------------------------------------------------------------
# Clean terminal-coherent basin discovery at selected K. No residual time.
# -----------------------------------------------------------------------------
res=perK[selected_K]['res']
adm=[]
for row,lab0,cent0,d0 in res:
    if row['passes']:
        adm.append(dict(row=row,lab=lab0.astype(int),cent=cent0.astype(float),detect=d0))
if len(adm)<3: raise RuntimeError('Too few source-admissible starts for coherent-basin discovery')

pd.DataFrame([a['row'] for a in adm]).to_csv(TAB/'K16_ADMISSIBLE_PARTITION_INDEX.csv',index=False)
n=len(adm); Sim=np.eye(n,dtype=float)
for i in range(n):
    for j in range(i+1,n):
        Sim[i,j]=Sim[j,i]=terminal_set_similarity(adm[i]['lab'],adm[i]['detect']['terms'],adm[j]['lab'],adm[j]['detect']['terms'])
Dist=1-Sim

cut_audit=[]; chosen=None; chosen_cut=None
cut_min=int(CFG['consensus'].get('mode_cut_min',2)); cut_max=min(int(CFG['consensus'].get('mode_cut_max',12)),n)
for kcut in range(cut_min,cut_max+1):
    grp=AgglomerativeClustering(n_clusters=kcut,metric='precomputed',linkage='average').fit_predict(Dist)
    qualifying=[]
    for c in range(kcut):
        idx=np.where(grp==c)[0]
        if len(idx)<3: continue
        members=[adm[i] for i in idx]
        cons=consensus_from_labels([m['lab'] for m in members],selected_K)
        sup=consensus_terminal_support(members,cons)
        row=dict(hierarchical_cut_k=kcut,cluster=int(c),n_starts=len(members),
                 starts=';'.join(str(m['row']['start']) for m in members),
                 consensus_source_admissible=bool(cons['detect']['ok']),
                 consensus_n_terminal=len(cons['detect']['terms']),
                 consensus_terminal_states=';'.join(map(str,cons['detect']['terms'])),
                 min_terminal_member_support_fraction=float(np.min(sup['fraction'])) if len(sup['fraction']) else 0.0,
                 min_terminal_median_matched_jaccard=float(np.min(sup['median'])) if len(sup['median']) else 0.0,
                 terminal_coherent=bool(cons['detect']['ok'] and sup['qualifies']))
        cut_audit.append(row)
        if row['terminal_coherent']:
            qualifying.append((members,cons,sup,row))
    if qualifying and chosen is None:
        qualifying.sort(key=lambda x:(len(x[0]),float(np.min(x[2]['fraction'])),float(np.mean(x[2]['median']))),reverse=True)
        chosen=qualifying[0]; chosen_cut=kcut
        break
if chosen is None: raise RuntimeError('No terminal-coherent basin found in declared cut range')

members,cons,sup,chosen_row=chosen
pd.DataFrame(cut_audit).to_csv(TAB/'START_MODE_CUT_AUDIT.csv',index=False)
selected_starts=[int(m['row']['start']) for m in members]
pd.DataFrame([m['row'] for m in members]).to_csv(TAB/'SELECTED_TERMINAL_COHERENT_STARTS.csv',index=False)

support_rows=[]
for j,t in enumerate(cons['detect']['terms']):
    support_rows.append(dict(consensus_terminal_state=int(t),
                             fraction_member_starts_with_matched_Jaccard_gt_threshold=float(sup['fraction'][j]),
                             median_matched_jaccard=float(sup['median'][j]),
                             min_matched_jaccard=float(sup['matrix'][:,j].min()),
                             max_matched_jaccard=float(sup['matrix'][:,j].max()),
                             threshold=JTHR,majority_fraction=MAJ))
pd.DataFrame(support_rows).to_csv(TAB/'SELECTED_BASIN_TERMINAL_SUPPORT.csv',index=False)
json.dump(dict(selected_cut=int(chosen_cut),selected_starts=selected_starts,n_selected_starts=len(selected_starts),
               selection_rule=CFG['consensus']['basin_rule'],terminal_identity_jaccard_threshold=JTHR,
               terminal_identity_majority_fraction=MAJ,residual_time_used=False,longitudinal_used=False,
               external_used=False,blood_used=False,ME_used=False,histology_used=False,ARI_used=False,
               requested_terminal_count_used=False),open(TAB/'TERMINAL_COHERENT_BASIN_POLICY.json','w'),indent=2)

# Final consensus from selected basin.
S=cons['S']; lab=cons['lab']; cent=cons['cent']; d=cons['detect']
all_adm_labels=np.stack([a['lab'] for a in adm]); sel_labels=np.stack([m['lab'] for m in members])
np.savez_compressed(EXP/'ADMISSIBLE_COASSIGNMENT.npz',
                    coassignment=S,ids=np.asarray(ids,dtype='U'),labels=lab.astype(np.int16),centroids=cent.astype(np.float32),
                    terminal_states=np.array(d['terms'],dtype=np.int16),K=np.array([selected_K]),
                    admissible_labels=all_adm_labels.astype(np.int16),admissible_starts=np.array([a['row']['start'] for a in adm],dtype=np.int16),
                    selected_basin_labels=sel_labels.astype(np.int16),selected_basin_starts=np.array(selected_starts,dtype=np.int16),
                    selected_hierarchical_cut=np.array([chosen_cut],dtype=np.int16),
                    consensus_occupancy_q025=d['q025'],consensus_occupancy_q975=d['q975'],
                    consensus_occupancy_draw_seed=np.array([d['draw_seed']]),consensus_occupancy_draws=np.array([d['draws']]))

pd.DataFrame(cent,columns=[f'program_{i+1:02d}' for i in range(Z.shape[1])]).assign(state=np.arange(selected_K)).to_csv(MOD/'STATE_CENTROIDS.csv',index=False)
G95=d['plans'][PRIMARY]
pd.DataFrame(G95,index=[f'S{i}' for i in range(selected_K)],columns=[f'S{i}' for i in range(selected_K)]).to_csv(MOD/'OT_PLAN.csv')

rows=[]
for s in range(selected_K):
    nstate=int(d['cc'][s]+d['pp'][s]); cf=float(d['cc'][s]/nstate) if nstate else np.nan
    gain=float(d['pp'][s]/Np-d['cc'][s]/Nc)
    role='Terminal' if s in d['terms'] else ('Source' if gain<=0 else 'Intermediate')
    mm=pred&(lab==s); types=meta.loc[mm,'type'].astype(str)
    rows.append(dict(state=s,n=nstate,prediagnostic_n=int(d['pp'][s]),control_n=int(d['cc'][s]),control_fraction=cf,
                     time_median=float(pd.to_numeric(meta.loc[mm,'time'],errors='coerce').median()) if mm.sum() else np.nan,
                     LUAD_fraction=float((types=='LUAD').mean()) if len(types) else np.nan,
                     LUSC_fraction=float((types=='LUSC').mean()) if len(types) else np.nan,
                     control_occupancy=float(d['cc'][s]/Nc),prediagnostic_occupancy=float(d['pp'][s]/Np),
                     occupancy_gain_mean=gain,occupancy_gain_q025=float(d['q025'][s]),occupancy_gain_q975=float(d['q975'][s]),
                     posterior_prob_gain=float(d['prob'][s]),gain_95_supported=bool(d['q025'][s]>0),
                     destination_only_partial_OT=bool(d['destination_primary'][s]),primary_partial_ot_mass=PRIMARY,
                     final_microstate_terminal=bool(s in d['terms']),macrostate_role=role,
                     partial_cross_out_primary=float(d['export_by_mass'][PRIMARY][s]),
                     visualization_sparse=bool(nstate<5),include_in_primary_state_visualization=bool(nstate>=5 or s in d['terms'])))
R=pd.DataFrame(rows)
R.to_csv(TAB/'FINAL_STATE_ROLE_CLASSIFICATION.csv',index=False); R.to_csv(MOD/'FINAL_STATE_ROLE_CLASSIFICATION.csv',index=False)
pd.DataFrame({'state':R.state,'macrostate':R.macrostate_role,
              'definition_note':['95% occupancy-gain interval > 0 + destination-only in primary exact 95% partial OT' if x=='Terminal' else ('posterior mean occupancy change <= 0' if x=='Source' else 'all remaining states') for x in R.macrostate_role]}).to_csv(TAB/'MACROSTATE_MAPPING_FINAL.csv',index=False)
pd.DataFrame({'state':R.state,'n':R.n,'prediagnostic_n':R.prediagnostic_n,'macrostate_role':R.macrostate_role,
              'visualization_sparse':R.visualization_sparse,'include_in_primary_state_visualization':R.include_in_primary_state_visualization}).to_csv(TAB/'STATE_VISUALIZATION_INDEX.csv',index=False)
pd.DataFrame([
    dict(macrostate='Source',operational_rule='posterior mean prediagnostic-control occupancy change <= 0',transport_interpretation='depleting/source-like population sector'),
    dict(macrostate='Intermediate',operational_rule='states that are neither Source-depleting nor supported primary-LP95 destination terminal sinks',transport_interpretation='relay/uncertain population sector'),
    dict(macrostate='Terminal',operational_rule='95% occupancy-gain interval > 0 + destination-only in primary exact 95% partial OT',transport_interpretation='population-gaining destination sinks')
]).to_csv(TAB/'MACROSTATE_RATIONALE_AUDIT.csv',index=False)

# 90/95/99 sensitivity on the frozen final consensus; only 95 participates in discovery.
sens=[]
for mass,G in sorted(d['plans'].items()):
    export=np.array([G[s].sum()-G[s,s] for s in range(selected_K)])
    terms_m=[s for s in range(selected_K) if d['q025'][s]>0 and export[s]<=1e-9]
    pn=[int(d['pp'][t]) for t in terms_m]; cf=[float(d['cc'][t]/(d['cc'][t]+d['pp'][t])) for t in terms_m]
    sens.append(dict(partial_mass=mass,primary=bool(abs(mass-PRIMARY)<1e-12),n_terminal=len(terms_m),terminal_states=';'.join(map(str,terms_m)),
                     min_pred=min(pn,default=0),max_ctrl=max(cf,default=1)))
pd.DataFrame(sens).to_csv(TAB/'TERMINAL_PARTIAL_OT_SENSITIVITY.csv',index=False)

print('selected K',selected_K)
print('source-admissible starts',len(adm))
print('coarsest coherent cut',chosen_cut)
print('selected starts',selected_starts)
print('final terminals',d['terms'])
print(pd.DataFrame(support_rows).to_string(index=False))
