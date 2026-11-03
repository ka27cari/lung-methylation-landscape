from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from .nmf import rank_consensus, cluster_components
from .utils import l2_rows

def _loading_stats(h):
    v=np.asarray(h,float); v=np.maximum(v,0); p=v/(v.sum()+1e-12); order=np.sort(p)[::-1]; cs=np.cumsum(order)
    n90=int(np.searchsorted(cs,.90)+1); n50=int(np.searchsorted(cs,.50)+1); entropy=float(-np.sum(p[p>0]*np.log(p[p>0]))); eff=float(np.exp(entropy))
    return n50,n90,eff,float(order[0] if len(order) else 0),float(np.sum(order[:min(50,len(order))]))

def discover_context(context, X, ranks, seeds, max_iter, representation, strategy, cpg):
    idx=np.asarray(context['indices'],int); Xi=X[idx]
    weights=context.get('weights')
    if weights is not None: weights=np.asarray(weights)[idx]
    comps=[]; loads=[]
    for r in ranks:
        if r>=min(Xi.shape): continue
        W,H,stab,err=rank_consensus(Xi,r,seeds,max_iter,weights)
        for j in range(r):
            comps.append(dict(strategy=strategy,representation=representation,context=context['name'],context_kind=context['kind'],
                              context_n=len(idx),rank=int(r),component=int(j),seed_stability=float(stab[j]),reconstruction_err=err))
            loads.append(H[j])
    return pd.DataFrame(comps), np.asarray(loads,np.float32)

def recurrent_families(comp, H, thresholds, min_distinct_ranks=2):
    fam_rows=[]; member_rows=[]; family_loads=[]
    for thr in thresholds:
        lab=cluster_components(comp,H,thr)
        for cl in sorted(set(lab)):
            ii=np.where(lab==cl)[0]; g=comp.iloc[ii]
            nr=g['rank'].nunique()
            if nr<min_distinct_ranks: continue
            h=l2_rows(np.mean(H[ii],axis=0,keepdims=True))[0]
            fid=f"{g.iloc[0].strategy}__{g.iloc[0].representation}__{g.iloc[0].context}__cos{int(round(100*thr))}__F{cl:03d}"
            n50,n90,eff,top1,top50=_loading_stats(h)
            fam_rows.append(dict(family_id=fid,cosine_threshold=thr,strategy=g.iloc[0].strategy,
                                 representation=g.iloc[0].representation,context=g.iloc[0].context,
                                 context_kind=g.iloc[0].context_kind,context_n=int(g.iloc[0].context_n),
                                 n_components=len(ii),n_ranks=int(nr),min_rank=int(g['rank'].min()),max_rank=int(g['rank'].max()),
                                 mean_rank=float(g['rank'].mean()),rank_span=int(g['rank'].max()-g['rank'].min()),
                                 low_rank_present=bool((g['rank']<=6).any()),mid_rank_present=bool(((g['rank']>=8)&(g['rank']<=12)).any()),
                                 high_rank_present=bool((g['rank']>=16).any()),seed_stability_mean=float(g['seed_stability'].mean()),
                                 loading_n50=int(n50),loading_n90=int(n90),loading_effective_cpg=float(eff),loading_top1_mass=top1,loading_top50_mass=top50))
            family_loads.append(h)
            for k in ii:
                member_rows.append(dict(family_id=fid,component_row=int(k),rank=int(comp.iloc[k]['rank']),component=int(comp.iloc[k]['component'])))
    return pd.DataFrame(fam_rows),pd.DataFrame(member_rows),np.asarray(family_loads,np.float32)

