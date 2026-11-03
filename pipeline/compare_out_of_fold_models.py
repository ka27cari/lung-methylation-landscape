#!/usr/bin/env python3
"""Repeated paired out-of-fold AUC table used by r/figures/terminal_blood_associations.R.

This is descriptive model-comparison support, not part of landscape selection.
Matched control/prediagnostic pairs are never split across train/test folds.
All standardization, one-hot encoding, and blood PCA are fitted within each training fold.
"""
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.impute import SimpleImputer

ROOT=Path(__file__).resolve().parents[1] / "rebuild"
TAB=ROOT/'analysis/tables'; MOD=ROOT/'cache/model'
meta=pd.read_csv(ROOT/'data/nested.csv')
scores=pd.read_csv(MOD/'SOURCE_PROGRAM_SCORES.csv')
post=pd.read_csv(MOD/'SOURCE_POSTERIOR.csv')
ids=scores.sentrix_id.astype(str)
meta=meta.set_index(meta.sentrix_id.astype(str),drop=False).loc[ids].reset_index(drop=True)
post=post.set_index(post.sentrix_id.astype(str),drop=False).loc[ids].reset_index(drop=True)
assert (meta.sentrix_id.astype(str).to_numpy()==ids.to_numpy()).all()

pcols=[c for c in scores.columns if c.startswith('program_')]
mecols=[f'ME{i}' for i in range(7)]
cell=['Bas','Bmem','Bnv','CD4mem','CD4nv','CD8mem','CD8nv','Eos','Mono','Neu','NK','Treg']
# Per-sample CLR has no population fitting. Fixed source-wide pseudocount is label-free.
A=meta[cell].apply(pd.to_numeric,errors='coerce').to_numpy(float)
med=np.nanmedian(A,axis=0); bad=~np.isfinite(A); A[bad]=np.take(med,np.where(bad)[1])
pos=A[A>0]; pc=float(np.min(pos)/2) if len(pos) else 1e-6
A=np.where(A>0,A,pc); B=np.log(A); B-=B.mean(axis=1,keepdims=True)
for j,c in enumerate(cell): meta['CLR_'+c]=B[:,j]
for c in pcols: meta[c]=scores[c].to_numpy(float)
meta['TM']=post.TM.to_numpy(float)
meta['y']=meta.group.astype(str).str.lower().eq('prediagnostic').astype(int)
meta['pair_group']=meta.origin.astype(str)+'::'+meta.matched_index.astype(str)

num_epi=['age','smk_amnt']; cat_epi=['sex','smk_status','origin']; blood=['CLR_'+c for c in cell]

def make_pre(features):
    trs=[]
    num=[c for c in features if c in num_epi+mecols+pcols+['TM']]
    cat=[c for c in features if c in cat_epi]
    has_blood=any(c in features for c in blood)
    if num:
        trs.append(('num',Pipeline([('imp',SimpleImputer(strategy='median')),('scale',StandardScaler())]),num))
    if cat:
        trs.append(('cat',Pipeline([('imp',SimpleImputer(strategy='most_frequent')),('oh',OneHotEncoder(handle_unknown='ignore'))]),cat))
    if has_blood:
        trs.append(('blood',Pipeline([('imp',SimpleImputer(strategy='median')),('scale',StandardScaler()),('pca',PCA(n_components=3,random_state=1))]),blood))
    return ColumnTransformer(trs,remainder='drop')

models={
 'Epidemiology': num_epi+cat_epi,
 'Epidemiology+Blood': num_epi+cat_epi+blood,
 'Epidemiology+Blood+ME': num_epi+cat_epi+blood+mecols,
 'Epidemiology+Blood+ME+TM': num_epi+cat_epi+blood+mecols+['TM'],
 'Epidemiology+Blood+ME+Programs': num_epi+cat_epi+blood+mecols+pcols,
 'Programs': pcols,
 'TM': ['TM'],
}
X=meta; y=meta.y.to_numpy(int); groups=meta.pair_group.to_numpy(str)
unique=np.unique(groups); rng=np.random.default_rng(20261004)
rows=[]; nrep=20; nfold=5
for rep in range(nrep):
    ug=unique.copy(); rng.shuffle(ug)
    fmap={g:i%nfold for i,g in enumerate(ug)}
    fold=np.array([fmap[g] for g in groups])
    for name,features in models.items():
        pred=np.full(len(X),np.nan)
        for f in range(nfold):
            tr=fold!=f; te=fold==f
            pipe=Pipeline([('pre',make_pre(features)),('logit',LogisticRegression(max_iter=3000,C=1.0,solver='liblinear',random_state=rep))])
            pipe.fit(X.loc[tr,features],y[tr])
            pred[te]=pipe.predict_proba(X.loc[te,features])[:,1]
        rows.append(dict(repeat=rep+1,model=name,AUC=roc_auc_score(y,pred)))
raw=pd.DataFrame(rows); raw.to_csv(TAB/'INCREMENTAL_INFORMATION_OOF_REPEATS.csv',index=False)
sumdf=raw.groupby('model').AUC.agg(['mean','std']).reset_index().rename(columns={'mean':'AUC_mean','std':'AUC_sd'})
sumdf['n_repeats']=nrep
sumdf['CV_unit']='matched_pair'; sumdf['n_folds']=nfold
sumdf.to_csv(TAB/'INCREMENTAL_INFORMATION_OOF_SUMMARY.csv',index=False)
print(sumdf.to_string(index=False))
