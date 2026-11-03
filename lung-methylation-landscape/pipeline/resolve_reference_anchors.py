#!/usr/bin/env python3
import os,sys,json,yaml
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler
from sklearn.neighbors import NearestNeighbors
ROOT=Path(__file__).resolve().parents[1] / "rebuild";sys.path.insert(0,str(ROOT.parent/'src'))
from pfl.io_source import load_source
from pfl.representations import derive_orientation,stable_anchor_indices,counterfactual_representation
CFG=yaml.safe_load((ROOT.parent/'config/model.yaml').read_text());DATA=ROOT/'data';OUT=ROOT/'analysis/tables';OUT.mkdir(parents=True,exist_ok=True)
X,ids,cpg,meta=load_source(DATA);ori=derive_orientation(X,meta);fracs=list(map(float,CFG['anchor_audit']['fractions']))

def mae(anchors,k=10):
 ctrl=np.where(meta.group.eq('control').to_numpy())[0];sc=RobustScaler(quantile_range=(10,90)).fit(X[ctrl][:,anchors]);C=sc.transform(X[ctrl][:,anchors]);nn=NearestNeighbors(n_neighbors=min(k+1,len(ctrl))).fit(C);loc=nn.kneighbors(C,return_distance=False);e=[]
 for ii,row in enumerate(loc):
  js=[j for j in row if j!=ii][:k];ref=X[ctrl[js]].mean(0);e.append(np.mean(np.abs(X[ctrl[ii]]-ref)))
 return float(np.mean(e))
def pr30(A):
 p=PCA(n_components=min(30,A.shape[0]-1,A.shape[1]),svd_solver='randomized',random_state=CFG['random_seed']).fit(A);v=p.explained_variance_;return float(v.sum()**2/(np.square(v).sum()+1e-12))
rows=[]; prev=None; reps={}
for f in fracs:
 anchors=stable_anchor_indices(X,meta,ori,f,CFG['representation']['anchor_mode']);_,A,_,_=counterfactual_representation(X,meta,ori,f,CFG['representation']['control_neighbors'],CFG['representation']['anchor_mode']);A=A.astype(np.float32);reps[f]=A
 corr=np.nan
 if prev is not None:
  a=reps[prev].ravel();b=A.ravel();corr=float(np.corrcoef(a,b)[0,1]);del reps[prev]
 rows.append(dict(anchor_fraction=f,anchor_n=len(anchors),control_reconstruction_MAE=mae(anchors,CFG['representation']['control_neighbors']),counterfactual_PCA30_participation_ratio=pr30(A),adjacent_representation_correlation=corr))
 prev=f
R=pd.DataFrame(rows)
# Selection: choose the prespecified anchor fraction with the minimum effective
# dimensionality of the counterfactual representation.  Participation ratio
# is computed from the first 30 source PCA eigenvalues.  Adjacent-representation
# correlation and control-reconstruction MAE are retained as sensitivity
# diagnostics only and do not enter the primary selection rule.
pr=R.counterfactual_PCA30_participation_ratio.to_numpy(float)
sel_idx=int(np.nanargmin(pr))
selected=float(R.loc[sel_idx,'anchor_fraction'])
R['selected']=R.anchor_fraction.eq(selected);R.to_csv(OUT/'ANCHOR_RESOLUTION_SOURCE_ONLY_AUDIT.csv',index=False)
json.dump(dict(
 selected_anchor_fraction=selected,
 rule='minimum counterfactual PCA30 participation ratio across the prespecified anchor-fraction grid',
 criterion='counterfactual_PCA30_participation_ratio',
 direction='minimize',
 selection_uses_adjacent_representation_correlation=False,
 selection_uses_control_reconstruction_MAE=False,
 prespecified_grid=[float(x) for x in fracs],
 selected_participation_ratio=float(R.loc[sel_idx,'counterfactual_PCA30_participation_ratio'])
),open(OUT/'ANCHOR_RESOLUTION_POLICY.json','w'),indent=2)
print(R.to_string(index=False));print('selected',selected)
