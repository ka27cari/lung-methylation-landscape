from __future__ import annotations
import numpy as np
from scipy.stats import rankdata
from sklearn.preprocessing import RobustScaler
from sklearn.neighbors import NearestNeighbors
from .utils import row_percentile_positive

def matched_pairs(meta):
    out=[]
    for mid,g in meta.groupby('matched_index',sort=False):
        c=g.index[g.group.eq('control')].tolist();p=g.index[g.group.eq('prediagnostic')].tolist()
        if len(c)==1 and len(p)==1:out.append((mid,c[0],p[0]))
    return out

def derive_orientation(X,meta):
    pairs=matched_pairs(meta);D=np.vstack([X[p]-X[c] for _,c,p in pairs]);med=np.nanmedian(D,axis=0);direction=np.where(med>=0,1.0,-1.0).astype(np.float32);reliability=np.mean(D*direction>0,axis=0).astype(np.float32)
    return dict(pairs=pairs,delta=D.astype(np.float32),direction=direction,reliability=reliability,median_delta=med.astype(np.float32))

def stable_anchor_indices(X,meta,orient,frac,mode='stability_null'):
    ctrl=meta.group.eq('control').to_numpy();mad=np.nanmedian(np.abs(X[ctrl]-np.nanmedian(X[ctrl],axis=0)),axis=0)+1e-6;a=rankdata(mad,method='average')/len(mad)
    if mode=='control_stability':score=a
    else:
        b=rankdata(np.abs(orient['median_delta']),method='average')/len(mad);score=.5*a+.5*b
    n=max(20,int(round(frac*X.shape[1])));return np.argsort(score)[:n]

def counterfactual_representation(X,meta,orient,anchor_frac=.05,k=30,mode='stability_null'):
    ctrl_idx=np.where(meta.group.eq('control').to_numpy())[0];anchors=stable_anchor_indices(X,meta,orient,anchor_frac,mode);scaler=RobustScaler(quantile_range=(10,90)).fit(X[ctrl_idx][:,anchors]);C=scaler.transform(X[ctrl_idx][:,anchors]);Z=scaler.transform(X[:,anchors]);kk=min(k,len(ctrl_idx)-1);nn=NearestNeighbors(n_neighbors=min(kk+1,len(ctrl_idx))).fit(C);loc=nn.kneighbors(Z,return_distance=False);ref=np.zeros_like(X,dtype=np.float32);ctrl_pos={int(v):j for j,v in enumerate(ctrl_idx)}
    for i in range(len(X)):
        js=list(loc[i])
        if i in ctrl_pos:js=[j for j in js if j!=ctrl_pos[i]]
        ref[i]=np.mean(X[ctrl_idx[js[:kk]]],axis=0)
    raw=np.maximum((X-ref)*orient['direction'],0).astype(np.float32);return raw,row_percentile_positive(raw),anchors,ref

def fit_counterfactual_model(X_source,meta,orient,anchor_frac=.05,k=30,mode='stability_null'):
    ctrl_idx=np.where(meta.group.eq('control').to_numpy())[0];anchors=stable_anchor_indices(X_source,meta,orient,anchor_frac,mode);scaler=RobustScaler(quantile_range=(10,90)).fit(X_source[ctrl_idx][:,anchors]);C=scaler.transform(X_source[ctrl_idx][:,anchors]);kk=min(k,len(ctrl_idx));nn=NearestNeighbors(n_neighbors=kk).fit(C)
    return dict(anchors=anchors,scaler=scaler,nn=nn,k=kk,controls=X_source[ctrl_idx].astype(np.float32),direction=orient['direction'])

def apply_counterfactual_model(X_query,model,rank=False):
    Z=model['scaler'].transform(X_query[:,model['anchors']]);loc=model['nn'].kneighbors(Z,return_distance=False)[:,:model['k']];ref=np.stack([model['controls'][js].mean(axis=0) for js in loc]);A=np.maximum((X_query-ref)*model['direction'],0);return row_percentile_positive(A) if rank else A.astype(np.float32)
