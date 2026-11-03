from __future__ import annotations
import numpy as np
from sklearn.decomposition import NMF
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import AgglomerativeClustering
from .utils import l2_rows

def fit_one(X, rank, seed, max_iter=600, weights=None):
    A=np.asarray(X,dtype=np.float32)
    if weights is not None:
        w=np.sqrt(np.clip(np.asarray(weights,float),0,None))[:,None]
        A=A*w
    model=NMF(n_components=rank,init="nndsvda",solver="cd",beta_loss="frobenius",
              max_iter=max_iter,random_state=int(seed),tol=1e-4)
    W=model.fit_transform(A); H=model.components_
    return W.astype(np.float32), l2_rows(H).astype(np.float32), float(model.reconstruction_err_)

def rank_consensus(X, rank, seeds, max_iter=600, weights=None):
    fits=[fit_one(X,rank,s,max_iter,weights) for s in seeds]
    H0=fits[0][1]; aligned=[H0]; st=[]
    for _,H,_ in fits[1:]:
        S=H0@H.T
        r,c=linear_sum_assignment(-S)
        inv=np.empty(rank,int); inv[r]=c
        Ha=H[inv]
        aligned.append(Ha); st.append(np.diag(H0@Ha.T))
    H=l2_rows(np.mean(np.stack(aligned),axis=0))
    if st: stability=np.mean(np.stack(st),axis=0)
    else: stability=np.ones(rank)
    err=np.mean([x[2] for x in fits])
    W=np.maximum(0,X@H.T)
    return W.astype(np.float32),H.astype(np.float32),stability.astype(float),err

def cluster_components(component_df, loadings, cosine_threshold=.85):
    H=l2_rows(loadings)
    if len(H)==1: labels=np.array([0])
    else:
        dist=np.clip(1-H@H.T,0,2)
        cl=AgglomerativeClustering(n_clusters=None,metric='precomputed',linkage='average',
                                   distance_threshold=1-cosine_threshold)
        labels=cl.fit_predict(dist)
    return labels
