from __future__ import annotations
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def kmeans_contexts(X, meta, k_grid, min_n=50, seed=1):
    """Define local unsupervised discovery contexts in the counterfactual rank space."""
    out=[]
    d=min(12,X.shape[1],max(2,len(X)-2))
    Z=StandardScaler().fit_transform(PCA(d,random_state=seed).fit_transform(X))
    for k in k_grid:
        if k>=len(X): continue
        lab=KMeans(k,n_init=10,random_state=seed+k).fit_predict(Z)
        for j in range(k):
            idx=np.where(lab==j)[0]
            if len(idx)>=min_n:
                out.append({'name':f'kmeans_k{k}_c{j}','indices':idx,'weights':None,'kind':'kmeans_local'})
    return out
