from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

def read_m_matrix(path: str | Path, sample_ids=None):
    df=pd.read_csv(path)
    cpg=df.iloc[:,0].astype(str).to_numpy()
    if sample_ids is None:
        cols=list(df.columns[1:])
    else:
        wanted=set(map(str,sample_ids)); cols=[c for c in df.columns[1:] if str(c) in wanted]
    X=df[cols].to_numpy(dtype=np.float32).T
    return X,np.asarray(cols,dtype=str),cpg

def load_source(data_dir: str | Path):
    d=Path(data_dir); meta=pd.read_csv(d/'nested.csv')
    X,ids,cpg=read_m_matrix(d/'m_tr.csv',meta.sentrix_id.astype(str))
    mi=meta.set_index(meta.sentrix_id.astype(str),drop=False);meta2=mi.loc[ids].reset_index(drop=True)
    med=np.nanmedian(X,axis=0);bad=~np.isfinite(X)
    if bad.any():
        X=X.copy();X[bad]=np.take(med,np.where(bad)[1])
    return X,ids,cpg,meta2
