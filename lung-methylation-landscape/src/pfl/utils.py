from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata

def l2_rows(H):
    H = np.asarray(H, float)
    n = np.linalg.norm(H, axis=1, keepdims=True)
    return H / np.maximum(n, 1e-12)

def row_percentile_positive(X):
    X = np.asarray(X, float)
    out = np.zeros_like(X, dtype=np.float32)
    for i in range(X.shape[0]):
        v = X[i]
        pos = v > 0
        if pos.any():
            r = rankdata(v[pos], method="average") / pos.sum()
            out[i, pos] = r
    return out

def file_sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()

def write_json(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path,'w',encoding='utf-8') as f:
        json.dump(obj,f,ensure_ascii=False,indent=2,default=str)
