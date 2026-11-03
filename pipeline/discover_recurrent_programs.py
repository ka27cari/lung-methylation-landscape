#!/usr/bin/env python3
import os,sys,json,yaml,argparse
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd
from joblib import Parallel,delayed
ROOT=Path(__file__).resolve().parents[1] / "rebuild";sys.path.insert(0,str(ROOT.parent/'src'))
from pfl.io_source import load_source
from pfl.representations import derive_orientation,counterfactual_representation
from pfl.contexts import kmeans_contexts
from pfl.discovery import discover_context,recurrent_families
CFG=yaml.safe_load((ROOT.parent/'config/model.yaml').read_text())

def main():
 out=ROOT/'cache/discovery';out.mkdir(parents=True,exist_ok=True);X,ids,cpg,meta=load_source(ROOT/'data');ori=derive_orientation(X,meta);rep=CFG['representation'];_,A,anchors,_=counterfactual_representation(X,meta,ori,rep['anchor_fraction'],rep['control_neighbors'],rep['anchor_mode'])
 pd.DataFrame({'cpg':cpg,'direction':ori['direction'],'reliability':ori['reliability'],'median_pair_delta':ori['median_delta']}).to_csv(out/'SOURCE_CPG_ORIENTATION.csv',index=False);pd.DataFrame({'anchor_cpg':cpg[anchors]}).to_csv(out/'SOURCE_COUNTERFACTUAL_ANCHORS.csv',index=False)
 contexts=kmeans_contexts(A,meta,CFG['local_contexts']['k_grid'],CFG['local_contexts']['min_context_n'],CFG['random_seed'])
 def work(i,c):
  comp,H=discover_context(c,A,CFG['nmf']['ranks'],CFG['nmf']['seeds'],CFG['nmf']['max_iter'],rep['name'],'kmeans_local',cpg);fam,mem,FH=recurrent_families(comp,H,CFG['recurrence']['within_context_cosines'],CFG['recurrence']['min_distinct_ranks']);return i,comp,fam,mem,FH
 res=Parallel(n_jobs=CFG['nmf']['parallel_jobs'],backend='loky',max_nbytes='8M')(delayed(work)(i,c) for i,c in enumerate(contexts));af=[];ah=[];am=[];ac=[];off=0
 for _,comp,fam,mem,FH in sorted(res):
  if fam.empty:continue
  comp=comp.copy();comp['global_component_row']=np.arange(len(comp))+off;mem=mem.copy();mem['component_row_local']=mem['component_row'];mem['component_row']=mem['component_row']+off;off+=len(comp);ac.append(comp);am.append(mem);af.append(fam);ah.append(FH)
 fam=pd.concat(af,ignore_index=True);fam['family_row']=np.arange(len(fam));Hfam=np.vstack(ah).astype(np.float32);mem=pd.concat(am,ignore_index=True);comp=pd.concat(ac,ignore_index=True);fam.to_csv(out/'RECURRENT_FAMILIES_ALL.csv',index=False);mem.to_csv(out/'FAMILY_MEMBERS.csv',index=False);comp.to_csv(out/'RANK_CONSENSUS_COMPONENTS.csv',index=False);np.savez_compressed(out/'RECURRENT_FAMILY_LOADINGS.npz',loadings=Hfam,cpg=np.asarray(cpg,dtype='U'))
 json.dump({'n_source_samples':len(meta),'n_cpg':len(cpg),'n_anchors':len(anchors),'n_contexts':len(contexts),'n_recurrent_families':len(fam)},open(out/'DISCOVERY_MANIFEST.json','w'),indent=2)
if __name__=='__main__':main()
