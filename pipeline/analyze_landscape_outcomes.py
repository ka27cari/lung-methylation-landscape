#!/usr/bin/env python3
import os,sys,json
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.preprocessing import RobustScaler
from sklearn.metrics import roc_auc_score,r2_score
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LinearRegression
from scipy.spatial.distance import cdist,pdist
from scipy.stats import spearmanr,linregress,pearsonr,rankdata
ROOT=Path(__file__).resolve().parents[1] / "rebuild";sys.path.insert(0,str(ROOT.parent/'src'))
from pfl.io_source import load_source,read_m_matrix
from pfl.external_projection import frozen_projection
from pfl.representations import derive_orientation,fit_counterfactual_model,apply_counterfactual_model
MOD=ROOT/'cache/model';LONG=ROOT/'cache/longitudinal';EXT=ROOT/'cache/external';TAB=ROOT/'analysis/tables';FIG=ROOT/'analysis/figure_data';DRAW=ROOT/'cache/posterior_draws'
for p in [MOD,LONG,EXT,TAB,FIG,DRAW]:p.mkdir(parents=True,exist_ok=True)
CELL12=['Bas','Bmem','Bnv','CD4mem','CD4nv','CD8mem','CD8nv','Eos','Mono','Neu','NK','Treg']; ME=[f'ME{i}' for i in range(7)]
X,ids,cpg,meta=load_source(ROOT/'data');ori=derive_orientation(X,meta);score=pd.read_csv(MOD/'SOURCE_PROGRAM_SCORES.csv');Pcols=[c for c in score if c.startswith('program_')];Z=score[Pcols].to_numpy(float);rt=np.load(MOD/'MODEL_RUNTIME_OUTPUT.npz',allow_pickle=True);cent=rt['centroids'];s=rt['state_scales'];tau=float(rt['tau'][0]);terms=rt['terminal_states'].astype(int).tolist();lab=rt['hard_state'].astype(int);P=rt['posterior'];TM=rt['TM'];K=len(cent);ctrl=meta.group.eq('control').to_numpy();pred=meta.group.eq('prediagnostic').to_numpy();tt=pd.to_numeric(meta.time,errors='coerce').to_numpy(float)
load=np.load(MOD/'SELECTED_PROGRAM_LOADINGS.npz',allow_pickle=True);H=load['loadings'];pcen=load['center'];pscale=load['scale']
def post(Q):
 d=cdist(Q,cent,'sqeuclidean');a=-d/(2*tau*s[None,:]);a-=a.max(1,keepdims=True);e=np.exp(a);return e/e.sum(1,keepdims=True)
def robust_affine(query,A,B):
 medA=np.median(A,0);medB=np.median(B,0);qA=np.quantile(A,[.25,.75],axis=0);qB=np.quantile(B,[.25,.75],axis=0);den=qA[1]-qA[0];sc=np.where(den>1e-8,(qB[1]-qB[0])/den,1.0);sc=np.clip(sc,.05,20);return (query-medA)*sc+medB

def _pearson_fast(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float);x=x-x.mean();y=y-y.mean();den=np.sqrt(np.dot(x,x)*np.dot(y,y));return float(np.dot(x,y)/den) if den>0 else np.nan
def fast_spearman(x,y):
 return _pearson_fast(rankdata(np.asarray(x,float),method='average'),rankdata(np.asarray(y,float),method='average'))
def fast_auc_binary(y,scores):
 y=np.asarray(y,int);scores=np.asarray(scores,float);n1=int(y.sum());n0=len(y)-n1
 if n0==0 or n1==0:return np.nan
 r=rankdata(scores,method='average');return float((r[y==1].sum()-n1*(n1+1)/2)/(n0*n1))
cf=fit_counterfactual_model(X,meta,ori,.02,10,'stability_null');fill=np.nanmedian(X[ctrl],axis=0);Ssrc=apply_counterfactual_model(X,cf,rank=True)@H.T
# Longitudinal projection. m_te contains the source samples as an internal platform bridge;
# only source samples determine the robust affine map from m_te program coordinates to source coordinates.
lm=pd.read_csv(ROOT/'data/longitudinal.csv');Xte,idte,cpgte=read_m_matrix(ROOT/'data/m_te.csv',None);pos={str(x):i for i,x in enumerate(cpgte)};keep=np.array([pos.get(str(x),-1) for x in cpg],int)
if np.min(keep)<0:raise RuntimeError('m_te missing source CpGs')
Xte=Xte[:,keep];Xte=np.where(np.isfinite(Xte),Xte,fill);Ste=apply_counterfactual_model(Xte,cf,rank=True)@H.T;tp={str(s0):i for i,s0 in enumerate(idte)};srcidx=np.array([tp[str(s0)] for s0 in ids]);lidx=np.array([tp[str(s0)] for s0 in lm.sentrix_id.astype(str)]);Sl=robust_affine(Ste,Ste[srcidx],Ssrc)[lidx];Zl=(Sl-pcen)/pscale;Lm=lm.set_index(lm.sentrix_id.astype(str),drop=False).loc[np.asarray(idte)[lidx]].reset_index(drop=True);lids=Lm.sentrix_id.astype(str).to_numpy();Pl=post(Zl);tml=Pl[:,terms].sum(1)
Lscore=pd.DataFrame({'sentrix_id':lids});
for j,c in enumerate(Pcols):Lscore[c]=Zl[:,j]
Lscore.to_csv(LONG/'LONGITUDINAL_PROGRAM_SCORES.csv',index=False);Lpost=Lm.copy();Lpost['hard_state']=Pl.argmax(1);Lpost['TM']=tml
for j,c in enumerate(Pcols):Lpost[c]=Zl[:,j]
for k in range(K):Lpost[f'P_state_{k:02d}']=Pl[:,k]
Lpost.to_csv(LONG/'LONGITUDINAL_POSTERIOR.csv',index=False);Lpost.to_csv(LONG/'LONGITUDINAL_MICROSTATE_POSTERIOR_WITH_CONTEXT.csv',index=False)
# External lower-platform projection: source-only masked-to-full linear adapter.
em=pd.read_csv(ROOT/'data/ph-ext.csv');Xe,eids,ecpg=read_m_matrix(ROOT/'data/m-ext.csv',em['sample'].astype(str));Em=em.set_index(em['sample'].astype(str),drop=False).loc[eids].reset_index(drop=True);pe={str(x):i for i,x in enumerate(ecpg)};ks=np.array([i for i,x in enumerate(cpg) if str(x) in pe],int);ke=np.array([pe[str(cpg[i])] for i in ks],int);Xm=X[:,ks];Xem=Xe[:,ke];om={'direction':ori['direction'][ks],'median_delta':ori['median_delta'][ks]};cfm=fit_counterfactual_model(Xm,meta,om,.02,10,'stability_null');Hmask=H[:,ks]/np.maximum(np.linalg.norm(H[:,ks],axis=1,keepdims=True),1e-12);Ssm=apply_counterfactual_model(Xm,cfm,rank=True)@Hmask.T;Sem=apply_counterfactual_model(Xem,cfm,rank=True)@Hmask.T
# adapter audit keeps matched source pairs in the same fold
groups=meta.matched_index.astype(str).to_numpy();predcv=np.zeros_like(Ssrc);gkf=GroupKFold(5)
for tr,va in gkf.split(Ssm,groups=groups):predcv[va]=LinearRegression().fit(Ssm[tr],Ssrc[tr]).predict(Ssm[va])
ar=[]
for j in range(Ssrc.shape[1]):ar.append(dict(program=j+1,cv_R2=r2_score(Ssrc[:,j],predcv[:,j]),cv_pearson_r=float(pearsonr(Ssrc[:,j],predcv[:,j]).statistic),cv_MAE=float(np.mean(np.abs(Ssrc[:,j]-predcv[:,j]))),full_score_SD=float(np.std(Ssrc[:,j],ddof=0))))
pd.DataFrame(ar).to_csv(TAB/'EXTERNAL_ADAPTER_SOURCE_CV.csv',index=False);pd.DataFrame([dict(n_source_samples=len(Ssrc),n_source_pairs=meta.matched_index.nunique(),overlap_probes=len(ks),full_model_probes=len(cpg),variance_weighted_cv_R2=r2_score(Ssrc,predcv,multioutput='variance_weighted'),median_program_cv_R2=pd.DataFrame(ar).cv_R2.median(),min_program_cv_R2=pd.DataFrame(ar).cv_R2.min(),median_program_cv_pearson_r=pd.DataFrame(ar).cv_pearson_r.median())]).to_csv(TAB/'EXTERNAL_ADAPTER_SOURCE_CV_SUMMARY.csv',index=False)
adapter=LinearRegression().fit(Ssm,Ssrc);Sext=adapter.predict(Sem);Ze_raw=(Sext-pcen)/pscale
projection=frozen_projection(Ze_raw, Z, ctrl, Em.group.eq('Control').to_numpy(), cent, s, tau, terms)
Ze=projection['scores'];Pe=projection['posterior'];tme=projection['terminal_mass']
np.savez_compressed(DRAW/'EXTERNAL_SOURCE_ADAPTER.npz',coef=adapter.coef_,intercept=adapter.intercept_,overlap_source_indices=ks,overlap_external_indices=ke)
# Export both the uncalibrated reference and the primary calibrated projection.
Escore=Em[['sample','group']].copy()
Escore_raw=Escore.copy()
for j,c in enumerate(Pcols):
    Escore[c]=Ze[:,j]
    Escore_raw[c]=projection['raw_scores'][:,j]
Escore.to_csv(EXT/'EXTERNAL_PROGRAM_SCORES.csv',index=False)
Escore_raw.to_csv(EXT/'EXTERNAL_PROGRAM_SCORES_RAW.csv',index=False)
Epost=Em[['sample','group']].copy()
Epost['hard_state']=Pe.argmax(1);Epost['MAP_state']=Epost['hard_state'];Epost['TM']=tme
for k in range(K): Epost[f'P_state_{k:02d}']=Pe[:,k]
Epost_raw=Em[['sample','group']].copy()
Epost_raw['TM']=projection['raw_terminal_mass']
for k in range(K): Epost_raw[f'P_state_{k:02d}']=projection['raw_posterior'][:,k]
Epost_raw.to_csv(EXT/'EXTERNAL_POSTERIOR_RAW.csv',index=False)
# source-anchored expected context, explicitly not a direct external measurement
ctxcols=['Mye','Lym','Neu','Mono','CD8tot','Btot','NK','Treg']+ME;spctx=[]
for st in range(K):
 g=meta[lab==st];row={'state':st,'n':len(g)}
 for c in ctxcols:row[c]=float(pd.to_numeric(g[c],errors='coerce').mean()) if c in g else np.nan
 spctx.append(row)
spctx=pd.DataFrame(spctx).set_index('state');spctx.reset_index().to_csv(EXT/'SOURCE_STATE_CONTEXT_MEANS.csv',index=False)
for c in ctxcols:Epost[f'source_expected_{c}']=Pe@spctx[c].to_numpy(float)
Epost.to_csv(EXT/'EXTERNAL_POSTERIOR.csv',index=False);Epost.to_csv(EXT/'EXTERNAL_STATE_POSTERIOR.csv',index=False)
# paired longitudinal table
pairs=[]
for mid,g in Lpost.assign(_i=np.arange(len(Lpost))).groupby('matched_index'):
 p=g[g.group.eq('prediagnostic')];c=g[g.group.eq('case')]
 if len(p)==1 and len(c)==1:pairs.append(dict(matched_index=mid,time=float(p.time.iloc[0]),pred_TM=float(p.TM.iloc[0]),case_TM=float(c.TM.iloc[0]),pred_row=int(p._i.iloc[0]),case_row=int(c._i.iloc[0])))
pairs=pd.DataFrame(pairs);pairs.to_csv(FIG/'LONGITUDINAL_TERMINAL_MASS_PAIRS.csv',index=False);scat=pd.concat([pairs[['matched_index','time','pred_TM']].rename(columns={'pred_TM':'TM'}).assign(group='prediagnostic'),pairs[['matched_index','case_TM']].rename(columns={'case_TM':'TM'}).assign(time=0.0,group='case')],ignore_index=True);scat.to_csv(FIG/'LONGITUDINAL_SINK_WEIGHTED_TM_SCATTER_DATA.csv',index=False);pairs.assign(delta_TM=pairs.case_TM-pairs.pred_TM).to_csv(TAB/'LONGITUDINAL_SINK_WEIGHTED_TM.csv',index=False)
# external bootstrap summary
rng=np.random.default_rng(20260928);y=Epost.group.eq('LUAD').astype(int).to_numpy();ec=Epost[y==0].reset_index(drop=True);ea=Epost[y==1].reset_index(drop=True);boot=np.empty((5000,2))
ec_tm=ec.TM.to_numpy(float);ea_tm=ea.TM.to_numpy(float);n0=len(ec_tm);n1=len(ea_tm);yy_boot=np.r_[np.zeros(n0,dtype=int),np.ones(n1,dtype=int)]
for b in range(len(boot)):
 ctm=ec_tm[rng.integers(0,n0,n0)];atm=ea_tm[rng.integers(0,n1,n1)];ss=np.r_[ctm,atm];boot[b]=[np.median(atm)-np.median(ctm),fast_auc_binary(yy_boot,ss)]
ext_summary=dict(n_control=len(ec),n_LUAD=len(ea),probe_overlap=len(ks),source_model_probes=len(cpg),projection='source_control_median_IQR',raw_AUC=float(fast_auc_binary(y,projection['raw_terminal_mass'])),control_median=float(np.median(ec_tm)),LUAD_median=float(np.median(ea_tm)),median_gap=float(np.median(ea_tm)-np.median(ec_tm)),gap_ci025=float(np.quantile(boot[:,0],.025)),gap_ci975=float(np.quantile(boot[:,0],.975)),AUC=float(fast_auc_binary(y,tme)),AUC_ci025=float(np.quantile(boot[:,1],.025)),AUC_ci975=float(np.quantile(boot[:,1],.975)));pd.DataFrame([ext_summary]).to_csv(TAB/'EXTERNAL_VALIDATION_SUMMARY.csv',index=False);pd.DataFrame({'metric':['TM_median_gap','TM_AUC'],'estimate':[ext_summary['median_gap'],ext_summary['AUC']],'ci025':[ext_summary['gap_ci025'],ext_summary['AUC_ci025']],'ci975':[ext_summary['gap_ci975'],ext_summary['AUC_ci975']]}).to_csv(TAB/'EXTERNAL_BOOTSTRAP_SUMMARY.csv',index=False);np.savez_compressed(DRAW/'EXTERNAL_BOOTSTRAP_DRAWS.npz',median_gap=boot[:,0],auc=boot[:,1])
# macrostate roles and coarse transport: final interpretive flow uses Sinkhorn net95; LP95 is discovery-only
roles=pd.read_csv(TAB/'FINAL_STATE_ROLE_CLASSIFICATION.csv');role=roles.set_index('state').macrostate_role.to_dict();G=pd.read_csv(MOD/'SINKHORN_NET_CURRENT_95.csv',index_col=0).to_numpy(float);macro=['Source','Intermediate','Terminal'];coarse=[];M=np.zeros((3,3))
for ai,a in enumerate(macro):
 for bi,b in enumerate(macro):
  ia=[i for i in range(K) if role[i]==a];ib=[i for i in range(K) if role[i]==b];mass=float(G[np.ix_(ia,ib)].sum()) if ia and ib else 0.0;M[ai,bi]=mass;coarse.append(dict(origin_role=a,destination_role=b,mass=mass))
pd.DataFrame(coarse).to_csv(TAB/'LANDSCAPE_COARSE_FLOW.csv',index=False);pd.DataFrame(M,index=macro,columns=macro).to_csv(TAB/'MACROSTATE_OT_PLAN_FINAL.csv')
# terminal route decomposition from Sinkhorn net-current completion
src=[i for i in range(K) if role[i]=='Source'];inter=[i for i in range(K) if role[i]=='Intermediate'];routes=[]
for t in terms:
 ds=float(G[src,t].sum());di=float(G[inter,t].sum());den=ds+di;frac=ds/den if den>0 else np.nan;rtype='Direct' if frac>=2/3 else ('Intermediate-mediated' if frac<=1/3 else 'Mixed');routes.append(dict(terminal_state=t,source_direct_inflow=ds,intermediate_inflow=di,direct_fraction=frac,relay_fraction=(1-frac if np.isfinite(frac) else np.nan),route_type=rtype,transport_basis='Sinkhorn_net95'))
routes=pd.DataFrame(routes);routes.to_csv(TAB/'TERMINAL_ROUTE_SUMMARY.csv',index=False);routes.to_csv(TAB/'TERMINAL_OT_ROUTE_SUMMARY.csv',index=False);routes.to_csv(TAB/'MACROSTATE_OT_ROUTES_FINAL.csv',index=False)
# blood CLR helper
def clr(df,med=None,pc=None):
 A=df[CELL12].apply(pd.to_numeric,errors='coerce').to_numpy(float);md=np.nanmedian(A,axis=0) if med is None else med;bad=~np.isfinite(A);A[bad]=np.take(md,np.where(bad)[1]);
 if pc is None:
  pos=A[np.isfinite(A)&(A>0)];pc=float(pos.min()/2)
 A=np.where(A<=0,pc,A);R=np.log(A);R-=R.mean(1,keepdims=True);return R,md,pc
B,bdmed,bdpc=clr(meta);bs=RobustScaler(quantile_range=(10,90)).fit(B[pred]);Bs=bs.transform(B);Bl,_,_=clr(Lm,bdmed,bdpc);Bls=bs.transform(Bl)
# state program-geometry versus blood-composition geometry, excluding only extremely sparse states n_pred<5
states=[k for k in range(K) if (pred&(lab==k)).sum()>=5];PC=np.vstack([np.median(Z[pred&(lab==k)],axis=0) for k in states]);BC=np.vstack([np.median(Bs[pred&(lab==k)],axis=0) for k in states]);dp=pdist(PC);db=pdist(BC);geom_rho=float(fast_spearman(dp,db));perm=np.empty(5000);pred_idx=np.where(pred)[0];pred_states=lab[pred];dp_rank=rankdata(dp,method='average')
for b in range(len(perm)):
 shuffled=rng.permutation(pred_idx);Bc=[np.median(Bs[shuffled[pred_states==k]],axis=0) for k in states];perm[b]=_pearson_fast(dp_rank,rankdata(pdist(np.vstack(Bc)),method='average'))
geom_p=float((1+(perm>=geom_rho).sum())/(1+len(perm)));np.savez_compressed(DRAW/'BLOOD_GEOMETRY_PERMUTATION_DRAWS.npz',rho_permutation=perm,observed=np.array([geom_rho]));pd.DataFrame([dict(rho=geom_rho,permutation_p=geom_p,n_states=len(states))]).to_csv(TAB/'PROGRAM_BLOOD_GEOMETRY_ASSOCIATION.csv',index=False);pd.DataFrame([dict(rho=geom_rho,permutation_p=geom_p,n_states=len(states))]).to_csv(TAB/'STATE_GEOMETRY_BLOOD_CONCORDANCE.csv',index=False)
# terminal-aware longitudinal methylation/blood displacement
progsc=RobustScaler(quantile_range=(10,90)).fit(Z);Zs=progsc.transform(Z);Zls=progsc.transform(Zl);sources=[k for k in range(K) if role[k]=='Source'];sm=pred&np.isin(lab,sources);fn=np.median(Zs[sm],axis=0);hn=np.median(Bs[sm],axis=0);fd=[];hd=[]
for t in terms:
 m=pred&(lab==t);fv=np.median(Zs[m],axis=0)-fn;hv=np.median(Bs[m],axis=0)-hn;fd.append(fv/(np.linalg.norm(fv)+1e-12));hd.append(hv/(np.linalg.norm(hv)+1e-12))
fd=np.asarray(fd);hd=np.asarray(hd);poolf=np.median(Zs[pred&np.isin(lab,terms)],axis=0)-fn;poolh=np.median(Bs[pred&np.isin(lab,terms)],axis=0)-hn;poolf/=np.linalg.norm(poolf)+1e-12;poolh/=np.linalg.norm(poolh)+1e-12
aware_x=[];aware_y=[];pool_x=[];pool_y=[]
for _,r in pairs.iterrows():
 ip=int(r.pred_row);ic=int(r.case_row);w=Pl[ip,terms];w/=max(w.sum(),1e-12);af=(w[:,None]*fd).sum(0);ah=(w[:,None]*hd).sum(0);af/=np.linalg.norm(af)+1e-12;ah/=np.linalg.norm(ah)+1e-12;dz=Zls[ic]-Zls[ip];dh=Bls[ic]-Bls[ip];aware_x.append(float(dz@af));aware_y.append(float(dh@ah));pool_x.append(float(dz@poolf));pool_y.append(float(dh@poolh))
aware_x=np.asarray(aware_x,float);aware_y=np.asarray(aware_y,float);pool_x=np.asarray(pool_x,float);pool_y=np.asarray(pool_y,float);aware_r=float(fast_spearman(aware_x,aware_y));pool_r=float(fast_spearman(pool_x,pool_y));cdf=pairs[['matched_index','time']].copy();cdf['terminal_aware_fate_displacement']=aware_x;cdf['terminal_aware_hemo_displacement']=aware_y;cdf['pooled_fate_displacement']=pool_x;cdf['pooled_hemo_displacement']=pool_y;cdf.to_csv(FIG/'LONGITUDINAL_TERMINAL_AWARE_COUPLING.csv',index=False);cdf.to_csv(TAB/'LONGITUDINAL_FATE_HEMO_DISPLACEMENT.csv',index=False)
n=len(cdf);boot_a=[];boot_p=[];perm_a=[];aware_x_rank=rankdata(aware_x,method='average');aware_y_rank=rankdata(aware_y,method='average')
for b in range(5000):
 ii=rng.integers(0,n,n);boot_a.append(fast_spearman(aware_x[ii],aware_y[ii]));boot_p.append(fast_spearman(pool_x[ii],pool_y[ii]));jj=rng.permutation(n);perm_a.append(_pearson_fast(aware_x_rank,aware_y_rank[jj]))
boot_a=np.asarray(boot_a);boot_p=np.asarray(boot_p);perm_a=np.asarray(perm_a);np.savez_compressed(DRAW/'LONGITUDINAL_COUPLING_DRAWS.npz',aware_bootstrap=boot_a,pooled_bootstrap=boot_p,aware_permutation=perm_a,aware_observed=np.array([aware_r]),pooled_observed=np.array([pool_r]));coupling_p=float((1+(perm_a>=aware_r).sum())/(1+len(perm_a)));csummary=pd.DataFrame([dict(pooled_rho=pool_r,terminal_aware_rho=aware_r,aware_boot_ci025=float(np.nanquantile(boot_a,.025)),aware_boot_ci975=float(np.nanquantile(boot_a,.975)),permutation_p=coupling_p,n_pairs=n)]);csummary.to_csv(TAB/'LONGITUDINAL_COUPLING_SUMMARY.csv',index=False);csummary.to_csv(TAB/'LONGITUDINAL_FATE_HEMO_SUMMARY.csv',index=False)
# source state fingerprints and compatibility STATE_TABLE
ctx=[]
for k in range(K):
 m=pred&(lab==k);row=dict(state=k,role=role[k],n=int((lab==k).sum()),n_pred=int(m.sum()),control_n=int((ctrl&(lab==k)).sum()),control_fraction=float((ctrl&(lab==k)).sum()/max(1,(lab==k).sum())),control_occupancy=float(np.mean(lab[ctrl]==k)),far_occupancy=float(np.mean(lab[pred&(tt>=5)]==k)),near_occupancy=float(np.mean(lab[pred&(tt<3)]==k)),median_residual_time=float(np.nanmedian(tt[m])) if m.sum() else np.nan,LUAD_fraction=float((meta.loc[m,'type'].astype(str)=='LUAD').mean()) if m.sum() else np.nan,LUSC_fraction=float((meta.loc[m,'type'].astype(str)=='LUSC').mean()) if m.sum() else np.nan,SCLC_fraction=float((meta.loc[m,'type'].astype(str)=='SCLC').mean()) if m.sum() else np.nan,is_terminal=bool(k in terms))
 for c in ['Mye','Lym','Neu','Mono','CD8tot','Btot','Treg','NK']+ME:row[c]=float(pd.to_numeric(meta.loc[m,c],errors='coerce').median()) if m.sum() else np.nan
 ctx.append(row)
ctx=pd.DataFrame(ctx);ctx.to_csv(TAB/'SOURCE_STATE_CONTEXT_SUMMARY.csv',index=False);ctx.to_csv(MOD/'STATE_TABLE.csv',index=False);ctx.to_csv(TAB/'STATE_TABLE.csv',index=False);ctx[ctx.role.eq('Terminal')].to_csv(TAB/'TERMINAL_CONTEXT_FINGERPRINTS.csv',index=False);ctx.groupby('role',as_index=False).median(numeric_only=True).to_csv(TAB/'SOURCE_MACROSTATE_CONTEXT_PROFILE.csv',index=False)
# macrostate posterior exports
role_states={r:[k for k in range(K) if role[k]==r] for r in macro}
sp=pd.DataFrame({'sentrix_id':ids,'group':meta.group.astype(str),'time':tt})
lp=Lpost[['sentrix_id','matched_index','group','time']].copy()
for r,ss in role_states.items():sp[r+'_macrostate']=P[:,ss].sum(1);lp[r+'_macrostate']=Pl[:,ss].sum(1)
sp.to_csv(TAB/'SOURCE_MACROSTATE_POSTERIOR_FINAL.csv',index=False);lp.to_csv(TAB/'LONGITUDINAL_MACROSTATE_POSTERIOR_FINAL.csv',index=False)
# Route-aware scores.  Terminal route decomposition is data-driven; scalar alternatives are sensitivity outputs, not model-selection criteria.
stage=np.array([0.0 if role[k]=='Source' else (0.5 if role[k]=='Intermediate' else 1.0) for k in range(K)]);stage_source=P@stage;stage_long=Pl@stage
Q=G.copy();np.fill_diagonal(Q,0);rs=Q.sum(1,keepdims=True);Tmat=np.divide(Q,rs,out=np.zeros_like(Q),where=rs>0)
for t in terms:Tmat[t,:]=0;Tmat[t,t]=1
trans=[k for k in range(K) if k not in terms];Qtt=Tmat[np.ix_(trans,trans)];Rtt=Tmat[np.ix_(trans,terms)];I=np.eye(len(trans))
try:
 Nmat=np.linalg.inv(I-Qtt+1e-8*I);absorb=Nmat@Rtt;steps=Nmat@np.ones(len(trans))
except Exception:
 absorb=np.zeros((len(trans),len(terms)));steps=np.full(len(trans),np.nan)
Abs=np.zeros((K,len(terms)));Abs[trans]=absorb
for j,t in enumerate(terms):Abs[t,j]=1
potential=np.ones(K);potential[trans]=1/(1+np.maximum(steps,0));route_source=P@potential;route_long=Pl@potential;event_src=P@Abs;event_long=Pl@Abs
direct_frac=np.array([float(routes.set_index('terminal_state').loc[t,'direct_fraction']) for t in terms]);direct_tm_src=event_src@direct_frac;relay_tm_src=event_src@(1-direct_frac);direct_tm_long=event_long@direct_frac;relay_tm_long=event_long@(1-direct_frac)
route_src_df=pd.DataFrame({'sentrix_id':ids,'TM':TM,'role_stage':stage_source,'route_potential':route_source,'eventual_direct_affinity':direct_tm_src,'eventual_relay_affinity':relay_tm_src});route_long_df=pd.DataFrame({'sentrix_id':lids,'TM':tml,'role_stage':stage_long,'route_potential':route_long,'eventual_direct_affinity':direct_tm_long,'eventual_relay_affinity':relay_tm_long})
for j,t in enumerate(terms):route_src_df[f'route_to_T{t}']=event_src[:,j];route_long_df[f'route_to_T{t}']=event_long[:,j]
route_src_df.to_csv(MOD/'SOURCE_ROUTE_SCORES.csv',index=False);route_long_df.to_csv(LONG/'LONGITUDINAL_ROUTE_SCORES.csv',index=False);pd.DataFrame(Abs,columns=[f'T{t}' for t in terms]).assign(state=np.arange(K)).to_csv(TAB/'STATE_TO_TERMINAL_ABSORPTION.csv',index=False)
route_sum=[]
for dataset,mask,timev,scoreset in [('source',pred,-tt[pred],{'TM':TM[pred],'role_stage':stage_source[pred],'route_potential':route_source[pred],'eventual_direct_affinity':direct_tm_src[pred],'eventual_relay_affinity':relay_tm_src[pred]}),('longitudinal_prediagnostic',Lm.group.eq('prediagnostic').to_numpy(),-pd.to_numeric(Lm.loc[Lm.group.eq('prediagnostic'),'time'],errors='coerce').to_numpy(float),None)]:
 if scoreset is None:
  m=Lm.group.eq('prediagnostic').to_numpy();scoreset={'TM':tml[m],'role_stage':stage_long[m],'route_potential':route_long[m],'eventual_direct_affinity':direct_tm_long[m],'eventual_relay_affinity':relay_tm_long[m]}
 for name,yv in scoreset.items():
  lr=linregress(timev,yv);route_sum.append(dict(dataset=dataset,score=name,spearman_rho=float(spearmanr(yv,timev).statistic),linear_r2=float(lr.rvalue**2),slope=float(lr.slope)))
pd.DataFrame(route_sum).to_csv(TAB/'ROUTE_AWARE_TEMPORAL_COMPARISON.csv',index=False)
# Route-stratified terminal occupancy.  These are descriptive outputs from the frozen landscape.
direct_terms=[t for t,f in zip(terms,direct_frac) if f>=2/3]
relay_terms=[t for t,f in zip(terms,direct_frac) if f<=1/3]
route_group_rows=[]
def _route_group_summary(label,src_score,long_score):
    row=dict(route_group=label,
             source_control_median=float(np.median(src_score[ctrl])),
             source_far_median=float(np.median(src_score[pred&(tt>=5)])),
             source_near_median=float(np.median(src_score[pred&(tt<3)])),
             source_proximity_rho=float(spearmanr(src_score[pred],-tt[pred]).statistic))
    tmp=Lm[['matched_index','group','time']].copy();tmp['_s']=long_score
    pr=[]
    for mid,g in tmp.groupby('matched_index'):
        a=g[g.group.eq('prediagnostic')];b=g[g.group.eq('case')]
        if len(a)==1 and len(b)==1:pr.append((float(a.time.iloc[0]),float(a._s.iloc[0]),float(b._s.iloc[0])))
    q=np.asarray(pr,float)
    row.update(long_far_median=float(np.median(q[q[:,0]>=5,1])),long_near_median=float(np.median(q[q[:,0]<3,1])),long_case_median=float(np.median(q[:,2])),long_caseward_fraction=float(np.mean(q[:,2]>q[:,1])))
    return row
route_group_rows.append(_route_group_summary('All terminals',TM,tml))
if direct_terms: route_group_rows.append(_route_group_summary('Direct terminals',P[:,direct_terms].sum(1),Pl[:,direct_terms].sum(1)))
if relay_terms: route_group_rows.append(_route_group_summary('Intermediate-mediated terminals',P[:,relay_terms].sum(1),Pl[:,relay_terms].sum(1)))
pd.DataFrame(route_group_rows).to_csv(TAB/'ROUTE_STRATIFIED_TERMINAL_MASS.csv',index=False)
# Terminal-specific source and longitudinal horizon summaries.
tr_rows=[]
for t in terms:
    ss=P[:,t];ls=Pl[:,t];tmp=Lm[['matched_index','group','time']].copy();tmp['_s']=ls;pr=[]
    for mid,g in tmp.groupby('matched_index'):
        a=g[g.group.eq('prediagnostic')];b=g[g.group.eq('case')]
        if len(a)==1 and len(b)==1:pr.append((float(a.time.iloc[0]),float(a._s.iloc[0]),float(b._s.iloc[0])))
    q=np.asarray(pr,float)
    tr_rows.append(dict(terminal_state=int(t),route_type=str(routes.set_index('terminal_state').loc[t,'route_type']),
        source_far_median=float(np.median(ss[pred&(tt>=5)])),source_near_median=float(np.median(ss[pred&(tt<3)])),source_proximity_rho=float(spearmanr(ss[pred],-tt[pred]).statistic),
        long_far_median=float(np.median(q[q[:,0]>=5,1])),long_near_median=float(np.median(q[q[:,0]<3,1])),long_case_median=float(np.median(q[:,2])),long_caseward_fraction=float(np.mean(q[:,2]>q[:,1]))))
pd.DataFrame(tr_rows).to_csv(TAB/'TERMINAL_ROUTE_TEMPORAL_PROFILES.csv',index=False)
# Convex route-aware score grid: secondary sensitivity check for whether route information improves a scalar temporal summary.
grid=[]
long_pred=Lm.group.eq('prediagnostic').to_numpy(); long_time=-pd.to_numeric(Lm.loc[long_pred,'time'],errors='coerce').to_numpy(float)
for route_name,src_alt,long_alt in [('role_stage',stage_source,stage_long),('route_potential',route_source,route_long)]:
    for lam in np.linspace(0,1,21):
        src_y=(1-lam)*TM+lam*src_alt; long_y=(1-lam)*tml+lam*long_alt
        # longitudinal pair medians/caseward
        tmp=Lm[['matched_index','group','time']].copy();tmp['_s']=long_y;pr=[]
        for mid,g in tmp.groupby('matched_index'):
            a=g[g.group.eq('prediagnostic')];b=g[g.group.eq('case')]
            if len(a)==1 and len(b)==1:pr.append((float(a.time.iloc[0]),float(a._s.iloc[0]),float(b._s.iloc[0])))
        q=np.asarray(pr,float)
        grid.append(dict(route_score=route_name,lambda_route=float(lam),source_spearman=float(spearmanr(src_y[pred],-tt[pred]).statistic),longitudinal_pred_spearman=float(spearmanr(long_y[long_pred],long_time).statistic),long_far_median=float(np.median(q[q[:,0]>=5,1])),long_near_median=float(np.median(q[q[:,0]<3,1])),long_case_median=float(np.median(q[:,2])),long_caseward_fraction=float(np.mean(q[:,2]>q[:,1]))))
pd.DataFrame(grid).to_csv(TAB/'ROUTE_WEIGHTED_SCORE_GRID.csv',index=False)
# terminal-specific destination affinity vs source residual time
hrows=[]
cond=P[:,terms]/np.maximum(P[:,terms].sum(1,keepdims=True),1e-12)
for j,t in enumerate(terms):
 hrows.append(dict(terminal_state=t,source_affinity_vs_proximity_rho=float(spearmanr(cond[pred,j],-tt[pred]).statistic),route_type=str(routes.set_index('terminal_state').loc[t,'route_type']),direct_fraction=float(direct_frac[j])))
pd.DataFrame(hrows).to_csv(TAB/'TERMINAL_ROUTE_HORIZON_SLOPES.csv',index=False)
# final summaries
summary=dict(source_control_median=float(np.median(TM[ctrl])),source_far_median=float(np.median(TM[pred&(tt>=5)])),source_near_median=float(np.median(TM[pred&(tt<3)])),source_proximity_rho=float(spearmanr(TM[pred],-tt[pred]).statistic),long_far_median=float(pairs.loc[pairs.time>=5,'pred_TM'].median()),long_near_median=float(pairs.loc[pairs.time<3,'pred_TM'].median()),long_case_median=float(pairs.case_TM.median()),long_caseward_fraction=float((pairs.case_TM>pairs.pred_TM).mean()),terminal_aware_coupling_rho=aware_r,pooled_coupling_rho=pool_r,blood_geometry_rho=geom_rho,blood_geometry_permutation_p=geom_p,external_auc=ext_summary['AUC'])
pd.DataFrame([summary]).to_csv(TAB/'ANALYSIS_SUMMARY.csv',index=False);json.dump(summary,open(TAB/'ANALYSIS_SUMMARY.json','w'),indent=2)
print(summary);print('routes');print(routes.to_string(index=False));print('route temporal');print(pd.DataFrame(route_sum).to_string(index=False))
