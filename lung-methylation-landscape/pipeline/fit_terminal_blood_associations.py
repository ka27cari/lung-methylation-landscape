#!/usr/bin/env python3
"""Terminal-conditioned Bayesian models for risk-associated methylation eigengenes.

Primary covariates: age, sex, smoking status and amount, and source cohort.
Cell-composition PC variants are available as secondary configurations.
"""
import json, yaml
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1] / "rebuild"
import sys
sys.path.insert(0, str(ROOT.parent / 'src'))
from pfl.bayesian import nig_fit as fit_bayesian
TAB=ROOT/'analysis/tables'; DRAW=ROOT/'cache/posterior_draws'; MOD=ROOT/'cache/model'; SPEC=ROOT/'analysis/spec'
for p in [TAB,DRAW,MOD,SPEC]: p.mkdir(parents=True,exist_ok=True)
CFG=yaml.safe_load((ROOT.parent/'config/model.yaml').read_text())
BCFG=CFG['bayesian_me']
meta=pd.read_csv(ROOT/'data/nested.csv')
post=pd.read_csv(MOD/'SOURCE_POSTERIOR.csv')
meta.sentrix_id=meta.sentrix_id.astype(str); post.sentrix_id=post.sentrix_id.astype(str)
meta=meta.set_index('sentrix_id').loc[post.sentrix_id].reset_index()
pred=meta.group.eq('prediagnostic').to_numpy()
terms=list(map(int,json.load(open(TAB/'FINAL_MODEL_SPEC.json'))['terminal_states']))
ME=[f'ME{i}' for i in range(7)]
rng=np.random.default_rng(20261004)

# -----------------------------------------------------------------------------
# Covariates. Continuous variables are imputed/scaled in prediagnostic samples.
# Source cohort is a fixed nuisance adjustment in every declared configuration.
# -----------------------------------------------------------------------------
def numeric(v):
    return pd.to_numeric(v,errors='coerce').to_numpy(float)

def impute_and_z(v,fit_mask=pred):
    raw=numeric(v); miss=~np.isfinite(raw)
    med=float(np.nanmedian(raw[fit_mask])); x=np.where(np.isfinite(raw),raw,med)
    mu=float(x[fit_mask].mean()); sd=float(x[fit_mask].std(ddof=0)); z=(x-mu)/(sd+1e-12)
    return z,dict(median=med,mean=mu,sd=sd,n_imputed=int(miss[fit_mask].sum()))

def binary_z(x,fit_mask=pred):
    x=np.asarray(x,float); mu=float(x[fit_mask].mean()); sd=float(x[fit_mask].std(ddof=0))
    return (x-mu)/(sd+1e-12)

proximity,prox_meta=impute_and_z(-pd.to_numeric(meta.time,errors='coerce'))
age,age_meta=impute_and_z(meta.age)
smk_amnt,smk_meta=impute_and_z(meta.smk_amnt)
sex_M=binary_z(meta.sex.astype(str).str.upper().eq('M').astype(float).to_numpy())
origin_NCC=binary_z(meta.origin.astype(str).str.upper().eq('NCC').astype(float).to_numpy())
smk_status=meta.smk_status.astype(str).str.lower()
smk_current=binary_z(smk_status.eq('current').astype(float).to_numpy())
smk_quit=binary_z(smk_status.eq('quit').astype(float).to_numpy())

# 12-cell composition: zero replacement -> CLR -> standardized PCA.
# PCA is fit only in the prediagnostic population used by the ME regression.
CELL12=['Bas','Bmem','Bnv','CD4mem','CD4nv','CD8mem','CD8nv','Eos','Mono','Neu','NK','Treg']
C=meta[CELL12].apply(pd.to_numeric,errors='coerce').to_numpy(float)
med=np.nanmedian(C[pred],axis=0); bad=~np.isfinite(C); C[bad]=np.take(med,np.where(bad)[1])
pos=C[pred][(C[pred]>0)&np.isfinite(C[pred])]; pc=float(pos.min()/2)
C=np.where(C<=0,pc,C); CLR=np.log(C); CLR-=CLR.mean(1,keepdims=True)
clr_scaler=StandardScaler().fit(CLR[pred]); CLRz=clr_scaler.transform(CLR)
pca=PCA(n_components=3,random_state=20261004).fit(CLRz[pred]); BPC=pca.transform(CLRz)
for j in range(3):
    mu=float(BPC[pred,j].mean()); sd=float(BPC[pred,j].std(ddof=0)); BPC[:,j]=(BPC[:,j]-mu)/(sd+1e-12)

MEZ={}; me_scale={}
for me in ME:
    z,md=impute_and_z(meta[me]); MEZ[me]=z; me_scale[me]=md

# Frozen final-terminal posteriors. Raw P_state_t is used as the soft
# terminal-conditioning weight. Terminal shares are used only in secondary
# terminal-affinity regressions below.
Pterm=np.column_stack([post[f'P_state_{t:02d}'].to_numpy(float) for t in terms])
Pshare=Pterm/np.maximum(Pterm.sum(axis=1,keepdims=True),1e-15)
raw_to_label={t:f'T{i+1}' for i,t in enumerate(terms)}

# -----------------------------------------------------------------------------
# Conjugate Normal-Inverse-Gamma regression with fractional/power likelihood.
# beta | sigma2 ~ N(0, sigma2 V0); sigma2 ~ InvGamma(a0,b0).
# For terminal fits, raw terminal posterior weights are rescaled so sum(w*)
# equals Kish effective sample size: (sum w)^2/sum(w^2).
# -----------------------------------------------------------------------------
PRIOR_VAR=float(BCFG['prior_variance']); A0=float(BCFG['inverse_gamma_a0']); B0=float(BCFG['inverse_gamma_b0'])
NDRAW=int(BCFG['draws']); SUPPORT=float(BCFG['support_probability'])

def nig_fit(X, y, w=None, n_draw=NDRAW):
    return fit_bayesian(X, y, w, n_draw=n_draw, prior_variance=PRIOR_VAR,
                        a0=A0, b0=B0, rng=rng)

def design_for_config(config):
    mask=pred
    cols=[np.ones(mask.sum()),proximity[mask]]; names=['Intercept','proximity']
    if config.get('age'): cols.append(age[mask]); names.append('age')
    if config.get('sex'): cols.append(sex_M[mask]); names.append('sex_M')
    if config.get('smoking_status'):
        cols.extend([smk_current[mask],smk_quit[mask]]); names.extend(['smk_status_current','smk_status_quit'])
    if config.get('smoking_amount'): cols.append(smk_amnt[mask]); names.append('smk_amnt')
    if config.get('origin'): cols.append(origin_NCC[mask]); names.append('origin_NCC')
    npc=int(config.get('blood_pcs',0))
    for j in range(npc): cols.append(BPC[mask,j]); names.append(f'bloodPC{j+1}')
    return np.column_stack(cols),names

# Fit each declared covariate configuration; select the configured primary model.
all_rows=[]; draw_store={}; config_rows=[]
for cfg in BCFG['covariate_configs']:
    cid=cfg['id']; X,names=design_for_config(cfg); ix=names.index('proximity')
    pooled_supported=0; terminal_supported=0; terminal_supported_features=[]
    for me in ME:
        y=MEZ[me][pred]
        fit=nig_fit(X,y)
        medb=float(fit['mn'][ix]); lo=float(fit['qlo'][ix]); hi=float(fit['qhi'][ix]); pgt=float(fit['pgt'][ix]); plt=float(fit['plt'][ix]); pdir=max(pgt,plt)
        supported=bool((lo>0 or hi<0) and pdir>=SUPPORT)
        pooled_supported+=int(supported)
        all_rows.append(dict(config=cid,condition='Pooled',terminal_state=-1,terminal_label='Pooled',terminal_affinity='unit_weight',
                             feature=me,outcome=me,predictor='proximity',n=int(pred.sum()),effective_n=float(fit['effective_n']),
                             posterior_median=medb,**{'CrI_2.5':lo,'CrI_97.5':hi},Pr_gt0=pgt,Pr_lt0=plt,Pr_direction=pdir,
                             supported=supported,df=float(fit['df']),interval_method='analytic_marginal_Student_t_Normal_Inverse_Gamma_power_likelihood',
                             selection_adjusted=False))
        draw_store[f'{cid}__Pooled__{me}']=fit['draws'][:,ix].astype(np.float32)
        for j,t in enumerate(terms):
            w=Pterm[pred,j]
            fit=nig_fit(X,y,w=w)
            medb=float(fit['mn'][ix]); lo=float(fit['qlo'][ix]); hi=float(fit['qhi'][ix]); pgt=float(fit['pgt'][ix]); plt=float(fit['plt'][ix]); pdir=max(pgt,plt)
            supported=bool((lo>0 or hi<0) and pdir>=SUPPORT)
            if supported:
                terminal_supported+=1; terminal_supported_features.append(f'T{j+1}:{me}')
            all_rows.append(dict(config=cid,condition=f'T{j+1}',terminal_state=int(t),terminal_label=f'T{j+1}',
                                 terminal_affinity='raw_terminal_posterior_power_weight',feature=me,outcome=me,predictor='proximity',
                                 n=int(pred.sum()),effective_n=float(fit['effective_n']),posterior_median=medb,
                                 **{'CrI_2.5':lo,'CrI_97.5':hi},Pr_gt0=pgt,Pr_lt0=plt,Pr_direction=pdir,supported=supported,
                                 df=float(fit['df']),interval_method='analytic_marginal_Student_t_Normal_Inverse_Gamma_power_likelihood',
                                 selection_adjusted=False))
            draw_store[f'{cid}__T{j+1}__{me}']=fit['draws'][:,ix].astype(np.float32)
    complexity=(int(cfg.get('age',False))+int(cfg.get('sex',False))+int(cfg.get('smoking_status',False))*2+
                int(cfg.get('smoking_amount',False))+int(cfg.get('origin',False))+int(cfg.get('blood_pcs',0)))
    config_rows.append(dict(config=cid,pooled_supported_MEs=int(pooled_supported),terminal_supported_slopes=int(terminal_supported),
                            terminal_supported_labels=';'.join(terminal_supported_features),covariate_complexity=int(complexity),
                            age=bool(cfg.get('age')),sex=bool(cfg.get('sex')),smoking_status=bool(cfg.get('smoking_status')),
                            smoking_amount=bool(cfg.get('smoking_amount')),blood_pcs=int(cfg.get('blood_pcs',0)),origin=bool(cfg.get('origin'))))

all_df=pd.DataFrame(all_rows); audit=pd.DataFrame(config_rows)
primary=str(BCFG['primary_config'])
if primary not in set(audit.config):
    raise RuntimeError(f'Locked Bayesian ME primary_config {primary!r} is not in the declared covariate grid.')
audit['exhibits_target_pattern']=(audit.pooled_supported_MEs==0)&(audit.terminal_supported_slopes>=1)
audit['primary_config']=audit.config.eq(primary)
audit['selection_rule']=BCFG['primary_selection_rule']
audit['primary_selected_from_ME_results']=False
audit.to_csv(TAB/'ME_COVARIATE_CONFIG_AUDIT.csv',index=False)

all_df['primary_config']=all_df.config.eq(primary)
all_df['primary_selected_from_ME_results']=False
all_df.to_csv(TAB/'TERMINAL_CONDITIONED_ME_HORIZON.csv',index=False)

primary_df=all_df[all_df.primary_config].copy()
compat=primary_df[['condition','terminal_state','terminal_label','feature','effective_n','posterior_median','CrI_2.5','CrI_97.5','Pr_gt0','Pr_lt0','Pr_direction','supported']].rename(columns={'CrI_2.5':'CrI_2_5','CrI_97.5':'CrI_97_5'})
compat.to_csv(TAB/'BAYESIAN_ME_PROXIMITY_BY_TERMINAL.csv',index=False)

# Legacy filenames retained only for path compatibility; contents remain Bayesian.
audit.to_csv(TAB/'ME_CONDITIONING_MODEL_COMPARISON.csv',index=False)
primary_df.assign(legacy_filename_note='Bayesian terminal-conditioned proximity slopes; not OLS interaction coefficients').to_csv(TAB/'ME_CONDITIONING_INTERACTION_COEFFICIENTS.csv',index=False)
pd.DataFrame([dict(method='Bayesian Normal-Inverse-Gamma fractional/power-likelihood regression',
                        primary_config=primary,
                        primary_selected_from_ME_results=False,
                        pooled_supported_MEs=int(audit.loc[audit.config.eq(primary),'pooled_supported_MEs'].iloc[0]),
                        terminal_supported_slopes=int(audit.loc[audit.config.eq(primary),'terminal_supported_slopes'].iloc[0]),
                        posterior_support='95% CrI excludes 0 and posterior sign probability >=0.975',
                        FDR_or_Wald_used=False,
                        terminal_conditioning='raw frozen terminal posterior weights scaled to Kish effective sample size',
                        selection_adjusted=False,
                        note='Primary covariates exclude blood-composition PCs; the complete declared grid is retained')]).to_csv(TAB/'ME_CONDITIONING_SUMMARY.csv',index=False)

# -----------------------------------------------------------------------------
# Secondary Bayesian terminal-affinity regressions retained for existing figures.
# -----------------------------------------------------------------------------
def logit(p):
    p=np.clip(np.asarray(p,float),1e-5,1-1e-5); return np.log(p/(1-p))

def global_z(arr):
    arr=np.asarray(arr,float); medv=np.nanmedian(arr[pred]); arr=np.where(np.isfinite(arr),arr,medv)
    mu=arr[pred].mean(); sd=arr[pred].std(); return (arr-mu)/(sd+1e-12)

primary_cfg=next(c for c in BCFG['covariate_configs'] if c['id']==primary)
Xbase,nbase=design_for_config(primary_cfg)
blood_names=['Mye/Lym','Neu','Mono','CD8tot','Btot','Treg']
myl=np.log((numeric(meta.Mye)+1e-5)/(numeric(meta.Lym)+1e-5))
B=np.column_stack([global_z(myl),global_z(numeric(meta.Neu)),global_z(numeric(meta.Mono)),global_z(numeric(meta.CD8tot)),global_z(numeric(meta.Btot)),global_z(numeric(meta.Treg))])
# Avoid double-adjusting for blood PCs in the terminal-blood association model.
keep_ix=[i for i,nm in enumerate(nbase) if not nm.startswith('bloodPC')]
Xadj=Xbase[:,keep_ix]
blood_rows=[]; term_me_rows=[]
for j,t in enumerate(terms):
    y=logit(Pshare[pred,j])
    X=np.column_stack([Xadj,B[pred]]); fit=nig_fit(X,y); start=Xadj.shape[1]
    for k,nm in enumerate(blood_names):
        ix=start+k; pdir=max(float(fit['pgt'][ix]),float(fit['plt'][ix]))
        blood_rows.append(dict(terminal_state=t,terminal_label=f'T{j+1}',feature=nm,posterior_median=float(fit['mn'][ix]),
                               CrI_2_5=float(fit['qlo'][ix]),CrI_97_5=float(fit['qhi'][ix]),Pr_direction=pdir,
                               supported=bool((fit['qlo'][ix]>0 or fit['qhi'][ix]<0) and pdir>=SUPPORT),config=primary))
    Xme=np.column_stack([Xbase,np.column_stack([MEZ[m][pred] for m in ME])]); fitm=nig_fit(Xme,y); st=Xbase.shape[1]
    for k,nm in enumerate(ME):
        ix=st+k; pdir=max(float(fitm['pgt'][ix]),float(fitm['plt'][ix]))
        term_me_rows.append(dict(terminal_state=t,terminal_label=f'T{j+1}',feature=nm,posterior_median=float(fitm['mn'][ix]),
                                 CrI_2_5=float(fitm['qlo'][ix]),CrI_97_5=float(fitm['qhi'][ix]),Pr_direction=pdir,
                                 supported=bool((fitm['qlo'][ix]>0 or fitm['qhi'][ix]<0) and pdir>=SUPPORT),config=primary))
pd.DataFrame(blood_rows).to_csv(TAB/'BAYESIAN_TERMINAL_BLOOD_ASSOCIATIONS.csv',index=False)
pd.DataFrame(term_me_rows).to_csv(TAB/'BAYESIAN_TERMINAL_ME_ASSOCIATIONS.csv',index=False)

np.savez_compressed(DRAW/'TERMINAL_CONDITIONED_ME_DRAWS.npz',**draw_store)
np.savez_compressed(DRAW/'BAYESIAN_ME_PROXIMITY_DRAWS.npz',**{k:v for k,v in draw_store.items() if k.startswith(primary+'__')})
np.savez_compressed(DRAW/'TERMINAL_CONDITIONED_ME_PARAMETERS.npz',blood_pca_components=pca.components_,
                    blood_pca_explained_variance_ratio=pca.explained_variance_ratio_,terminal_states=np.array(terms),
                    blood_pca_fit_population=np.array(['prediagnostic']))
# Explicit compatibility objects; inference summaries live in CSV and primary draw object above.
np.savez_compressed(DRAW/'BAYESIAN_TERMINAL_BLOOD_DRAWS.npz',terminal_states=np.array(terms),note=np.array(['analytic Student-t posterior summaries in BAYESIAN_TERMINAL_BLOOD_ASSOCIATIONS.csv']))
np.savez_compressed(DRAW/'BAYESIAN_TERMINAL_ME_DRAWS.npz',terminal_states=np.array(terms),note=np.array(['analytic Student-t posterior summaries in BAYESIAN_TERMINAL_ME_ASSOCIATIONS.csv']))

spec=dict(result_version='3.0.0',selected_config=primary,
          formula='ME_z ~ proximity_z + declared covariates; pooled uses unit likelihood weights; each terminal-conditioned fit uses frozen raw P_state_t as a fractional likelihood weight',
          family='Conjugate Normal-Inverse-Gamma Bayesian linear regression with posterior-weighted power likelihood',
          prior='beta|sigma2 ~ Normal(0, sigma2*V0), slope prior variance 10, intercept prior variance 100; sigma2 ~ InvGamma(2,1)',
          proximity='standardized -residual_time in prediagnostic samples',interval='analytic marginal Student-t 95% credible interval',
          posterior_support_probability=SUPPORT,support_rule='95% CrI excludes 0 and posterior probability of the estimated sign >= support_probability',
          terminal_conditioning='raw frozen terminal posterior weight; weights rescaled so total likelihood weight equals Kish effective sample size',
          terminal_state_labels={str(t):f'T{i+1}' for i,t in enumerate(terms)},covariate_configurations=BCFG['covariate_configs'],
          primary_selection_rule=BCFG['primary_selection_rule'],selection_adjusted=False,
          model_scope='Pointwise posterior intervals across the declared model configurations',
          blood_pca_fitted_on='408 prediagnostic source samples; 12-cell CLR after zero replacement',blood_pca_components=3,
          zero_replacement=pc,draw_count_per_fit=NDRAW,table='analysis/tables/TERMINAL_CONDITIONED_ME_HORIZON.csv',
          draws='cache/posterior_draws/TERMINAL_CONDITIONED_ME_DRAWS.npz')
json.dump(spec,open(SPEC/'TERMINAL_CONDITIONED_ME_MODEL.json','w'),indent=2)
json.dump(dict(primary_config=primary,selection_rule=BCFG['primary_selection_rule'],primary_selected_from_ME_results=False,
               development_phase_exploratory=True,FDR_or_Wald_used=False,
               interval_multiplicity_adjusted=False,terminal_conditioning='soft raw terminal posterior power likelihood'),
          open(TAB/'ME_PRIMARY_CONFIG.json','w'),indent=2)

print('primary ME covariate config:',primary)
print(audit.to_string(index=False))
print('\nPrimary supported slopes:')
print(primary_df[primary_df.supported][['condition','feature','effective_n','posterior_median','CrI_2.5','CrI_97.5','Pr_direction']].to_string(index=False))
