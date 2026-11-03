# File rename map

Names describe scientific purpose; analytical formulas remain unchanged.

## Python stages (execution order is in run.py)

| Previous | New |
|---|---|
| `pipeline/maintenance/00_clean_outputs.py` | `pipeline/clean_generated_outputs.py` |
| `pipeline/discovery/01_anchor_resolution.py` | `pipeline/resolve_reference_anchors.py` |
| `pipeline/discovery/02_discover_programs.py` | `pipeline/discover_recurrent_programs.py` |
| `pipeline/discovery/03_select_programs.py` | `pipeline/select_stable_programs.py` |
| `pipeline/landscape/04_fit_landscape.py` | `pipeline/fit_microstate_landscape.py` |
| `pipeline/landscape/05_calibrate_posterior.py` | `pipeline/calibrate_terminal_posterior.py` |
| `pipeline/landscape/05b_sinkhorn_net_current.py` | `pipeline/infer_sinkhorn_transport.py` |
| `pipeline/validation/06_downstream.py` | `pipeline/analyze_landscape_outcomes.py` |
| `pipeline/validation/07_context_models.py` | `pipeline/fit_terminal_blood_associations.py` |
| `pipeline/validation/12_finalize_extensions.py` | `pipeline/summarize_validation_extensions.py` |
| `pipeline/validation/12a_incremental_oof.py` | `pipeline/compare_out_of_fold_models.py` |
| `pipeline/reporting/08_export_figure_data.py` | `pipeline/export_figure_tables.py` |
| `pipeline/reporting/09_render_core_figures.py` | `pipeline/render_analysis_diagnostics.py` |
| `pipeline/reporting/10_validate_build.py` | `pipeline/check_analysis_outputs.py` |

## R study analyses

| Previous | New |
|---|---|
| `r/preprocessing/cohorts/00_environment.r` | `r/analysis/analysis_environment.R` |
| `r/preprocessing/cohorts/01_manipulation.r` | `r/analysis/prepare_methylation_cohorts.R` |
| `r/preprocessing/cohorts/02_heterogeniety.r` | `r/analysis/estimate_blood_cell_composition.R` |
| `r/preprocessing/cohorts/03_batch-correction.r` | `r/analysis/correct_methylation_batch_effects.R` |
| `r/preprocessing/cohorts/04_cohort.r` | `r/analysis/analyze_cohort_risk_and_screening.R` |
| `r/preprocessing/cohorts/05_nested.r` | `r/analysis/analyze_nested_case_control.R` |
| `r/preprocessing/epigenome/06_reproduction.r` | `r/analysis/analyze_epigenomic_replication.R` |
| `r/preprocessing/epigenome/07_validation.r` | `r/analysis/validate_methylation_predictors.R` |
| `r/preprocessing/epigenome/08_ewas.r` | `r/analysis/run_epigenome_wide_associations.R` |
| `r/preprocessing/epigenome/09_inflection.r` | `r/analysis/analyze_diagnosis_proximity.R` |
| `r/preprocessing/summaries/11_interim.r` | `r/analysis/summarize_methylation_networks.R` |
| `r/preprocessing/summaries/12_aov.r` | `r/analysis/analyze_program_cell_type_context.R` |

## R figure scripts

| Previous | New |
|---|---|
| `r/plots/01_programs.R` | `r/figures/methylation_programs.R` |
| `r/plots/02_states.R` | `r/figures/microstate_landscape.R` |
| `r/plots/03_transport.R` | `r/figures/state_transport_routes.R` |
| `r/plots/04_associations.R` | `r/figures/terminal_blood_associations.R` |
| `r/plots/05_longitudinal.R` | `r/figures/longitudinal_dynamics.R` |
| `r/plots/06_external.R` | `r/figures/external_pbmc_validation.R` |
| `r/plots/07_supplementary.R` | `r/figures/supplementary_figures.R` |
| `r/plots/common.R` | `r/figures/figure_environment.R` |
