# Cohort and epigenome R analyses

The original study-specific cohort scripts are grouped in this **one directory**
with descriptive names rather than multiple numbered folders.
They are retained for traceability and use private cohort/project assets.
They are **not** invoked by `python run.py run`.

| Script | Purpose |
|---|---|
| `analysis_environment.R` | R libraries, analysis functions, project paths |
| `prepare_methylation_cohorts.R` | Assemble and preprocess cohort methylation |
| `estimate_blood_cell_composition.R` | Estimate inferred cell composition |
| `correct_methylation_batch_effects.R` | Handle batch effects |
| `analyze_cohort_risk_and_screening.R` | Cohort risk, LDCT and epidemiology |
| `analyze_nested_case_control.R` | Nested case-control comparisons |
| `analyze_epigenomic_replication.R` | Epigenome and module replication |
| `validate_methylation_predictors.R` | Methylation marker validation |
| `run_epigenome_wide_associations.R` | EWAS |
| `analyze_diagnosis_proximity.R` | Time-to-diagnosis patterns |
| `summarize_methylation_networks.R` | Module summaries |
| `analyze_program_cell_type_context.R` | External program and cell-type context |

Set `LUNG_METHYLATION_REPO` and `LUNG_METHYLATION_PROJECT_ROOT`
before running a script; study-specific private files and R packages are required.
Legacy relative paths to private study resources remain unchanged.
