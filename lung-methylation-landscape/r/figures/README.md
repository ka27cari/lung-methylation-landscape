# Figure-generation scripts

These scripts create component plots, not the final arranged manuscript figures.
Each file is named after its analysis topic, not a historical execution number.

| Script | What it draws |
|---|---|
| `methylation_programs.R` | Recurrent programs, chromatin annotation, embeddings |
| `microstate_landscape.R` | Microstate geometry and program profiles |
| `state_transport_routes.R` | Transport currents and route maps |
| `terminal_blood_associations.R` | Terminal associations and out-of-fold comparisons |
| `longitudinal_dynamics.R` | Paired changes and routes over time |
| `external_pbmc_validation.R` | External PBMC projection |
| `supplementary_figures.R` | Supplementary plots |
| `figure_environment.R` | Shared setup and project-local object loading |

For a study-specific run, configure `LUNG_METHYLATION_PROJECT_ROOT`
to the private study directory containing `05_ncc/`, and
`LUNG_METHYLATION_REPO` to the cloned code repository.
For example, in an R session with the approved R packages and private objects:

```r
Sys.setenv(LUNG_METHYLATION_PROJECT_ROOT = "C:/path/to/Project")
Sys.setenv(LUNG_METHYLATION_REPO = "C:/path/to/lung-methylation-landscape")
source(file.path(Sys.getenv("LUNG_METHYLATION_REPO"), "r", "figures", "methylation_programs.R"))
```

`figure_environment.R` reads private `05_ncc/object/*.rds` objects,
which are not in this public repository. These R scripts have **not**
been run end-to-end in the packaging environment. The Python pipeline
is independently tested; an R rendering validation must use the
authorized study environment.
