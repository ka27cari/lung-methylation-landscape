# Export figure-associated numerical tables from local, reduced R objects.
# Usage: Rscript r/source_data/export_figure_rds.R <figure1.rds> <figure2.rds> <output_directory>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop('Provide: figure1_session_summary.rds figure2_minimal_review.rds output_directory')
}

out <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

save_columns <- function(x, columns, name, subdir) {
  if (!is.data.frame(x) || !all(columns %in% names(x))) {
    stop('Unexpected columns in ', name)
  }
  dat <- x[, columns, drop = FALSE]
  if (any(vapply(dat, is.list, logical(1)))) stop('Nested columns in ', name)
  dir.create(file.path(out, subdir), showWarnings = FALSE, recursive = TRUE)
  write.csv(dat, file.path(out, subdir, paste0(name, '.csv')),
            row.names = FALSE, na = '', fileEncoding = 'UTF-8')
}

f1 <- readRDS(args[[1]])
f2 <- readRDS(args[[2]])
if (!is.list(f1) || !is.list(f2) || !is.list(f2$fig2)) {
  stop('Expected reduced figure RDS objects')
}

save_columns(f1$risk_surface,
             c('smk_dur', 'smk_amnt', 'mu', 'lower', 'age_band'),
             'risk_surface', 'Fig1')
save_columns(f1$survival_halfyear_grid,
             c('group', 'time_years', 'survival'),
             'survival_curves', 'Fig1')

fields <- list(
  ewas_selected_cpg = c('cpg', 'logFC', 'P.Value', 'adj.P.Val', 'bonferroni'),
  ewas_counts = c('tested', 'nominal_p05', 'bonferroni_p05'),
  roc_coordinates = c('model', 'fpr', 'tpr'),
  model_performance = c('model', 'ro_auc', 'pr_auc', 'n', 'n_used'),
  histology_recall = c('type', 'n', 'n_correct', 'n_incorrect', 'recall'),
  module_cpg_map = c('cpg', 'module')
)
for (name in names(fields)) {
  save_columns(f2$fig2[[name]], fields[[name]], name, 'Fig2')
}

cat('Wrote aggregate figure tables to:', out, '\n')
