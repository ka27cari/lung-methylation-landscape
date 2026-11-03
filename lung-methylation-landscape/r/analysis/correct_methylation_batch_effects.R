# 1. Set Environment ───────────────────────────────────────────────────────────
library(umap)
library(minfi)
library(ComBatMet)
library(missMethyl)

# Load source
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Load data
beta <- file.path(data_dir, "beta-raw.rds") %>% read_rds
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Control conflicts
conflicts_prefer(base::intersect, matrixStats::rowVars, matrixStats::rowMedians)

# 2. Generate Functions ────────────────────────────────────────────────────────
# Dimensional reduction on partial matrix
partial <- function(
    mat, meta, anno,
    n_pcs = 50, n_neighbors = 10, n_target = 100, frac_min = 0.70, seed = NULL
) {
  idx <- is_auto(anno$CHR)  # keep autosomes only
  tgt <- intersect(rownames(mat), anno$IlmnID[idx])
  m   <- logit(mat[tgt, meta$sentrix_id, drop = FALSE])
  
  # NA rates by origin
  ok <- lapply(
    unique(meta$origin), \(i) {
      idx    <- which(meta$origin == i)
      non_na <- rowMeans(is.finite(m[, idx, drop = FALSE]))
      non_na >= frac_min
    }
  )

  keep <- Reduce("&", ok)
  m    <- impute(m[keep, , drop = FALSE])
  
  # Keep highly variable CpGs
  if (nrow(m) > n_target) {
    var <- rowVars(m, na.rm = TRUE)
    m   <- m[order(var, decreasing = TRUE)[seq_len(n_target)], , drop = FALSE]
  }

  sc <- t(scale(t(m), center = TRUE, scale = TRUE))

  # Do UMAP based on PCs
  pca <- prcomp(t(sc), center = TRUE, scale. = FALSE)
  pcs <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]

  if (!is.null(seed)) set.seed(seed)
  umap <- umap(pcs, n_components = 2, n_neighbors = n_neighbors, metric = "euclidean")

  return(umap)
}

# 3. Batch Correction ──────────────────────────────────────────────────────────
coh <- meta$nested %>% mutate(across(starts_with("date"), as.character))
coh$ctr_age <- scale(coh$age, center = TRUE, scale = FALSE)

# Run embedding on raw data
bef <- partial(meta = coh, mat = beta, anno = anno, n_target = 1e3)

# Assign batches
coh$batch <- case_when(
  coh$date_analyzed %in% c("2023-03-29", "2023-04-25") ~ 1L,
  coh$date_analyzed %in% c("2025-07-24", "2025-07-30") ~ 2L,
  coh$date_analyzed %in% c("2025-08-01", "2025-08-06", "2025-08-07", "2025-08-08") ~ 3L,
  .default = 0L
)

# Prepare model matrix for adjustment
mod <- model.matrix(
  ~ sex + ctr_age + smoking_status + Bmem + Bnv + CD4mem + CD4nv + Treg + CD8mem + CD8nv + Bas + Eos + Mono + Neu, data = coh
)

# Apply ComBat-met()
new_beta <- ComBat_met(
  bv        = beta[, coh$sentrix_id],
  batch     = coh$batch,
  group     = coh$histology,  # not dissort gradual drift too much
  covar_mod = mod,
  full_mod  = TRUE
)

# Run embedding on adjusted data
aft <- partial(meta = coh, mat = new_beta, anno = anno, n_target = 1e3)

emb <- list(before = bef, after = aft)

# Save objects
write_rds(emb, "object/embed.rds")
file.path(data_dir, "beta-adjusted.rds") %>% write_rds(new_beta, .)
