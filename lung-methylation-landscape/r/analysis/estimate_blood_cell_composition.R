# 1. Set Environment ───────────────────────────────────────────────────────────
library(EpiDISH)
library(ExperimentHub)

# Load source
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Load beta matrix
beta <- file.path(data_dir, "beta-raw.rds") %>% read_rds; rownames(beta) %<>% strip_sfx

# Control conflicts
conflicts_prefer(base::intersect)

# 2. Estimate Cell Heterogeneity ───────────────────────────────────────────────
# Search references
eh <- ExperimentHub()
qh <- query(eh, c("FlowSorted", "Blood", "EPIC"))
mcols(qh) %>% as_tibble %>% View

# Load a reference
ref <- eh[["EH5425"]]

# Make reference beta matrix
ref_beta <- getBeta(preprocessNoob(ref))

ct  <- as.character(colData(ref)$CellType)
set <- split(seq_len(ncol(ref_beta)), ct)

# Remove artifact
set <- set[names(set) != "MIX"]

mat <- sapply(set, \(i) rowMeans(ref_beta[, i, drop = FALSE], na.rm = TRUE))
tgt <- intersect(rownames(mat), rownames(beta))

na_rate <- rowMeans(is.na(beta[tgt, , drop = FALSE]))  # check NA rates
keep    <- tgt[na_rate <= 0.3]

# Estimate heterogeneity
hetg <- epidish(beta.m = beta[keep, ], ref.m = mat[keep, ], method = "RPC")$estF
hetg %>% as.data.frame %>% rownames_to_column("sentrix_id") %>% write_rds("object/heterogeneity.rds")
