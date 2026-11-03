library(WGCNA)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

beta <- file.path(data_dir, "beta-raw.rds") %>% read_rds
net  <- read_rds("net-modules.rds")

# Redefine eigengenes
ids  <- colnames(beta); cpgs <- names(net$colors)
expr <- beta[cpgs, ] %>% logit %>% impute %>% t

egs <- moduleEigengenes(expr, colors = net$colors)$eigengenes
egs$sentrix_id <- ids

egs %<>% as_tibble %>% relocate(sentrix_id, .before = 0)

meta %<>% lapply(left_join, egs)

# Save
meta$nested %>% write_csv("model/data/nested.csv")
meta$longitudinal %>% write_csv("model/data/longitudinal.csv")

# Load model data
load("modeling-set.rda")

m_tr %>% as.data.frame %>% rownames_to_column("cpg") %>% write_csv("model/data/m_tr.csv")
m_te %>% as.data.frame %>% rownames_to_column("cpg") %>% write_csv("model/data/m_te.csv")

read_csv("result/table/aggregate.csv") %>% write_csv("model/data/agg.csv")

read_rds("cg-set.rds") %>% imap_dfr(
  \(x, nm) {
    imap_dfr(
      x, \(cpgs, axis) {
        tibble(
          set  = nm,
          axis = axis,
          cpg  = cpgs
        ) 
      }
    )
  }
) %>% write_csv(
  "model/data/cg-set.csv"
)

# Load external data
load("external-data.rda")

ph_ext %>% write_csv("model/data/ph-ext.csv")
m_ext %>% as.data.frame %>% rownames_to_column("cpg") %>% write_csv("model/data/m-ext.csv")

# Load trajectory
read_csv("result/table/cluster-assigned.csv") %>% write_csv("model/data/cluster-assigned.csv")
read_csv("result/table/cluster-summary.csv") %>% write_csv("model/data/cluster-summary.csv")
read_csv("result/table/cluster-spec.csv") %>% write_csv("model/data/cluster-spec.csv")
