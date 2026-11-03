# 1. Load Libraries ────────────────────────────────────────────────────────────
library(broom)
library(ggsci)
library(ggpubr)
library(ggside)
library(scales)
library(readxl)
library(useful)
library(cowplot)
library(forcats)
library(ggrepel)
library(ggforce)
library(pbapply)
library(svglite)
library(magrittr)
library(ggdendro)
library(ggplotify)
library(patchwork)
library(tidyverse)
library(colorspace)
library(conflicted)
library(ggnewscale)
library(geomtextpath)

# Set library priority
conflict_prefer_all("dplyr")
conflict_prefer_all("matrixStats")

# Control conflicts
conflicts_prefer(purrr::discard)

# 2. Generate Functions ────────────────────────────────────────────────────────
# Get root of analysis environment
get_root <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.getenv("USERPROFILE")
  } else {
    Sys.getenv("HOME")
  }
}

# Find working path
find_path <- function(name, max_depth = 2) {
  name <- paste0("/", name, "$")
  
  root <- get_root()  # get root
  curr <- list(root)  # level 0
  all_dirs <- character()
  
  for (depth in seq_len(max_depth)) {
    next_dir <- list()
    for (dir in curr) {
      sub_dirs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
      next_dir <- c(next_dir, sub_dirs)
    }
    all_dirs <- c(all_dirs, next_dir)
    curr <- next_dir
  }
  
  # All directories under max_depth
  candidates <- all_dirs
  paths <- grep(name, candidates, value = TRUE)
  
  if (length(paths) == 0) {
    message("Cannot find the path")
    return(NULL)
  } else {
    return(paths)
  }
}

# Strip suffix
strip_sfx <- function(x) sub("_.*$", "", x)

# Find suffix
add_sfx <- function(x) {
  sfx     <- anno$Name
  cpg_map <- sub("_(TC|BC)\\d{2}$", "", sfx)
  
  sfx %<>% setNames(cpg_map)
  sfx[x] %>% unname
}

# Logit transformation
logit <- function(x, eps = 1e-6, thr = 4, block = 25) {
  if (as.numeric(object.size(x)) <= thr * 1024 ^ 3) {
    x <- pmin(pmax(x, eps), 1 - eps)
    return(qlogis(x))
  }
  
  if (is.null(dim(x))) {
    x <- pmin(pmax(x, eps), 1 - eps)
    return(qlogis(x))
  }
  
  p <- dim(x)[2]
  
  for (j in seq.int(1L, p, by = block)) {
    idx <- j:min(j + block - 1L, p)
    b <- x[, idx, drop = FALSE]
    b <- pmin(pmax(b, eps), 1 - eps)
    
    x[, idx] <- qlogis(b)
    rm(b); gc()
  }
  
  x
}

# Inverse logit transformation
ilogit <- function(vec) {
  vec <- pmin(pmax(vec, -30), 30)  # prevent overflow
  exp(vec) / (1 + exp(vec))
}

# Split specimens
divide <- function(
    data, base = NULL, ...
) {
  cond_qs <- enquos(...)
  base_q  <- enquo(base)
  
  imap(
    cond_qs, function(cond, nm) {
      out <- data
      
      if (!rlang::quo_is_null(base_q)) {
        out <- filter(out, !!base_q)
      }
      filter(out, !!cond)
    }
  )
}

# Check if autosomal CpG exists in annotation table
is_auto <- function(chr) sub("^chr","", chr) %in% as.character(1:22)

# Impute matrix
impute <- function(mat) {
  require(matrixStats)
  
  mat %<>% as.matrix
  mat[is.infinite(mat)] <- NA
  
  # Get median for each row
  med <- rowMedians(mat, na.rm = TRUE)
  
  keep <- is.finite(med)
  mat  <- mat[keep, , drop = FALSE]
  med  <- med[keep]
  
  na_pos <- which(is.na(mat), arr.ind = TRUE)
  if (length(na_pos)) mat[na_pos] <- med[na_pos[, 1]]
  
  # Rule out zero variance CpG
  v <- rowVars(mat)
  mat[v > 0, , drop = FALSE]
}

# Convert significance to asterisk
p_to_star <- function(p) {
  star <- "*"
  if (is.na(p)) "n.s"
  else if (p < 0.001) paste0(star, star, star)
  else if (p < 0.01)  paste0(star, star)
  else if (p < 0.05)  star
  else if (p < 0.1)   "\u0387"
  else "n.s"
}

# Fast AUC
fast_auc <- function(score, y) {
  y <- as.integer(y)
  
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# AUC Checker
eval_auc <- function(pos, neg) {
  score <- c(pos, neg)
  y     <- c(rep(1L, length(pos)), rep(0L, length(neg)))
  
  fast_auc(score, y)
}

# 3. Set File Paths ────────────────────────────────────────────────────────────
# Get working directory path
base <- Sys.getenv("LUNG_METHYLATION_PROJECT_ROOT", unset = "")
if (!nzchar(base)) stop("Set LUNG_METHYLATION_PROJECT_ROOT to the study directory")
setwd(file.path(base, "05_ncc/"))

manifest_dir <- file.path(base, "00_manifest")
pheno_dir    <- file.path(base, "01_phenotype")
data_dir     <- file.path(base, "02_preprocess")
build_dir    <- file.path(Sys.getenv("LUNG_METHYLATION_REPO"), "rebuild")

# 4. Load Dataset ──────────────────────────────────────────────────────────────
meta <- read_rds("object/meta.rds")
hetg <- read_rds("object/heterogeneity.rds")

meta %<>% map(
  ~ {
    .x %<>% left_join(
      hetg
    ) %>% mutate(
      CD4tot = CD4mem + CD4nv,
      CD8tot = CD8mem + CD8nv,
      Btot   = Bmem + Bnv,
      Mye    = Neu + Mono + Eos + Bas,
      Lym    = CD4tot + CD8tot + Btot + NK + Treg
    )
    
    # Add histological types
    pre <- .x %>% filter(
      group == "prediagnostic"
    ) %>% mutate(
      type = if_else(
        morphology %in% c("Large cell carcinoma", "Neuroendocrine carcinoma", "NOS, unspecified"),
        "Other (rare, unspecific)", 
        morphology
      ),
      time_bin = cut(
        time,
        breaks = c(-Inf, 2, 5, 10, Inf),
        labels = c("<2", "2-5", "5-10", ">10"),
        right  = TRUE
      )
    )
    
    pre$type %<>% fct_recode(
      "LUAD"  = "Adenocarcinoma (incl. variants)",
      "LUSC"  = "Squamous cell carcinoma",
      "SCLC"  = "Small cell carcinoma",
      "NSCLC" = "Non-small cell lung cancer",
      "Other" = "Other (rare, unspecific)"
    )
    
    pre %<>% select(matched_index, group, type) %>% distinct
    
    .x %>% left_join(pre) %>% drop_na(sentrix_id)
  }
)

# 5. Visualization Set-up ──────────────────────────────────────────────────────
windowsFonts(
  Arial   = windowsFont("Arial"),
  Calibri = windowsFont("Calibri"),
  Times   = windowsFont("Times New Roman")
)

names(windowsFonts())

# Basic theme
basic <- theme_classic(
  base_size = 9
) + theme(
  strip.background = element_blank(),
  strip.text = element_text(size = 6),
  
  axis.line = element_line(linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.35),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 7),
  
  plot.margin = margin(l = 10, r = 5, b = 5, t = 5),
  
  plot.tag = element_text(
    family = "Arial", face = "bold", size = 14,
    margin = margin(l = -5)
  ),
  plot.tag.location = "plot",
  plot.tag.position = c(0, 0.95),
  
  legend.title.align = 0,
  legend.title = element_text(
    size = 6, margin = margin(b = 1.5)
  ),
  legend.spacing.y = unit(5, "mm"), 
  legend.key.spacing.y = unit(0.5, "mm"),
  legend.key.height = unit(2, "mm"),
  legend.key.width = unit(2, "mm"),
  legend.text = element_text(
    size = 6,
    margin = margin(l = 2)
  ),
  legend.box.spacing = unit(3, "mm"),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.key = element_rect(fill = "transparent", color = NA)
)
