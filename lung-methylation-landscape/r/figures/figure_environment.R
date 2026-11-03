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
library(magrittr)
library(ggdendro)
library(ggplotify)
library(patchwork)
library(tidyverse)
library(colorspace)
library(conflicted)
library(ggnewscale)
library(geomtextpath)

conflicted::conflicts_prefer(dplyr::filter)

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

# 3. Study paths
base <- Sys.getenv("LUNG_METHYLATION_PROJECT_ROOT", unset = "")
if (!nzchar(base)) stop("Set LUNG_METHYLATION_PROJECT_ROOT to the study directory")
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
setwd(file.path(base, "05_ncc"))
if (!dir.exists("result")) dir.create("result")

manifest_dir <- file.path(base, "00_manifest")
pheno_dir    <- file.path(base, "01_phenotype")
data_dir     <- file.path(base, "02_preprocess")
build_dir    <- file.path(Sys.getenv("LUNG_METHYLATION_REPO"), "rebuild")

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ NA_character_,
    p < 0.001   ~ "***",
    p < 0.01    ~ "**",
    p < 0.05    ~ "*",
    TRUE        ~ "n.s"
  )
}

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
if (.Platform$OS.type == "windows") {
  windowsFonts(
    Arial = windowsFont("Arial"),
    Calibri = windowsFont("Calibri"),
    Times = windowsFont("Times New Roman")
  )
}


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
  legend.spacing.y = unit(2, "mm"), 
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

library(svglite)
