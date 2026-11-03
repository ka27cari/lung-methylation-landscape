# 1. Set Environment ───────────────────────────────────────────────────────────
library(limma)
library(WGCNA)
library(rGREAT)
library(DMRcate)
library(msigdbr)
library(data.table)
library(missMethyl)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Set priority
conflicts_prefer(base::intersect, stats::cor)

# Load data
hr <- read_rds("object/hazard-risk.rds")

meta %<>% map2(hr, left_join)

# Expression matrix
beta <- file.path(data_dir, "beta-adjusted.rds") %>% read_rds; obs <- colnames(beta)
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# eForge
e_forge <- read_tsv("result/table/h3.tsv")

# Clonal haematopoiesis
esm <- file.path(manifest_dir, "41467_2025_59333_MOESM34_ESM.xlsx") %>% read_xlsx

# 2. Generate Functions ────────────────────────────────────────────────────────
# Correlation p-value
cor_p <- function(r, n) {
  if (is.na(r) || is.na(n) || n < 3) return(NA)
  t <- r * sqrt((n - 2) / pmax(1 - r ^ 2, .Machine$double.eps))
  2 * pt(-abs(t), df = n - 2)
}

# Fisher (one-sided greater) with Haldane OR
fisher <- function(query, bg, chip) {
  query <- intersect(unique(query), unique(bg))
  out   <- setdiff(unique(bg), query)
  
  a <- sum(query %in% chip)
  b <- length(query) - a
  c <- sum(out %in% chip)
  d <- length(out) - c
  
  mat <- matrix(c(a, b, c, d), 2, byrow = TRUE)
  
  p <- tryCatch(
    fisher.test(mat, alternative = "greater")$p.value, 
    error = function(e) NA
  )
  
  or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  tibble(a = a, b = b, c = c, d = d, or = or, p = p)
}

# 3. Risk Exploration of CpGs ──────────────────────────────────────────────────
nst <- meta$nested %>% filter(
  sentrix_id %in% obs
) %>% mutate(
  event    = if_else(group != "control", 1L, 0L),
  entry    = age,
  exit     = age + time,
  ever     = smk_status != "none",
  past     = smk_status == "quit",
  current  = smk_status == "current",
  smk_dur  = case_when(
    isFALSE(ever)   ~ 0, 
    !is.na(smk_dur) ~ smk_dur, 
    !is.na(smk_age) ~ pmax(0, age - smk_age), 
    TRUE ~ 0
  ),
  ces_dur  = case_when(
    isFALSE(past)   ~ 0,
    !is.na(ces_dur) ~ ces_dur,
    !is.na(ces_age) ~ pmax(0, age - ces_age),
    TRUE ~ 0
  ),
  na_amnt  = ever & is.na(smk_amnt),
  packyear = (smk_amnt / 20) * smk_dur,
  screen   = case_when(
    is.na(age) | is.na(smk_status) | is.na(smk_amnt) | is.na(smk_dur) ~ NA,
    age < 54 | age > 74                   ~ FALSE,
    packyear < 30                         ~ FALSE,
    smk_status == "none"                  ~ FALSE,
    smk_status == "current"               ~ TRUE,
    smk_status == "quit" & is.na(ces_dur) ~ NA,
    smk_status == "quit" & ces_dur > 15   ~ FALSE,
    smk_status == "quit" & ces_dur <= 15  ~ TRUE,
    TRUE ~ NA
  )
)

m <- beta[, nst$sentrix_id, drop = FALSE] %>% logit

# Filter CpGs if NA > 30%
na <- rowMeans(is.na(m), na.rm = TRUE)
m  <- m[na < 0.30, , drop = FALSE]

cpgs <- intersect(rownames(m), anno$Name)

m <- m[cpgs, , drop = FALSE]

# Save
file.path(data_dir, "m-adjusted.rds") %>% write_rds(m, .)

# EWAS in NCC examinee
ncc  <- nst %>% filter(origin != "KCPS")
expr <- m[, ncc$sentrix_id, drop = FALSE]

# Design matrix
design <- model.matrix(~ risk_score + Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + Neu + NK, ncc)

fit <- lmFit(expr, design)
fit %<>% eBayes

# Extract significance
ext <- topTable(fit, coef = "risk_score", number = Inf, sort.by = "p", adjust.method = "fdr")
ext %<>% rownames_to_column("cpg")
ext$bonferroni <- p.adjust(ext$P.Value, "bonferroni")

# Save
write_rds(ext, "object/ewas-risk.rds")
ext <- read_rds("object/ewas-risk.rds")

# Annotate CpGs
ann <- cpg.annotate(
  datatype      = "array",
  object        = expr,
  what          = "M",
  arraytype     = "EPICv2",
  analysis.type = "differential",
  design        = design,
  coef          = "risk_score",
  fdr           = 0.05,
  epicv2Filter  = "mean"
)

# Save
write_rds(ann, "object/cpg-annot.rds")
ann <- read_rds("object/cpg-annot.rds")

# DMR calling
dmr <- dmrcate(ann, lambda = 1000, C = 2)
rng <- extractRanges(dmr, genome = "hg38")

# Match extract with annotation
anno %<>% as.data.frame %>% transmute(
  cpg   = Name,
  base  = sub("_.*$", "", Name),
  chr   = chr,
  pos   = pos,
  gene  = UCSC_RefGene_Name,
  group = UCSC_RefGene_Group,
  prev  = Methyl450_Loci
)

ext %<>% left_join(anno, by = "cpg") %>% drop_na(pos)

ext$base <- sub("_.*$", "", ext$cpg)
ext$chr %<>% fct_relevel(paste0("chr", c(1:22, "X", "Y")))

# Match extract to gene set
sg <- ext %>% filter(bonferroni < 0.05)

sg_gr <- GRanges(
  seqnames = sg$chr,
  ranges   = IRanges(sg$pos, width = 1),
  cpg      = sg$cpg,
  p        = sg$P.Value
)

hits <- findOverlaps(sg_gr, rng)
rng  <- rng[sort(unique(subjectHits(hits)))]

# Reactome
reactome <- great(rng, "msigdb:C2:CP:REACTOME", "hg38")

# Top Reactome terms
tr <- getEnrichmentTable(reactome)
tr %<>% filter(p_adjust < 0.05)

tr$label <- tr$id %>% str_remove("^REACTOME_") %>% str_replace_all("_", " ") %>% str_to_lower %>% str_to_sentence

# Replace nomenclatures
tr$label %<>% str_replace_all(
  "\\bDna\\b", "DNA"
) %>% str_replace_all(
  "\\bdna\\b", "DNA"
) %>% str_replace_all(
  "\\bRna\\b", "RNA"
) %>% str_replace_all(
  "\\bMapk\\b", "MAPK"
) %>% str_replace_all(
  "\\bmapk\\b", "MAPK"
) %>% str_replace_all(
  "\\bMap2k\\b", "MAP2K"
) %>% str_replace_all(
  "\\bBraf\\b", "BRAF"
) %>% str_replace_all(
  "\\bRaf1\\b", "RAF1"
) %>% str_replace_all(
  "\\bExtracellular matrix\\b", "extracellular matrix"
) %>% str_replace_all(
  "\\by family\\b", "Y-family"
) %>% str_replace_all(
  "\\bbraf\\b", "BRAF"
) %>% str_replace_all(
  "\\braf1\\b", "RAF1"
)

tr$label %<>% recode(
  "Translesion synthesis by Y-family DNA polymerases bypasses lesions on DNA template" =
    "Translesion synthesis by Y-family DNA polymerases",
  "Metabolism of lipids" =
    "Lipid metabolism",
  "Signaling by BRAF and RAF1 fusions" =
    "BRAF/RAF1 fusion signaling",
  "MAP2K and MAPK activation" =
    "MAP2K/MAPK activation"
)

# Cumulative position for Manhattan
ext$num <- match(ext$chr, paste0("chr", c(1:22, "X", "Y")))

info <- ext %>% group_by(num, chr) %>% summarise(len = max(pos, na.rm = TRUE), .groups = "drop") %>% arrange(num)
info$offset <- lag(cumsum(as.numeric(info$len)), default = 0)

ext %<>% left_join(
  info, by = c("num", "chr")
) %>% mutate(
  c_pos = pos + offset,
  log_p = -log10(P.Value)
)

axis <- info %>% mutate(center = offset + len / 2)

# Connect CpGs to GRanges
ext_gr <- GRanges(
  seqnames = ext$chr,
  ranges   = IRanges(ext$pos, width = 1),
  cpg      = ext$cpg,
  P.Value  = ext$P.Value
)

# Find representative term by each annotation points
n <- nrow(tr)

lab <- map_dfr(
  1:n, function(i) {
    id    <- tr$id[i]
    label <- tr$label[i]
    
    assoc <- getRegionGeneAssociations(reactome, term_id = id)
    
    hit <- findOverlaps(ext_gr, assoc)
    idx <- sort(unique(queryHits(hit)))
    
    cand <- ext[idx, , drop = FALSE] %>% slice_min(
      P.Value, n = 1, with_ties = FALSE
    ) %>% mutate(
      id    = id,
      label = label
    )
    
    return(cand)
  }
)

# Position adjustment
spec <- tribble(
  ~label,                                              ~shift, ~y_txt, ~off, ~length,
  "Deubiquitination",                                  17e7,   34,     NA,   NA,
  "Innate immune system",                              23e7,    30,     1.2,  10,
  "Cytokine signaling in immune system",               1e8,    23,     1.2,  3,
  "Lipid metabolism",                                  19e7,   32,     1.2,  4,
  "Intrinsic pathway for apoptosis",                   3e8,   27,     1.2,  7,
  "DNA repair",                                        11e7,   23,     1.2,  3,
  "Interferon signaling",                              19e7,   28,     1.2,  8,
  "Cell extracellular matrix interactions",            37e7,   33,     1.2,  13,
  "BRAF/RAF1 fusion signaling",                        29e7,   23,     1.2,  3,
  "Oncogenic MAPK signaling",                          28e7,   27,     NA,   NA,
  "Translesion synthesis by Y-family DNA polymerases", NA,     NA,     NA,   NA,
  "DNA damage bypass",                                 24e7,   33,     1.2,  13,
)

pos <- lab %>% arrange(c_pos) %>% inner_join(
  spec, by = "label"
) %>% mutate(
  x_txt     = c_pos - shift,
  seg_start = y_txt - off,
  seg_end   = seg_start - length
) %>% select(
  -shift, -off, -length
)

# Pick frequent genes
d <- sg %>% transmute(
  cpg, chr, pos,
  p_value = P.Value,
  gene    = strsplit(gene, ";")
) %>% unnest(
  gene
) %>% mutate(
  gene = trimws(gene)
)

frq <- d %>% dplyr::count(gene, sort = TRUE) %>% filter(n > 5) %>% pull(gene)
d %<>% filter(gene %in% frq) %>% distinct

gn <- ext %>% filter(cpg %in% d$cpg) %>% select(cpg, c_pos, log_p) %>% left_join(d)
gn %<>% filter(chr != "chrX") %>% select(c_pos, log_p, gene) %>% group_by(gene) %>% slice_max(log_p) %>% ungroup

# Down-sampling
fg <- ext %>% filter(bonferroni < 0.05)
bg <- anti_join(
  ext, fg
) %>% mutate(
  bin = cut(
    log_p,
    breaks = c(-Inf, 1, 2, 3, 4, Inf),
    labels = c("<=1", "1~2", "2~3", "3~4", ">=4")
  )
)

bg %<>% group_by(bin) %>% group_modify(
  ~ {
    p <- case_when(
      .y$bin == "<=1" ~ 0.05,
      .y$bin == "1~2" ~ 0.10,
      .y$bin == "2~3" ~ 0.15,
      .y$bin == "3~4" ~ 0.30,
      .y$bin == ">=4" ~ 0.80
    )
    
    slice_sample(.x, prop = p)
  }
) %>% ungroup

down <- bind_rows(bg, fg); bof <- -log10(0.05 / nrow(ext))

# Plot manhattan
manh <- ggplot(
) + geom_hline(
  yintercept  = bof,
  linetype    = "dashed",
  linewidth   = 0.35,
  show.legend = FALSE
) + geom_point(
  data = fg,
  aes(x = c_pos, y = log_p),
  size  = 0.25,
  color = "#db6d5f",
  shape = 16
) + geom_point(
  data = bg,
  aes(x = c_pos, y = log_p, color = factor(num %% 2)),
  size  = 0.25,
  shape = 16
) + geom_text(
  data = gn,
  aes(x = c_pos, y = log_p, label = gene),
  fontface = "italic",
  hjust = -0.05,
  vjust = -0.2,
  size  = 2.15,
  check_overlap = TRUE,
  inherit.aes   = FALSE
) + geom_segment(
  data = pos,
  aes(x = c_pos, xend = c_pos, y = seg_start, yend = seg_end),
  linewidth   = 0.35,
  inherit.aes = FALSE
) + geom_text(
  data = pos,
  aes(x = x_txt, y = y_txt, label = label),
  size  = 2.15,
  color = "#db6d5f",
  hjust = 0,
  vjust = 0,
  inherit.aes = FALSE
) + geom_segment(
  data = pos,
  aes(x = c_pos, xend = c_pos, y = seg_start, yend = seg_end),
  arrow = arrow(length = unit(0.01, "npc"),),
  color = "#db6d5f",
  linewidth = 0.35,
  inherit.aes = FALSE
) + scale_x_continuous(
  expand = expansion(add = 0),
  breaks = axis$center,
  labels = ifelse(
    axis$chr %in% c("chr19", "chr21"),
    paste0("\n", gsub("^chr", "", axis$chr)),
    gsub("^chr", "", axis$chr)
  )
) + scale_y_continuous(
  trans  = "sqrt",
  breaks = c(2, 5, 10, 20),
  limits = c(1, 40),
  expand = expansion(mult = c(0, 0.05))
) + scale_color_manual(
  values = c("0" = "#dadada", "1" = "#839ba3"),
  guide = "none"
) + labs(
  x = "Chromosome",
  y = expression(-log[10](italic(p)*"-value"))
)

manh <- manh + basic + labs(tag = "a") + theme(
  axis.title.x = element_text(size = 7, margin = margin(t = -5)),
  plot.margin = margin(l = 11, r = 5, b = 2, t = 5)
)

# Save
ggsave("result/manhattan.png", manh, width = 7.2, height = 1.5, dpi = 2000, bg = "transparent")

# 5. Chromatin Accessibility ───────────────────────────────────────────────────
chr <- e_forge %>% filter(
  Tissue %in% "Blood"
) %>% mutate(
  sort = case_when(
    grepl("E029", Cell) ~ "Monocyte",
    grepl("E032", Cell) ~ "B cell",
    grepl("E033", Cell) ~ "T cell (cord)",
    grepl("E034", Cell) ~ "T cell (peripheral)",
    grepl("E046", Cell) ~ "NK cell (peripheral)",
    grepl("E050", Cell) ~ "Hematopoietic (E050)",
    grepl("E051", Cell) ~ "Hematopoietic (E051)"
  ),
  star = map_vec(Qvalue, p_to_star)
)

chr$log_p <- -log10(chr$Pvalue)
chr$type  <- factor(chr$Datatype, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3", "H3K36me3"))

chr$sort %<>% reorder(chr$log_p)
chr$type %<>% fct_relevel("H3K36me3", "H3K9me3", "H3K27me3", "H3K4me3", "H3K4me1")

# Draw bubble heatmap
pal <- c("#374e55", "#8c9Ca1", "#f2eCe7", "#d27c68", "#b32e1e")

bubble <- ggplot(
  chr, aes(x = sort, y = type)
) + geom_tile(
  fill = NA, color = NA
) + geom_point(
  aes(size = Zscore, color = log_p),
  shape = 16
) + scale_color_gradientn(
  colors = pal,
  name   = expression(-log[10](italic(p)*"-value")),
  breaks = c(0, 5, 10, 15)
) + scale_size_continuous(
  range  = c(1.5, 3.75),
  breaks = c(-2.58, 0, 2.58)
) + guides(
  size  = "none",
  color = guide_colorbar(
    title.position = "right",
    direction      = "vertical",
    barwidth       = unit(1.5, "mm"),
    barheight      = unit(15, "mm"),
    ticks          = FALSE
  )
) + labs(
  x = NULL, y = NULL
) + basic + theme(
  plot.tag.position = c(0, 0.95),
  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.4
  ),
  legend.position      = "right",
  legend.direction     = "vertical",
  legend.box           = "vertical",
  legend.box.just      = "left",
  legend.justification = "top",
  legend.title         = element_text(
    angle = 270, margin = margin(l = 2)
  ),
  legend.text          = element_text(
    margin = margin(l = 2)
  ),
  legend.box.spacing = unit(1, "mm")
) + labs(
  tag = "b"
)

# Save
ggsave("result/bubble.svg", device = svglite, bubble, width = 2.2, height = 1.75, dpi = 800, bg = "transparent")

# 6. CpG Modulation ────────────────────────────────────────────────────────────
cpgs <- ext %>% slice_min(adj.P.Val, n = 1e3) %>% pull(cpg)
sub  <- ncc %>% filter(sex != "F")

x <- expr[cpgs, sub$sentrix_id] %>% impute %>% t

# Extract eigengenes
allowWGCNAThreads()
pickSoftThreshold(x, powerVector = 1:20, networkType = "signed", corFnc = "bicor", verbose = 5)

pwr <- 20  # pick very first power that exceeds 0.8 SFT.R.sq, or the best

# Define modules
net <- blockwiseModules(
  x,
  power             = pwr,
  networkType       = "signed",
  corType           = "bicor",
  maxPOutliers      = 0.05,
  minModuleSize     = 10,
  reassignThreshold = 0,
  mergeCutHeight    = 0.4,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 5
)

# Save
write_rds(net, "object/net-modules.rds")
net <- read_rds("object/net-modules.rds")

# Extract WGCNA dendrogram
soft_power <- pwr

hc  <- net$dendrograms[[1]]
den <- as.dendrogram(hc)

dd  <- dendro_data(den, type = "rectangle")
seg_all <- dd$segments

# Determine module assignment
block_idx <- net$blockGenes[[1]]
leaf_idx  <- block_idx[hc$order]

cpg_order  <- colnames(x)[leaf_idx]
module_num <- net$colors[leaf_idx]

# Module order
module_levels <- paste0("ME", sort(unique(net$colors)))

module_cols <- c(
  "ME0" = "#D7D7D7",
  "ME1" = "#71c9dd",
  "ME2" = "#8480c2",
  "ME3" = "#6F758A",
  "ME4" = "#82617E",
  "ME5" = "#954D73",
  "ME6" = "#B44770"
)

module_cols <- module_cols[module_levels]

# ── terminal edge의 시작점만 추출 ──────────────────────────────
terminal_end <-  seg_all %>%
  filter(
    abs(x - xend) < 1e-12,
    abs(yend) < 1e-12
  ) %>%
  transmute(
    x = x,
    y = y,
    leaf_pos = as.integer(round(x))
  ) %>%
  filter(
    leaf_pos >= 1,
    leaf_pos <= length(module_num)
  ) %>%
  mutate(
    module = factor(
      paste0("ME", module_num[leaf_pos]),
      levels = module_levels
    )
  ) %>%
  filter(module != "ME0")

# 원래대로 긴 terminal vertical edge는 제거
seg <- seg_all %>%
  filter(
    !(abs(x - xend) < 1e-12 & abs(yend) < 1e-12)
  )

# Absolute biweight mid-correlation
abs_bicor <- function(z, y) {
  ok <- is.finite(z) & is.finite(y)
  
  if (sum(ok) < 4) return(NA_real_)
  if (sd(z[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  
  abs(as.numeric(
    WGCNA::bicor(
      z[ok],
      y[ok],
      maxPOutliers = 0.05
    )
  ))
}

# CpG correlation with risk score
risk_cor <- apply(
  x,
  2,
  abs_bicor,
  y = sub$risk_score
)

# CpG correlation with Time-to-diagnosis
case_idx <- sub$group == "prediagnostic" & !is.na(sub$time)

time_cor <- apply(
  x[case_idx, , drop = FALSE],
  2,
  abs_bicor,
  y = sub$time[case_idx]
)

# Module track
module_track <- tibble(
  x      = seq_along(leaf_idx),
  track  = "Module",
  module = factor(paste0("ME", module_num), levels = module_levels)
)

# Correlation tracks
cor_track <- bind_rows(
  tibble(
    x     = seq_along(leaf_idx),
    track = "Risk score",
    cor   = risk_cor[leaf_idx]
  ),
  tibble(
    x     = seq_along(leaf_idx),
    track = "Time-to-diagnosis",
    cor   = time_cor[leaf_idx]
  )
)

# Common correlation scale
cor_max <- max(cor_track$cor, na.rm = TRUE)

if (!is.finite(cor_max) || cor_max == 0) {
  cor_max <- 1
}

# Horizontal / vertical branches
seg <- seg %>% mutate(horizontal = abs(y - yend) < 1e-12)

# Draw dendrogram
tree <- ggplot() +
  geom_segment(
    data = seg %>% filter(!horizontal),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.1,
    color = "black",
    alpha = 0.5,
    lineend = "butt"
  ) +
  geom_segment(
    data = seg %>% filter(horizontal),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.20,
    color = "black",
    alpha = 0.85,
    lineend = "butt"
  ) +
  # terminal line-end를 ME 색으로 표시
  geom_point(
    data = terminal_end,
    aes(
      x = x,
      y = y,
      color = module
    ),
    shape = 16,
    size = 0.2,
    alpha = 0.5,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = module_cols,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0.5, length(module_num) + 0.5),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    trans = scales::trans_new(
      "power3",
      transform = function(x) x^3,
      inverse   = function(x) x^(1/3)
    ),
    breaks = c(0.6, 0.8, 1.0),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = NULL,
    y = "TOM dissimilarity"
  )

xr <- range(c(seg$x, seg$xend))
yr <- range(c(seg$y, seg$yend))

tree <- tree + annotate(
  "text",
  x = xr[2],
  y = yr[1],
  label = paste0("Biweight mid-correlation\nsoft-threshold \u03B2 = ", soft_power),
  # label = paste0(
  #   "'Soft-threshold'~italic(beta)==", soft_power
  # ),
  # parse = TRUE,
  hjust = 1,
  vjust = 0,
  size = 2.0,
  lineheight = 0.9
)

tree <- tree + basic + theme(
  # plot.margin  = margin(l = 10, r = 1, b = 2, t = 2),
  # plot.tag     = element_text(
  #   family = "Arial",
  #   face = "bold",
  #   size = 16,
  #   margin = margin(l = -18, t = -4)
  # ),
  plot.tag.position = c(0, 1),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.x  = element_blank(),
  axis.title.x = element_blank(),
  # axis.text.y  = element_text(size = 7),
  # axis.title.y = element_text(size = 8),
  axis.ticks.y = element_line(linewidth = 0.35),
  axis.line.y  = element_line(linewidth = 0.35)
) + labs(
  tag = "c"
)

# Multi-row annotation strip
strip <- ggplot() +
  geom_raster(
    data = module_track,
    aes(x = x, y = track, fill = module),
    interpolate = FALSE
  ) +
  scale_fill_manual(
    values = module_cols,
    breaks = setdiff(module_levels, "ME0"),
    drop = FALSE,
    name = NULL,
    guide = guide_legend(
      # nrow = 1,
      ncol = 1,
      byrow = TRUE,
      keywidth = unit(1.8, "mm"),
      keyheight = unit(1.8, "mm"),
      override.aes = list(colour = NA)
    )
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = cor_track,
    aes(x = x, y = track, fill = cor),
    interpolate = FALSE
  ) +
  scale_fill_gradient(
    low = "#ffffff",
    high = "#d1321f",
    limits = c(0, cor_max),
    oob = scales::squish,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0.5, length(module_num) + 0.5),
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(
    limits = c(
      "Time-to-diagnosis",
      "Risk score",
      "Module"
    ),
    labels = c(
      "Module" = "ME",
      "Risk score" = "RS",
      "Time-to-diagnosis" = "TTD"
    ),
    position = "left",
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(
    # plot.margin = margin(l = 10, r = 1, b = 15, t = 1),
    axis.text.y = element_text(
      size = 6,
      color = "black",
      hjust = 1,
      margin = margin(r = 2)
    ),
    legend.position = "bottom",
    # legend.position.inside = c(0.5, -0.75),
    # legend.justification = c(0.5, 0),
    legend.direction = "horizontal",
    
    legend.title = element_blank(),
    
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.text = element_text(
      size = 6,
      margin = margin(l = 2)
    ),
    legend.key.spacing.x = unit(2, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  )

# Combine
network <- (tree / strip) +
  plot_layout(
    heights = c(0.6, 0.15),
    guides = "collect"
  ) + plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = NA,
        colour = NA
      )
    )
  ) &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "top",
    legend.box.just = "top",
    legend.title = element_blank(),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.key.spacing.y = unit(0.75, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(
      size = 6,
      margin = margin(l = 2)
    ),
    legend.box.spacing = unit(3, "mm"),
    legend.spacing.y = unit(2, "mm")
  )

# Save
ggsave("result/wgc-net.svg", device = svglite, network, width = 3, height = 1.6, dpi = 800, bg = "transparent")

# ── Prepare CHIP overlap / enrichment ─────────────────────────────────────────
mod_map <- tibble(cpg = cpgs, module = paste0("ME", unname(net$colors)))

chip_info <- esm %>%
  select(group = chip_group, base = CpG)

agg <- ext %>%
  transmute(
    cpg,
    base   = str_remove(cpg, "_.*$"),
    log_fc = logFC,
    t,
    p_value = P.Value,
    p_adj   = adj.P.Val
  ) %>%
  left_join(mod_map, by = "cpg") %>%
  left_join(chip_info, by = "base", relationship = "many-to-many") %>%
  drop_na(module)

# Background and unique CpG × module pairs (same approach as original script)
bg <- anno %>%
  filter(!is.na(prev), prev != "") %>%
  pull(base) %>%
  unique()

u0 <- agg %>% distinct(base, module)
fg <- intersect(unique(u0$base), bg)
modules <- unique(u0$module)

chip_info <- chip_info %>%
  filter(base %in% bg) %>%
  distinct()

chip_sets <- split(chip_info$base, chip_info$group)

ovr <- imap_dfr(chip_sets, function(vec, nm) {
  map_dfr(modules, function(m) {
    module_cpgs <- u0 %>% filter(module == m) %>% pull(base)
    query <- intersect(module_cpgs, fg)
    fisher(query, bg, vec) %>%
      mutate(regulator = nm, module = m, .before = 1)
  })
}) %>%
  mutate(p_adj = p.adjust(p, method = "fdr")) %>%
  filter(regulator != "any") %>%
  mutate(
    n_bg  = a + b + c + d,
    n_mod = a + b,
    p_exp = (a + c) / n_bg,
    p_obs = a / n_mod
  )

# Make x-position and facet order deterministic
chip_genes <- c("ASXL1", "DNMT3A", "TET2")
ovr <- ovr %>%
  mutate(
    regulator = factor(regulator, levels = chip_genes),
    module    = factor(module, levels = module_levels)
  ) %>%
  filter(!is.na(regulator), !is.na(module))

# ── Bar chart data ─────────────────────────────────────────────────────────────
bar <- ovr %>%
  transmute(module, regulator,
            x = as.numeric(regulator), prop = p_obs,
            key = as.character(regulator))

exp <- ovr %>%
  transmute(module, regulator,
            x = as.numeric(regulator), prop = p_exp)

tr <- pseudo_log_trans(sigma = 0.05)

sig <- ovr %>%
  filter(!is.na(p_adj), p_adj < 0.05, b > 2) %>%
  mutate(
    x    = as.numeric(regulator),
    star = sapply(p_adj, p_to_star),
    base = tr$transform(pmax(p_exp, p_obs)),
    y    = tr$inverse(base + 0.12)
  )

pal <- c(
  "ASXL1"  = "#7cb8ba",
  "DNMT3A" = "#5b7e8a",
  "TET2"   = "#db6d5f"
)

# ── Compare EWAS direction with CHIP-study beta ──────────────────────────────
# Compare at unique CpG base level. A base with inconsistent signs among its
# probes is excluded from the direction denominator, not from Fisher counts.
unique_direction <- function(v) {
  signs <- unique(sign(v[is.finite(v) & v != 0]))
  if (length(signs) == 1L) signs[[1L]] else NA_real_
}

chip_direction <- esm %>%
  filter(chip_group %in% chip_genes) %>%
  transmute(
    regulator = as.character(chip_group),
    base      = as.character(CpG),
    effect    = as.numeric(beta)
  ) %>%
  filter(base %in% bg) %>%
  group_by(regulator, base) %>%
  summarise(chip_sign = unique_direction(effect), .groups = "drop")

ewas_direction <- agg %>%
  filter(base %in% fg) %>%
  distinct(base, module, cpg, log_fc) %>%
  group_by(base, module) %>%
  summarise(ewas_sign = unique_direction(log_fc), .groups = "drop")

chip_direction_summary <- ewas_direction %>%
  inner_join(chip_direction, by = "base", relationship = "many-to-many") %>%
  mutate(
    direction = case_when(
      is.na(ewas_sign) | is.na(chip_sign) ~ NA_character_,
      ewas_sign == chip_sign           ~ "Same",
      TRUE                            ~ "Opposite"
    )
  ) %>%
  group_by(module, regulator) %>%
  summarise(
    n_comparable = sum(!is.na(direction)),
    n_same       = sum(direction == "Same", na.rm = TRUE),
    n_opposite   = sum(direction == "Opposite", na.rm = TRUE),
    .groups = "drop"
  )

# n: number of CHIP-overlapping CpGs (ovr$a)
# Total: all eligible CpGs within the module (ovr$n_mod)
# S/O (%): concordant/opposite directions among comparable overlap CpGs
chip_table_wide <- ovr %>%
  transmute(module,
            regulator = as.character(regulator),
            n = a, Total = n_mod) %>%
  left_join(chip_direction_summary,
            by = c("module", "regulator")) %>%
  mutate(
    across(c(n_comparable, n_same, n_opposite), ~ replace_na(.x, 0L)),
    same_pct = if_else(
      n_comparable > 0,
      round(100 * n_same / pmax(n_comparable, 1)),
      NA_real_
    ),
    so_label = if_else(
      n_comparable > 0,
      sprintf("%d/%d", as.integer(same_pct),
              as.integer(100 - same_pct)),
      "–"
    )
  )

stopifnot(all(chip_table_wide$n_comparable <= chip_table_wide$n))

# One total value (one merged visual column) per module
module_totals <- chip_table_wide %>%
  group_by(module) %>%
  summarise(
    Total = first(Total),
    n_distinct_total = n_distinct(Total),
    .groups = "drop"
  )
stopifnot(all(module_totals$n_distinct_total == 1L))

# Unboxed table: only n and S/O values are repeated under each gene.
table_numbers <- chip_table_wide %>%
  mutate(
    x = match(regulator, chip_genes),
    module = factor(module, levels = module_levels)
  )

table_numbers <- table_numbers %>%
  mutate(
    so_label = if_else(so_label == "100/0", "All", so_label)
  )

module_totals <- module_totals %>%
  mutate(module = factor(module, levels = module_levels))

# ── Original bar plot (no facet headers) ──────────────────────────────────────
enr <- ggplot(bar, aes(x = x, y = prop, fill = key)) +
  geom_col(width = 0.45) +
  geom_segment(
    data = exp,
    aes(x = x - 0.27, xend = x + 0.27,
        y = prop, yend = prop, color = "Expected"),
    linewidth = 0.55,
    lineend = "butt",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = sig,
    aes(x = x, y = y, label = star),
    size = 2,
    vjust = 0,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ module, nrow = 1) +
  scale_y_continuous(
    trans  = tr,
    limits = c(0, 2),
    breaks = c(0, 0.25, 0.5, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    limits = c(0.5, 3.5),
    breaks = 1:3,
    labels = rep("", 3),
    expand = expansion(add = 0)
  ) +
  scale_fill_manual(
    values = alpha(pal, 0.9),
    name   = NULL,
    breaks = chip_genes,
    labels = c(
      expression(italic("ASXL1")),
      expression(italic("DNMT3A")),
      expression(italic("TET2"))
    ),
    guide = guide_legend(order = 2)
  ) +
  scale_color_manual(
    values = c("Expected" = "#5f5f5f"),
    name = NULL,
    guide = guide_legend(
      order = 1,
      override.aes = list(linewidth = 0.45)
    )
  ) +
  # coord_cartesian(ylim = c(0, 1.2), clip = "off") +
  labs(x = NULL, y = "CHIP-associated", tag = "d") +
  basic +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 6,
      lineheight = 0.85,
      margin = margin(t = 1, b = 1)
    ),
    panel.spacing.x = unit(1, "mm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6, margin = margin(l = 2)),
    legend.key.height = unit(2.3, "mm"),
    legend.key.width = unit(2.3, "mm"),
    legend.key.spacing.x = unit(2, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.margin = margin(0, 0, 0, b = -5)
  )

# ── Text-only table: no tile fills, no cell borders ───────────────────────────
chip_table_plot <- ggplot() +
  # Hairline horizontal row dividers, not individual cell outlines
  geom_hline(
    yintercept = c(2.5, 1.5, 0.5),
    linewidth = 0.18,
    colour = "black"
  ) +
  # Row 1: n under each CHIP gene
  geom_text(
    data = table_numbers,
    aes(x = x, y = 3, label = n),
    size = 1.8,
    family = "Arial"
  ) +
  # Row 2: one centered Total per ME
  geom_text(
    data = module_totals,
    aes(x = 2, y = 2, label = Total),
    size = 1.8,
    family = "Arial"
  ) +
  # Row 3: same/opposite directional percentage under each gene
  geom_text(
    data = table_numbers,
    aes(x = x, y = 1, label = so_label),
    size = 1.8,
    family = "Arial"
  ) +
  facet_wrap(~ module, nrow = 1) +
  scale_x_continuous(
    limits = c(0.5, 3.5),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, 3.45),
    breaks = c(1, 2, 3),
    labels = expression(
      "S/O (%)",
      "Total",
      "Overlap"
    ),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(x = NULL, y = NULL) +
  basic +
  theme(
    panel.spacing.x = unit(1, "mm"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(
      size = 6,
      hjust = 1,
      margin = margin(r = 2)
    ),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

# ── Combine and save ──────────────────────────────────────────────────────────
chip_enrichment_table <- (enr / chip_table_plot) +
  plot_layout(heights = c(1.0, 0.5), guides = "collect") + plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = NA,
        colour = NA
      )
    )
  ) &
  theme(
    legend.position = "top",
    plot.margin = margin(l = 5, t = 5, b = 1, r = 5)
  )

ggsave(
  "result/chip-enrichment.svg",
  plot = chip_enrichment_table,
  device = svglite::svglite,
  width = 5.15,
  height = 1.5,
  bg = "transparent"
)

# 7. Trait Exploration ─────────────────────────────────────────────────────────

egs <- moduleEigengenes(x, colors = net$colors)$eigengenes
egs %<>% orderMEs

# Compute correlation between traits and module eigengenes
tr <- imap_dfr(
  egs, \(eg, mod) {
    names(eg) <- rownames(egs)
    
    eg <- eg[sub$sentrix_id]
    ev <- as.numeric(sub$event)
    
    # Helper to organize correlation
    org_cor <- function(y, var, extra = NULL) {
      cond <- is.null(extra)
      if (cond) extra <- rep(TRUE, length(y))
      
      idx <- complete.cases(eg, y) & extra
      
      n <- sum(idx)
      r <- cor(eg[idx], y[idx])
      p <- cor_p(r, n)
      
      org <- data.frame(
        module = mod,
        variable = var,
        correlation = r,
        p_value = p,
        n = n
      )
      
      return(org)
    }
    
    # Compute correlations
    out <- sub %>% select(
      -sentrix_id
    ) %>% imap_dfr(
      \(x, nm) {
        is_num <- is.numeric(x) || is.integer(x)
        
        # Time-to-diagnosis in cases
        if (nm %in% "time") {
          res <- org_cor(
            x, "time",
            !is.na(ev) & ev != 0
          )
          
          return(res)
        }
        
        # Smoking traits
        if (nm == "smk_status") {
          targets <- intersect(
            c("quit", "current"), levels(x)
          )
          
          res <- map_dfr(targets, \(d) {
            y <- case_when(
              is.na(x)    ~ NA,
              x == "none" ~ 0,
              x == d      ~ 1,
              TRUE        ~ NA
            )
            
            t <- paste0("smk_status=", d)
            org_cor(y, t)
          })
          
          return(res)
        }
        
        # Each morphology to controls
        if (nm %in% "morphology") {
          res <- map_dfr(
            levels(x), \(d) {
              y <- if_else(
                is.na(ev),
                NA,
                if_else(
                  ev == 0, 0,
                  if_else(x == d, 1, NA)
                )
              )
              
              t <- paste0("morphology=", d)
              org_cor(y, t)
            }
          )
          
          return(res)
        }
        
        # Other numeric traits
        if (is_num) {
          if (nm == "smk_dur") {
            res <- org_cor(
              x, nm,
              !is.na(sub$smk_status) &
                sub$smk_status != "none"
            )
            return(res)
          }
          
          if (nm == "ces_dur") {
            res <- org_cor(
              x, nm,
              !is.na(sub$smk_status) &
                sub$smk_status == "quit"
            )
            return(res)
          }
          
          res <- org_cor(x, nm)
          return(res)
        }
        
        # Categorical traits
        res <- map_dfr(
          levels(x), \(d) {
            y <- if_else(
              is.na(x),
              NA,
              as.numeric(x == d)
            )
            
            t <- paste0(nm, "=", d)
            org_cor(y, t)
          }
        )
        
        return(res)
      }
    )
    
    out
  }
)

morph_labels <- c(
  "Adenocarcinoma (incl. variants)" = "LUAD (incl. variants)",
  "Squamous cell carcinoma" = "LUSC (squamous)",
  "Small cell carcinoma" = "SCLC (small cell)",
  "Non-small cell lung cancer" = "NSCLC",
  "Large cell carcinoma" = "LCLC (large cell)",
  "Neuroendocrine carcinoma" = "NEC (neuroendocrine)",
  "NOS, unspecified" = "NOS (unspecified)"
)

morph_short <- function(x) {
  dplyr::recode(x, !!!morph_labels)
}

# ── Organize traits ───────────────────────────────────────────────────────────

sup <- c("Mye", "Lym")
ty  <- sub %>% select(Bas:Lym) %>% names %>% setdiff(sup)

tr %<>% drop_na %>% mutate(
  base = sub("=.*$", "", variable),
  level = if_else(
    grepl("=", variable),
    sub("^.*=", "", variable),
    NA
  ),
  group = case_when(
    grepl("^smk|ces|cot", base) ~ "smoke",
    grepl("^morph", base)       ~ "morph",
    base %in% sup               ~ "super",
    base %in% ty                ~ "type",
    TRUE                        ~ "any"
  ),
  p_adjust = p.adjust(p_value, method = "bonferroni"),
  variable = case_when(
    base == "event"      ~ "Lung cancer",
    base == "time"       ~ "Time-to-diagnosis",
    base == "risk_score" ~ "Risk score",
    base == "morphology" ~ morph_short(level),
    base == "age"        ~ "Age",
    base == "sex=M"      ~ "Sex",
    base == "smk_dur"    ~ "Duration: smoke",
    base == "ces_dur"    ~ "Duration: quit",
    base == "ces_age"    ~ "Quit age",
    base == "smk_status" & level == "none"    ~ "Smoking: none",
    base == "smk_status" & level == "quit"    ~ "Smoking: quit",
    base == "smk_status" & level == "current" ~ "Smoking: current",
    base == "cotinine"   ~ "Cotinine",
    base == "Lym"        ~ "Lymphoid",
    base == "Mye"        ~ "Myeloid",
    base == "CD4tot"     ~ "CD4 T (total)",
    base == "CD8tot"     ~ "CD8 T (total)",
    base == "Btot"       ~ "B (total)",
    base == "CD4nv"      ~ "CD4 T (naïve)",
    base == "CD8nv"      ~ "CD8 T (naïve)",
    base == "Bnv"        ~ "B (naïve)",
    base == "CD4mem"     ~ "CD4 T (memory)",
    base == "CD8mem"     ~ "CD8 T (memory)",
    base == "Bmem"       ~ "B (memory)",
    base == "Treg"       ~ "T reg",
    base == "Neu"        ~ "Neutrophil",
    base == "Mono"       ~ "Monocyte",
    base == "Eos"        ~ "Eosinophil",
    base == "Bas"        ~ "Basophil",
    base == "NK"         ~ "NK cell",
    TRUE                ~ base
  )
)

# ── Variable order ────────────────────────────────────────────────────────────

frm <- c("any", "smoke", "super", "type", "morph")
ord <- c(
  "Lung cancer", "Time-to-diagnosis", "Risk score",
  "Age", "Sex",
  morph_short(levels(sub$morphology)),
  "Smoking: none", "Smoking: quit",
  "Smoking: current", "Cotinine",
  "Lymphoid", "Myeloid",
  "Neutrophil", "Monocyte", "Eosinophil",
  "Basophil", "NK cell", "T reg",
  "CD4 T (total)", "CD8 T (total)", "B (total)",
  "CD4 T (naïve)", "CD8 T (naïve)", "B (naïve)",
  "CD4 T (memory)", "CD8 T (memory)", "B (memory)"
)

tr %<>% filter(variable %in% ord)

# ── Organize heatmap ──────────────────────────────────────────────────────────

# CHANGED 1: Include morphology in heatmap
ax <- tr %>% ungroup()

ax$group <- factor(ax$group, levels = frm)
ax$variable <- factor(ax$variable, levels = rev(ord))
ax$star <- sapply(ax$p_adjust, p_to_star)

# ── Clustering ────────────────────────────────────────────────────────────────

# Preserve original clustering by excluding morphology
mt <- ax %>%
  filter(base != "morphology") %>%
  group_by(module, variable) %>%
  summarise(
    correlation = mean(correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = correlation
  ) %>%
  column_to_rownames("module")

hc <- hclust(dist(mt), method = "average")
lv <- hc$labels[hc$order]
n_mod <- length(lv)

ax$module <- factor(ax$module, levels = lv)
ax$star <- sub("n.s", "", ax$star)

dd <- dendro_data(hc, type = "rectangle")

# ── Dendrogram ────────────────────────────────────────────────────────────────

dend <- ggplot(dd$segments) +
  geom_segment(
    aes(
      x = x, y = y,
      xend = xend, yend = yend
    ),
    linewidth = 0.25
  ) +
  scale_x_continuous(
    limits = c(0.5, n_mod + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  basic +
  theme(
    plot.margin = margin(l = 5, r = 5, b = 0, t = 5),
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0.05, 1.02),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(tag = "g")

# ── Heatmap ───────────────────────────────────────────────────────────────────

heat <- ggplot(
  ax,
  aes(
    x = module,
    y = variable,
    fill = correlation
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = star,
      vjust = if_else(star == "·", 0.7, 0.8)
    ),
    size = 2,
    na.rm = TRUE
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_discrete(
    expand = expansion(add = 0)
  ) +
  
  scale_y_discrete(
    position = "left",
    expand = expansion(add = 0)
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "#478a8f",
      "#88bec2",
      "#ccebed",
      "#F8FAFA",
      "#FFFFFF",
      "#FFF4EC",
      "#F0A06A",
      "#D66A32"
    ),
    values = scales::rescale(
      c(
        -1, -0.60, -0.30, -0.18,
        0,  0.18,  0.30,  0.60, 1
      )
    ),
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      direction = "horizontal",
      barwidth = unit(17, "mm"),
      barheight = unit(1.5, "mm"),
      title.hjust = 0.5
    )
  ) +
  basic +
  theme(
    
    plot.margin = margin(l = 5, r = 5, b = 5, t = 0),
    strip.text = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.title = element_text(
      angle = 0,
      size = 6,
      margin = margin(b = 3, t = -6)
    ),
    legend.text = element_text(
      size = 5.8,
      margin = margin(t = 2, b = -4, l = 2)
    )
  )

# ── Combine ───────────────────────────────────────────────────────────────────

trait <- (dend / heat) +
  plot_layout(
    heights = c(0.3, 6)
  ) + plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = NA,
        colour = NA
      )
    )
  )

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(
  "result/trait-map.svg",
  plot = trait,
  device = svglite::svglite,
  width = 2,
  height = 5,
  dpi = 800,
  bg = "transparent"
)

