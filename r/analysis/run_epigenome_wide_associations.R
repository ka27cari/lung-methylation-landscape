# 1. Set Environment ───────────────────────────────────────────────────────────
library(limma)
library(purrr)
library(parallel)
library(data.table)
library(missMethyl)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Set priority
conflicts_prefer(base::intersect, matrixStats::rowMedians, matrixStats::rowVars, purrr::reduce, dplyr::count)

# Load data
m <- file.path(data_dir, "m-adjusted.rds") %>% read_rds

coh  <- read_rds("object/obs-pred.rds")
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# CHIP information
esm <- file.path(manifest_dir, "41467_2025_59333_MOESM34_ESM.xlsx") %>% read_xlsx

# 2. Generate Functions ────────────────────────────────────────────────────────
# Build paired data
make_pair <- function(data, origin = NULL, time_min = NULL, time_max = NULL) {
  if (!is.null(origin)) data %<>% filter(origin %in% !!origin)
  
  prs <- data %>% group_by(
    origin, matched_index
  ) %>% summarise(
    matched_index = matched_index[group == "prediagnostic"][1],
    time = time[group == "prediagnostic"][1],
    type = type[group == "prediagnostic"][1],
    hat  = .hat[group == "prediagnostic"][1],
    .groups = "drop"
  )
  
  if (!is.null(time_max)) prs %<>% filter(time < time_max)
  if (!is.null(time_min)) prs %<>% filter(time >= time_min)
  
  prs %<>% mutate(
    axis = case_when(
      type %in% "LUAD" ~ "Ad",
      type %in% c("LUSC", "SCLC") ~ "Sc",
      TRUE ~ "Other"
    )
  )
  
  prs$axis %<>% fct_relevel("Sc", "Ad", "Other")
  
  return(prs)
}

# Contrast generator
make_contrasts <- function(
    nms, axs = c("Sc", "Ad")
) {
  n <- length(nms)
  z <- setNames(rep(0, n), nms)
  
  nm <- function(ax, grp) paste(ax, grp, sep = "_")
  
  basis <- function(ax, grp, coef = 1) {
    v <- z
    key <- nm(ax, grp)
    
    if (key %in% names(v)) v[key] <- coef
    v
  }
  
  delta <- function(ax) {
    basis(ax, "prediagnostic", 1) + basis(ax, "control", -1)
  }
  
  gap <- function(grp) {
    basis("Ad", grp, 1) + basis("Sc", grp, -1)
  }
  
  contr <- list()
  
  contr[["dAd"]] <- delta("Ad")
  contr[["dSc"]] <- delta("Sc")
  
  contr[["dPool"]] <- (contr[["dAd"]] + contr[["dSc"]]) / 2
  contr[["dAx"]]   <- contr[["dAd"]] - contr[["dSc"]]
  
  contr[["cCase"]] <- gap("prediagnostic")
  contr[["cCtr"]]  <- gap("control")
  
  do.call(cbind, contr)
}

# Pooled run_axis
run_axis <- function(
    meta,
    origin   = NULL,
    time_min = NULL,
    time_max = NULL,
    alpha    = 1,
    covars   = c("ME1", "ME5", "ME6"),
    equal    = TRUE,
    block    = FALSE
) {
  if (!is.null(origin)) meta %<>% filter(origin %in% !!origin)
  
  prs <- make_pair(
    meta,
    origin   = origin,
    time_min = time_min,
    time_max = time_max
  )
  
  idx <- prs$matched_index
  
  meta %<>% filter(
    matched_index %in% idx
  ) %>% left_join(
    prs, by = c("origin", "matched_index")
  ) %>% filter(
    axis != "Other", sentrix_id %in% colnames(m)
  ) %>% droplevels
  
  meta$g <- interaction(meta$axis, meta$group, sep = "_") %>% droplevels
  expr <- m[, meta$sentrix_id, drop = FALSE]
  
  rhs <- c("0 + g")
  if (length(covars) > 0) rhs <- c(rhs, covars)
  
  term <- paste("~", paste(rhs, collapse = " + ")) %>% as.formula
  
  design <- model.matrix(term, data = meta)
  colnames(design) %<>% sub("^g", "", .)
  
  w_spl <- rep(1, nrow(meta))
  
  if (equal) {
    cl  <- interaction(meta$axis, drop = TRUE)
    tab <- table(cl)
    
    w_spl <- 1 / (tab[cl] ^ alpha)
    w_spl <- w_spl / mean(w_spl)
  }
  
  contr <- make_contrasts(
    nms = colnames(design)
  )
  
  if (block) {
    b   <- interaction(meta$matched_index, drop = TRUE)
    obs <- duplicateCorrelation(expr, design = design, block = b, weights = w_spl)
    fit <- lmFit(expr, design = design, block = b, correlation = obs$consensus, weights = w_spl)
  } else {
    fit <- lmFit(expr, design, weights = w_spl)
  }
  
  eb <- eBayes(contrasts.fit(fit, contr), robust = TRUE, trend = TRUE)
  
  coef  <- eb$coefficients
  p_val <- eb$p.value
  
  top <- data.frame(cpg = rownames(coef), stringsAsFactors = FALSE)
  
  for (nm in colnames(coef)) {
    top[[paste0(nm, "_eff")]] <- coef[, nm]
    top[[paste0(nm, "_p")]]   <- p_val[, nm]
  }
  
  top %<>% mutate(
    origin = origin,
    window = if_else(
      time_min != 0, "ds", sprintf("imm_%sy", time_max)
    ),
    alpha  = alpha,
    .after = cpg
  )
  
  top
}

# Compute signed z-statistic
signed_z <- function(b, p) {
  p <- pmin(pmax(p, 1e-300), 1 - 1e-16)
  sign(b) * qnorm(1 - p / 2)
}

# Merge z-statistic
stouffer <- function(b1, p1, b2, p2, w1, w2) {
  z1 <- signed_z(b1, p1)
  z2 <- signed_z(b2, p2)
  
  z  <- (w1 * z1 + w2 * z2) / sqrt(w1 ^ 2 + w2 ^ 2)
  p  <- 2 * pnorm(-abs(z))
  
  tibble(z = z, p = p)
}

# Directional coherence
same <- function(x, y) sign(x) == sign(y) & sign(x) != 0 & sign(y) != 0

# Transform log into score
lp <- function(x) -log10(pmax(x, 1e-300))

# Assign axis
assign <- function(data, origin, prm = params) {
  dAd_eff   <- data[[paste0("dAd_eff_", origin)]]
  dSc_eff   <- data[[paste0("dSc_eff_", origin)]]
  dPool_eff <- data[[paste0("dPool_eff_", origin)]]
  dAx_eff   <- data[[paste0("dAx_eff_", origin)]]
  
  dAd_p     <- data[[paste0("dAd_p_", origin)]]
  dSc_p     <- data[[paste0("dSc_p_", origin)]]
  dPool_p   <- data[[paste0("dPool_p_", origin)]]
  dAx_p     <- data[[paste0("dAx_p_", origin)]]
  
  score_ad <- lp(dAd_p)
  score_sc <- lp(dSc_p)
  
  sc <- dSc_p < prm$p_org &
    dAx_eff < -prm$eff &
    (dAx_p < prm$p_ax | score_sc > score_ad)
  
  ad <- dAd_p < prm$p_org &
    dAx_eff > prm$eff &
    (dAx_p < prm$p_ax | score_ad > score_sc)
  
  bth <- dPool_p < prm$p_org &
    same(dAd_eff, dSc_eff) &
    abs(dAx_eff) <= prm$eff
  
  case_when(
    sc   ~ "Sc",
    ad   ~ "Ad",
    bth  ~ "Both",
    TRUE ~ "None"
  ) %>% factor(
    levels = c("Ad", "Sc", "Both", "None")
  )
}

# Background thinning
thin_to <- function(data, n) {
  
  if (nrow(data) <= n) {
    return(data)
  }
  
  slice_sample(
    data,
    n = n
  )
}

# 3. EWAS ──────────────────────────────────────────────────────────────────────
org   <- levels(coh$origin)
times <- tibble(
  start = c(0, 0, 5),
  end   = c(3, 5, Inf)
)

alpha <- seq(0, 1, by = 0.25)

# Runs
axis <- expand_grid(
  origin = org,
  times,
  alpha = alpha
) %>% pmap_dfr(
  \(origin, start, end, alpha) {
    run_axis(
      meta     = coh,
      origin   = origin,
      alpha    = alpha,
      time_min = start,
      time_max = end,
      covars   = character(0),  # c("ME1", "ME5", "ME6")
      equal    = TRUE,
      block    = FALSE
    )
  }
)

# Save
write_csv(axis, "result/table/top-axis.csv")

# Weighted Stouffer meta-analysis
n <- length(times$start)

wgt <- expand_grid(
  origin = org,
  idx = seq_len(n)
) %>% mutate(
  start  = times$start[idx],
  end    = times$end[idx],
  window = if_else(
    start != 0, "ds", sprintf("imm_%sy", end)
  ),
  dataset = pmap(
    list(origin, start, end),
    \(origin, start, end) {
      make_pair(
        data     = coh,
        origin   = origin,
        time_min = start,
        time_max = end
      ) %>%
        filter(axis != "Other")
    }
  ),
  n_pairs = map_int(dataset, nrow),
  weight  = sqrt(n_pairs)
) %>% select(
  window, origin, weight
) %>% pivot_wider(
  names_from   = origin,
  values_from  = weight,
  names_prefix = "w_"
)

# Widen by cohort
w <- axis %>% pivot_wider(
  id_cols     = c(cpg, window, alpha),
  names_from  = origin,
  values_from = c(
    dAd_eff, dSc_eff, dPool_eff, dAx_eff, cCase_eff, cCtr_eff,
    dAd_p,   dSc_p,   dPool_p,   dAx_p,   cCase_p,   cCtr_p
  ),
  names_glue = "{.value}_{origin}"
) %>% left_join(
  wgt, by = "window"
)

# Make meta beta, z and p for each contrast
cons <- c("dAd", "dSc", "dPool", "dAx", "cCase", "cCtr")

for (c in cons) {
  b1 <- paste0(c, "_eff_KCPS")
  p1 <- paste0(c, "_p_KCPS")
  b2 <- paste0(c, "_eff_NCC")
  p2 <- paste0(c, "_p_NCC")
  
  stf <- stouffer(
    b1 = w[[b1]],
    p1 = w[[p1]],
    b2 = w[[b2]],
    p2 = w[[p2]],
    w1 = w$w_KCPS,
    w2 = w$w_NCC
  )
  
  w[[paste0(c, "_eff_meta")]] <-
    (w$w_KCPS * w[[b1]] + w$w_NCC * w[[b2]]) /
    (w$w_KCPS + w$w_NCC)
  
  w[[paste0(c, "_z_meta")]] <- stf$z
  w[[paste0(c, "_p_meta")]] <- stf$p
}

# BH-correction on meta p-values
for (c in cons) {
  p <- paste0(c, "_p_meta")
  q <- paste0(c, "_q_meta")
  
  w %<>% group_by(
    window, alpha
  ) %>% mutate(
    !!q := p.adjust(.data[[p]], method = "fdr")
  ) %>% ungroup
}

# Glance
pfx <- c("dSc", "dAd", "dPool")

f <- w %>% filter(
  if_any(
    all_of(paste0(pfx, "_q_meta")),
    ~ .x < 0.2
  )
) %>% group_by(
  window, alpha
)

f$grp <- paste(f$window, f$alpha, sep = "_")

f %>% summarise(
  across(
    all_of(paste0(pfx, "_p_NCC")), ~ {
      oth <- sub("_p_NCC$", "_p_KCPS", cur_column())
      
      sum(.x < 5e-2 & get(oth) < 5e-2, na.rm = TRUE)
    }
  ),
  .groups = "drop"
) %>% rename_with(
  ~ sub("_p_NCC$", "", .x), ends_with("_p_NCC")
)
  
f %>% ungroup %>% distinct(
  grp, cpg
) %>% inner_join(
  f, by = "cpg", relationship = "many-to-many"
) %>% count(
  grp.x, grp.y, name = "n"
) %>% pivot_wider(
  names_from = grp.y,
  values_from = n,
  values_fill = 0
) %>% arrange(
  grp.x
)

# Narrow down sets
w %<>% filter(
  window %in% c("ds", "imm_3y", "imm_5y"), alpha %in% c(0, 0.25, 0.5)
)

ds <- w %>% filter(
  !grepl("imm_", window)
) %>% select(
  cpg, alpha,
  matches("^(dAd|dSc|dPool|dAx|cCase|cCtr|hAd|hSc|hCase|hCtr)_")
) %>% rename_with(
  ~ paste0(.x, "_ds"), -c(cpg, alpha)
)

ext <- w %>% filter(
  window != "ds"
) %>% left_join(
  ds, by = c("cpg", "alpha")
) %>% group_by(
  window, alpha
) %>% mutate(
  across(
    matches("^d(Ad|Sc)_p$"),
    ~ p.adjust(.x, method = "fdr"),
    .names = "{sub('_p$', '_q', .col)}"
  )
) %>% ungroup

# Chosen operating point to limit histology imbalance
core <- ext %>% filter(window %in% "imm_3y", alpha == 0.25)

# Save
write_csv(core, "result/table/core.csv")

# 4. Selection of CpGs ─────────────────────────────────────────────────────────
ft <- core %>% filter(
  dAd_p_meta_ds   > 0.05 &
  dSc_p_meta_ds   > 0.05 &
  dPool_p_meta_ds > 0.05
) %>% distinct(
  cpg, .keep_all = TRUE
) %>% arrange(
  dPool_p_meta
)

# Set parameters
params <- list(
  p_org = 0.2,
  p_ax  = 0.05,
  eff   = 0.02,
  universe = nrow(core)
)

# Pair 
pd <- ft %>% mutate(
  keep  = dPool_p_meta < 0.01,
  KCPS  = assign(., "KCPS"),
  NCC   = assign(., "NCC"),
  match = paste(KCPS, NCC, sep = ":"),
  const = {
    bth <- same(dPool_eff_KCPS, dPool_eff_NCC)
    ad  <- same(dAd_eff_KCPS, dAd_eff_NCC)
    sc  <- same(dSc_eff_KCPS, dSc_eff_NCC)
    
    case_when(
      KCPS == "Ad"   & NCC == "Ad"   ~ ad,
      KCPS == "Sc"   & NCC == "Sc"   ~ sc,
      KCPS == "Both" & NCC == "Both" ~ ad & sc,
      KCPS == "Ad"   & NCC == "Both" ~ ad,
      KCPS == "Both" & NCC == "Ad"   ~ ad,
      KCPS == "Sc"   & NCC == "Both" ~ sc,
      KCPS == "Both" & NCC == "Sc"   ~ sc,
      TRUE ~ FALSE
    )
  },
  share = KCPS != "None" & NCC != "None"
)

# Add CHIP
info <- esm %>% filter(chip_group != "any") %>% select(group = chip_group, base = CpG)

pd$base <- str_remove(pd$cpg, "_.*$")
pd %<>% left_join(info, by = "base")

# Rule out insignificant targets
pd %<>% filter(keep)

# 5. Complete Dataset ──────────────────────────────────────────────────────────
ewas <- core %>% filter(cpg %in% pd$cpg)

assigned <- pd %>% select(cpg, base, NCC, KCPS, associated = group, match, const, share)

# Get training set
ids  <- intersect(colnames(m), meta$nested$sentrix_id)
cpgs <- unique(assigned$cpg)

m_tr <- m[cpgs, ids, drop = FALSE]

# Testing set
be <- file.path(data_dir, "beta-raw.rds") %>% read_rds
m  <- logit(be)

m_te <- m[cpgs, , drop = FALSE]

# Save
save(ewas, assigned, m_tr, m_te, file = "object/modeling-set.rda")

# 6. Visualization ─────────────────────────────────────────────────────────────
P_DS   <- 0.05
P_PROX <- 0.01

n_sel <- 7763

# Display only
N_BG_FAIL <- 100000
N_BG_PASS <- 300000

# Number of labels per category in panel B
N_LABEL <- 2

# Reconstruct exact selection universe
u <- core %>%
  distinct(cpg, .keep_all = TRUE) %>%
  mutate(
    p_ds_min = pmin(
      dAd_p_meta_ds,
      dSc_p_meta_ds,
      dPool_p_meta_ds,
      na.rm = FALSE
    ),
    ds_score = -log10(
      pmax(p_ds_min, 1e-300)
    ),
    prox_score = -log10(
      pmax(dPool_p_meta, 1e-300)
    ),
    pass_ds = coalesce(
      dAd_p_meta_ds   > P_DS &
        dSc_p_meta_ds   > P_DS &
        dPool_p_meta_ds > P_DS,
      FALSE
    ),
    selected =
      pass_ds &
      !is.na(dPool_p_meta) &
      dPool_p_meta < P_PROX,
    stage = case_when(
      selected ~ "Final selected",
      pass_ds   ~ "Passed distant gate",
      TRUE      ~ "Excluded by distant gate"
    ),
    stage = factor(
      stage,
      levels = c(
        "Excluded by distant gate",
        "Passed distant gate",
        "Final selected"
      )
    )
  )

bg_fail <- u %>%
  filter(
    stage == "Excluded by distant gate",
    is.finite(ds_score),
    is.finite(prox_score)
  ) %>%
  thin_to(N_BG_FAIL)

bg_pass <- u %>%
  filter(
    stage == "Passed distant gate",
    is.finite(ds_score),
    is.finite(prox_score)
  ) %>%
  thin_to(N_BG_PASS)

fg <- u %>%
  filter(
    selected,
    is.finite(ds_score),
    is.finite(prox_score)
  )

plot_sel <- bind_rows(
  bg_fail,
  bg_pass,
  fg
)

thr_ds   <- -log10(P_DS)
thr_prox <- -log10(P_PROX)

# X-axis transform: moderately compact large values
x_trans <- scales::pseudo_log_trans(base = 10, sigma = 0.5)

# Visual midpoint of the selected region after transformation
x_collected <- x_trans$inverse(
  mean(x_trans$transform(c(0, thr_ds)))
)

p_selection <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = thr_ds, ymin = thr_prox, ymax = Inf,
           fill = "#cfdbe3", alpha = 0.75) +
  geom_point(data = bg_fail, aes(x = ds_score, y = prox_score),
             shape = 16, color = "grey85", size = 0.2) +
  geom_point(data = bg_pass, aes(x = ds_score, y = prox_score),
             shape = 16, color = "grey50", size = 0.2) +
  geom_point(data = fg, aes(x = ds_score, y = prox_score),
             shape = 16, color = "#3a7399", size = 0.2) +
  geom_vline(xintercept = thr_ds, linewidth = 0.2, linetype = "dashed") +
  geom_hline(yintercept = thr_prox, linewidth = 0.2, linetype = "dashed") +
  annotate("text", x = Inf, y = Inf,
           label = expression(paste("Inter-cohort combined ", italic(p))),
           hjust = 1.0, vjust = 1.4, size = 2.15) +
  annotate("text", x = x_collected, y = Inf,
           label = paste0("Collected\nn = ", comma(n_sel)),
           hjust = 0.5, vjust = 1.35, size = 2.15) +
  annotate("text", x = thr_ds, y = -Inf,
           label = expression(italic(p)[min] == 0.05),
           hjust = -0.1, vjust = -0.8, size = 2.15) +
  annotate("text", x = Inf, y = thr_prox,
           label = expression(italic(p)[pool] == 0.01),
           hjust = 1.0, vjust = -0.5, size = 2.15) +
  labs(
    x = expression(
      "Distal " * -log * "[min(" *
        italic(p)[Ad] * ", " *
        italic(p)[Sc] * ", " *
        italic(p)[Pool] * ")]"
    ),
    y = expression(
      "Proximal " * -log * "(" *
        italic(p)[Pool] * ")"
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    trans = x_trans,
    breaks = c(0, 1, 2.5, 5, 10)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  )

p_selection <- p_selection + basic +
  theme(
    # plot.margin = margin(l = 15, r = 5, b = 5, t = 5),
    # plot.tag = element_text(
    #   family = "Arial", face = "bold", size = 16,
    #   margin = margin(l = -15)
    # )
  ) +
  labs(tag = "a")

ggsave(
  "result/select-map.png",
  p_selection,
  width = 2.2,
  height = 1.6,
  dpi = 2000
)

# ── Data ────────────────────────────────────────
net <- read_rds("object/net-modules.rds")

# Data
plot_xy <- fg %>%
  distinct(cpg, .keep_all = TRUE) %>%
  filter(is.finite(dSc_eff_meta), is.finite(dAd_eff_meta)) %>%
  mutate(dAx = abs(dAd_eff_meta) - abs(dSc_eff_meta))

net_cpg <- unique(names(net$colors))
net_pts <- plot_xy %>%
  filter(cpg %in% net_cpg) %>%
  mutate(label = sub("_.*$", "", cpg))

cat("Selected CpGs:", nrow(plot_xy),
    "| Network CpGs:", length(net_cpg),
    "| Matched:", nrow(net_pts), "\n")

# Axis limits: central 95%, including network CpGs
x_lim <- range(c(quantile(plot_xy$dSc_eff_meta, c(0.01, 0.99)),
                 net_pts$dSc_eff_meta), na.rm = TRUE)
y_lim <- range(c(quantile(plot_xy$dAd_eff_meta, c(0.01, 0.99)),
                 net_pts$dAd_eff_meta), na.rm = TRUE)

x_pad <- max(diff(x_lim) * 0.08, 1e-6)
y_pad <- max(diff(y_lim) * 0.08, 1e-6)
x_lim <- x_lim + c(-x_pad, x_pad)
y_lim <- y_lim + c(-y_pad, y_pad)
x_lim[2] <- x_lim[2] + diff(x_lim) * 0.08

# Effect-skew color limits
lim <- as.numeric(quantile(abs(plot_xy$dAx), 0.98, na.rm = TRUE))
if (!is.finite(lim) || lim <= 0) lim <- max(abs(plot_xy$dAx), na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- 1e-6

# KDE
kd <- MASS::kde2d(
  x = plot_xy$dSc_eff_meta, y = plot_xy$dAd_eff_meta,
  n = 500, lims = c(x_lim, y_lim)
)

cloud <- expand.grid(dSc_eff_meta = kd$x, dAd_eff_meta = kd$y) %>%
  mutate(
    density = as.vector(kd$z),
    density_rel = density / max(density),
    alpha_cloud = ifelse(density_rel < 0.002, 0,
                         pmin(0.96, 0.28 + 0.78 * density_rel^0.28)),
    dAx_display = pmax(pmin(abs(dAd_eff_meta) - abs(dSc_eff_meta), lim), -lim)
  )

# Plot
p_skew <- ggplot() +
  geom_raster(
    data = cloud,
    aes(x = dSc_eff_meta, y = dAd_eff_meta,
        fill = dAx_display, alpha = alpha_cloud * 0.8),
    interpolate = TRUE
  ) +
  
  scale_alpha_identity() +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey75") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey75") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.2,
              linetype = "22", color = "grey75") +
  geom_abline(slope = -1, intercept = 0, linewidth = 0.2,
              linetype = "22", color = "grey75") +
  geom_contour(
    data = cloud,
    aes(x = dSc_eff_meta, y = dAd_eff_meta, z = density_rel),
    breaks = c(0.01, 0.10, 0.25, 0.50, 0.75), color = "grey30",
    linewidth = 0.1
  ) +
  # ME-associated annotation (top-left)
  annotate(
    "point",
    x = x_lim[1] + diff(x_lim) * 0.04,
    y = y_lim[2] - diff(y_lim) * 0.07,
    shape = 16, size = 1,
    color = "#BF424A"
  ) +
  annotate(
    "text",
    x = x_lim[1] + diff(x_lim) * 0.065,
    y = y_lim[2] - diff(y_lim) * 0.07,
    label = "'ME-associated' ~ (italic(n) == '5 / 1,000')",
    parse = TRUE,
    hjust = 0, vjust = 0.5,
    size = 1.5, family = "Arial",
    color = "black"
  ) +
  geom_point(
    data = net_pts,
    aes(x = dSc_eff_meta, y = dAd_eff_meta),
    shape = 16, size = 1, color = "#BF424A"
  ) +
  geom_text(
    data = net_pts,
    aes(x = dSc_eff_meta, y = dAd_eff_meta, label = label),
    nudge_x = diff(x_lim) * 0.018,
    hjust = 0, size = 1.5, family = "Arial",
    color = "black", show.legend = FALSE
  ) +
  scale_fill_gradientn(
    colors = c("#467691", "#92AEBB", "#DFC48D", "#c2aa80"),
    # colors = c("#286F68", "#85B6AD", "#C3CDD1", "#88ACD1", "#3267A8"),
    values = scales::rescale(c(-lim, -0.35*lim, 0, 0.35*lim, lim)),
    limits = c(-lim, lim), oob = scales::squish,
    breaks = c(-lim, 0, lim),
    labels = c("Sc", "Equal", "Ad"),
    name = expression("|" * Delta*M[Ad] * "|" - "|" * Delta*M[Sc] * "|"),
    guide = guide_colorbar(
      direction = "vertical", title.position = "right",
      label.position = "right", ticks = TRUE,
      barwidth = unit(1.5, "mm"), barheight = unit(13, "mm")
    )
  ) +
  labs(
    x = expression("Proximal " * Delta*M[Sc]),
    y = expression("Proximal " * Delta*M[Ad]),
    tag = "b"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  basic +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    
    legend.title = element_text(
      family = "Arial", size = 6, angle = -90,
      hjust = 0.5, vjust = 0.5, margin = margin(l = 2, unit = "pt")
    ),
    legend.text = element_text(
      family = "Arial", size = 5.5, angle = -90,
      hjust = 0.5, vjust = 0.5, margin = margin(l = 2, unit = "pt")
    ),
    legend.box.spacing = unit(2, "mm")
  )

ggsave(
  "result/skewed-distribution.svg", device = svglite,
  p_skew,
  width = 2.8,
  height = 1.6,
  dpi = 800
)
