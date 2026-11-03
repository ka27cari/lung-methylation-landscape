# 1. Set Environment ───────────────────────────────────────────────────────────
library(WGCNA)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Load data
load("object/modeling-set.rda")

# 2. Generate Functions ────────────────────────────────────────────────────────
# Breakpoint finder
find_break <- function(time, delta, min_n = 5, top_n = 5) {
  idx <- complete.cases(time, delta)
  
  t <- time[idx]
  d <- delta[idx]
  
  if (length(t) < 2 * min_n || length(unique(t)) < 2) {
    brk <- tibble(
      tau    = double(),
      score  = double(),
      mu_imm = double(),
      mu_dis = double(),
      n_imm  = integer(),
      n_dis  = integer()
    )
    
    return(brk)
  }
  
  o <- order(t)
  t <- t[o]
  d <- d[o]
  
  # Midpoint between adjacent timepoints
  cand <- unique((t[-1] + t[-length(t)]) / 2)
  
  out <- map_dfr(
    cand, function(tau) {
      imm_idx <- t <= tau
      dis_idx <- t > tau
      
      n_imm <- sum(imm_idx, na.rm = TRUE)
      n_dis <- sum(dis_idx, na.rm = TRUE)
      
      if (n_imm < min_n || n_dis < min_n) return(NULL)
      
      d_imm <- d[imm_idx]
      d_dis <- d[dis_idx]
      
      mu_imm <- mean(d_imm, na.rm = TRUE)
      mu_dis <- mean(d_dis, na.rm = TRUE)
      
      sse   <- sum((d_imm - mu_imm) ^ 2, na.rm = TRUE) + sum((d_dis - mu_dis) ^ 2, na.rm = TRUE)
      score <- -sse
      
      tibble(
        tau    = tau,
        score  = score,
        mu_imm = mu_imm,
        mu_dis = mu_dis,
        n_imm  = n_imm,
        n_dis  = n_dis
      )
    }
  )
  
  out %>% slice_max(score, n = top_n, with_ties = FALSE)
}

# Evaluate breakpoint
scan <- function(cpg, axis = c("Pool", "Ad", "Sc"), pair, expr, delta, min_n = 5, top_n = 5) {
  axis <- match.arg(axis)
  
  typ <- switch(
    axis,
    Pool = levels(pair$type),
    Ad   = "LUAD", 
    Sc   = c("LUSC", "SCLC")
  )
  
  use <- pair$type %in% typ
  
  # Breakpoint discovery set
  tr <- use & pair$origin == "KCPS"
  
  bp <- find_break(
    time  = pair$time[tr],
    delta = delta[cpg, tr],
    min_n = min_n,
    top_n = top_n
  )
  
  if (nrow(bp) == 0) {
    return(tibble(
      cpg      = cpg,
      cpg_axis = axis,
      tau      = NA,
      score    = NA,
      gain     = NA
    ))
  }
  
  # Evaluate one breakpoint
  get_eval <- function(idx, tau, ref) {
    t <- pair$time[idx]
    
    pre <- expr[cpg, pair$prediagnostic[idx]]
    con <- expr[cpg, pair$control[idx]]
    
    sgn <- if_else(ref >= 0, 1, -1)
    
    pre <- sgn * pre
    con <- sgn * con
    
    imm_idx <- t <= tau
    dis_idx <- t > tau
    
    tibble(
      n_pair_imm = sum(imm_idx, na.rm = TRUE),
      n_pair_out = sum(dis_idx, na.rm = TRUE),
      auc_imm    = eval_auc(pre[imm_idx], con[imm_idx]),
      auc_dis    = eval_auc(pre[dis_idx], con[dis_idx]),
      mu_imm     = mean(delta[cpg, idx][imm_idx], na.rm = TRUE),
      mu_dis     = mean(delta[cpg, idx][dis_idx], na.rm = TRUE),
      acc_imm    = mean(sign(delta[cpg, idx][imm_idx]) == sign(ref), na.rm = TRUE),
      acc_dis    = mean(sign(delta[cpg, idx][dis_idx]) == sign(ref), na.rm = TRUE)
    )
  }
  
  # Evaluate all candidate breakpoints
  n <- nrow(bp)
  
  out <- map_df(
    1:n, function(i) {
      t   <- bp$tau[i]
      ref <- bp$mu_imm[i]
      
      ori  <- levels(pair$origin)
      metr <- map(
        ori, ~ {
          o <- .x
          
          res <- get_eval(
            idx = use & pair$origin == o,
            tau = t,
            ref = ref
          )
          
          res %>% rename_with(
            ~ paste(tolower(o), .x, sep = "_")
          )
        }
      )
      
      tibble(
        cpg      = cpg,
        cpg_axis = axis,
        tau      = t,
        score    = bp$score[i]
      ) %>% bind_cols(
        metr
      )
    }
  )
  
  # Pick best breakpoint in NCC
  out %<>% mutate(
    kcps_gain      = kcps_auc_imm - kcps_auc_dis,
    ncc_gain       = ncc_auc_imm  - ncc_auc_dis,
    consensus_gain = pmin(kcps_gain, ncc_gain),
    mean_gain      = rowMeans(cbind(kcps_gain, ncc_gain), na.rm = TRUE),
    gain_balance   = abs(kcps_gain - ncc_gain)
  )
  
  # Tau selection by KCPS → evaluate on NCC
  validation <- out %>% arrange(
    desc(score), desc(kcps_gain)
  ) %>% slice(1) %>% mutate(
    selection = "validation",
    gain = ncc_gain
  )
  
  # Tau selection by joint
  consensus <- out %>% arrange(
    desc(consensus_gain), gain_balance, desc(mean_gain), desc(score)
  ) %>% slice(1) %>% mutate(
    selection = "consensus",
    gain = consensus_gain
  )
  
  # Optimized for comparison
  optimized <- out %>% arrange(
    desc(ncc_gain), desc(score)
  ) %>% slice(1) %>% mutate(
    selection = "optimized",
    gain = ncc_gain
  )
  
  bind_rows(validation, consensus, optimized)
}

# CpG set selector
make_set <- function(ax, scope = c("pool", "only"), acc_gap = FALSE, coherent = FALSE) {
  scope <- match.arg(scope)
  axes  <- if (scope != "pool") ax else c("Pool", ax)
  
  # Define columns
  ord     <- paste0(ax, "_ord")
  support <- paste0(ax, "_support")
  coh     <- paste0(ax, "_coh")
  
  out <- agg %>% filter(
    selection == "consensus",
    cpg_axis %in% axes,
    
    valid,
    agree,
    
    gain > 0,
    auc_gap > 0,
    
    .data[[ord]],
    .data[[support]],
    
    ds_neg,
    ctr_neg
  )
  
  if (acc_gap)  out <- out %>% filter(acc_gap >= 0)
  if (coherent) out <- out %>% filter(.data[[coh]])
  
  out %>% distinct(cpg) %>% pull(cpg)
}

# 3. Find Inflection Points ────────────────────────────────────────────────────
pair <- meta$nested %>% filter(
  sentrix_id %in% colnames(m_tr)
) %>% group_by(
  origin, matched_index
) %>% summarise(
  matched_index = matched_index[group == "prediagnostic"][1],
  prediagnostic = sentrix_id[group == "prediagnostic"][1],
  control       = sentrix_id[group == "control"][1],
  time          = time[group == "prediagnostic"][1],
  type          = type[group == "prediagnostic"][1],
  .groups = "drop"
)

d <- m_tr[, pair$prediagnostic] - m_tr[, pair$control]; colnames(d) <- pair$matched_index
d[is.na(d)] <- 0

# Save
save(pair, d, file = "object/pair-src.rda")

# Find breakpoint
n <- length(ewas$cpg)

jobs <- crossing(
  idx  = seq_len(n), 
  axis = c("Pool", "Ad", "Sc")
)

bps <- pmap_dfr(
  jobs, function(idx, axis) {
    scan(
      cpg    = ewas$cpg[idx],
      axis   = axis,
      pair   = pair,
      expr   = m_tr,
      delta  = d,
      min_n  = 10,
      top_n  = 5
    )
  }, .progress = TRUE
)

# Save
write_csv(bps, "result/table/breakpoints.csv")

# Make consensus CpG sets
bps %<>% mutate(
  gain     = pmin(kcps_gain, ncc_gain, na.rm = TRUE),
  
  auc_imm  = pmin(kcps_auc_imm, ncc_auc_imm, na.rm = TRUE),
  auc_dis  = pmax(kcps_auc_dis, ncc_auc_dis, na.rm = TRUE),
  auc_gap  = auc_imm - auc_dis,
  
  acc_kcps = kcps_acc_imm - kcps_acc_dis,
  acc_ncc  = ncc_acc_imm - ncc_acc_dis,
  acc_gap  = pmin(acc_kcps, acc_ncc, na.rm = TRUE),
  
  mu       = rowMeans(
    cbind(kcps_mu_imm, ncc_mu_imm), na.rm = TRUE
  ),
  dir      = sign(mu),
  valid    = is.finite(tau) & is.finite(dir) & dir != 0,
  agree    = 
    is.finite(kcps_mu_imm) &
    is.finite(ncc_mu_imm) &
    sign(kcps_mu_imm) == sign(ncc_mu_imm) &
    sign(kcps_mu_imm) != 0 &
    sign(ncc_mu_imm) != 0
)

# Get information
ewas %<>% distinct(cpg, .keep_all = TRUE) %>% select(-window, -alpha)

assigned %<>% distinct(cpg, .keep_all = TRUE) %>% select(cpg, KCPS, NCC, const, share)

agg <- bps %>% left_join(
  assigned, by = "cpg"
) %>% left_join(
  ewas, by = "cpg"
) %>% mutate(
  ds_neg  = 
    dAd_p_meta_ds   > 0.05 &
    dSc_p_meta_ds   > 0.05 &
    dPool_p_meta_ds > 0.05,
  ctr_neg = 
    cCtr_p_meta > 0.05 |
    abs(cCtr_eff_meta) < 0.10,
  
  lp_Ad   = -log10(dAd_p_meta),
  lp_Sc   = -log10(dSc_p_meta),
  lp_Pool = -log10(dPool_p_meta),
  
  Ad_coh  = case_when(
    cpg_axis == "Pool" ~ KCPS != "None" & NCC != "None",
    cpg_axis == "Ad"   ~ KCPS %in% c("Ad", "Both") & NCC %in% c("Ad", "Both"),
    TRUE ~ FALSE
  ),
  
  Sc_coh  = case_when(
    cpg_axis == "Pool" ~ KCPS != "None" & NCC != "None",
    cpg_axis == "Sc"   ~ KCPS %in% c("Sc", "Both") & NCC %in% c("Sc", "Both"),
    TRUE ~ FALSE
  ),
  
  Ad_support = lp_Ad > 1 | lp_Pool > 2,
  Sc_support = lp_Sc > 1 | lp_Pool > 2
)

# Tau-stage order
nst <- meta$nested %>% filter(
  sentrix_id %in% colnames(m_tr)
) %>% mutate(
  axis = case_when(
    type == "LUAD" ~ "Ad",
    type %in% c("LUSC", "SCLC") ~ "Sc",
    TRUE ~ NA
  )
)

ids  <- nst %>% filter(group == "control") %>% pull(sentrix_id)
n    <- nrow(agg)
axis <- c("Ad", "Sc")

agg$Ad_ord <- FALSE 
agg$Sc_ord <- FALSE

for (i in 1:n) {
  cg  <- agg$cpg[i]
  tau <- agg$tau[i]
  dir <- agg$dir[i]

  if (!is.finite(tau) || !is.finite(dir) || dir == 0) next
  
  for (ax in axis) {
    pre <- nst %>% filter(group == "prediagnostic", axis == ax)
    
    imm <- pre %>% filter(time <= tau) %>% pull(sentrix_id)
    dis <- c(pre %>% filter(time > tau) %>% pull(sentrix_id), ids)
    
    if (length(imm) == 0 || length(dis) == 0) next
    
    hi <- dir * as.numeric(m_tr[cg, imm])
    lo <- dir * as.numeric(m_tr[cg, dis])
    
    x_auc  <- eval_auc(hi, lo)
    x_diff <- mean(hi, na.rm = TRUE) - mean(lo, na.rm = TRUE)
    
    pass <- is.finite(x_auc) && x_auc > 0.5 && is.finite(x_diff) && x_diff > 0
    
    if (ax == "Ad") agg$Ad_ord[i] <- pass
    if (ax == "Sc") agg$Sc_ord[i] <- pass
  }
}

# Save
write_csv(agg, "result/table/aggregate.csv")

# Generate CpG set
spec <- list(
  stage    = list(scope = "pool"),
  acc_gap  = list(scope = "pool", acc_gap = TRUE),
  only     = list(scope = "only"),
  coherent = list(scope = "pool", coherent = TRUE),
  robust   = list(scope = "pool", acc_gap = TRUE, coherent = TRUE)
)

cg_set <- map(
  spec, ~ list(
    Ad = exec(make_set, "Ad", !!!.x),
    Sc = exec(make_set, "Sc", !!!.x)
  )
)

# Save
write_rds(cg_set, "object/cg-set.rds")
