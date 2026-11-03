# 1. Set Environment ───────────────────────────────────────────────────────────
require(pROC)
require(PRROC)
library(WGCNA)
# library(tidymodels)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Load data
hr  <- read_rds("object/hazard-risk.rds")
net <- read_rds("object/net-modules.rds")

meta %<>% map2(hr, left_join)

m <- file.path(data_dir, "m-adjusted.rds") %>% read_rds

# 2. Generate Functions ────────────────────────────────────────────────────────
# Random combination generator
combinate <- function(vars, n_draw = 300, k_min = 1, k_max = 4, seed = 1) {
  set.seed(seed)
  k_max <- min(k_max, length(vars))
  
  cb <- replicate(
    n_draw, {
      k <- sample(k_min:k_max, 1)
      sort(sample(vars, k))
    }, simplify = FALSE
  )
  
  # Erase duplicates
  key <- vapply(cb, paste, collapse = "+", FUN.VALUE = character(1))
  cb[!duplicated(key)]
}

# Fit safe GLM
safe_glm <- function(fml, data) {
  tryCatch(
    suppressWarnings(glm(fml, data = data, family = binomial)),
    error = function(e) NULL
  )
}

# Make folds
make_folds <- function(y, max_k = 5, seed = 1) {
  set.seed(seed)
  
  fid <- integer(length(y))
  
  idx1 <- which(y != 0)
  idx0 <- which(y == 0)
  
  fid[idx1] <- sample(rep(1:max_k, length.out = length(idx1)))
  fid[idx0] <- sample(rep(1:max_k, length.out = length(idx0)))
  
  return(fid)
}

# Compute CV score
cv_score <- function(data, formula, max_k = 5, seed = 1) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  
  cv <- list(ro_auc = NA, pr_auc = NA, n = nrow(mf))
  if (nrow(mf) < max_k * 5) return(cv)
  
  # Split folds
  fid <- make_folds(mf$event, max_k = max_k, seed = seed)
  ro_auc <- numeric(max_k); pr_auc <- numeric(max_k)
  
  for (k in 1:max_k) {
    tr <- mf[fid != k, , drop = FALSE]; te <- mf[fid == k, , drop = FALSE]
    if (length(unique(tr$event)) < 2 || length(unique(te$event)) < 2) return(cv)
    
    fit <- safe_glm(formula, tr)
    if (is.null(fit)) next
    
    p <- suppressWarnings(predict(fit, newdata = te, type = "response"))
    
    ro_auc[k] <- auc(roc(te$event, p, quiet = TRUE, direction = "<")) %>% as.numeric
    pr_auc[k] <- pr.curve(
      scores.class0 = p[te$event == 1], 
      scores.class1 = p[te$event == 0], 
      curve = FALSE
    )$auc.integral
  }
  
  cv$ro_auc <- mean(ro_auc, na.rm = TRUE)
  cv$pr_auc <- mean(pr_auc, na.rm = TRUE)
  
  return(cv)
}

# Performance evaluation
evaluate <- function(fit, data, name) {
  mf <- model.frame(formula(fit), data = data, na.action = na.omit)
  
  event <- model.response(mf)
  pred  <- predict(fit, newdata = mf, type = "response")
  
  keep <- complete.cases(event, pred)
  
  event <- event[keep]
  pred  <- pred[keep]
  
  if (length(unique(event)) < 2) list(
    model = name,
    roc = tibble(fpr = numeric(), tpr = numeric(), model = character()),
    pr  = tibble(recall = numeric(), precision = numeric(), threshold = numeric(), model = character()),
    n   = length(event),
    ro_auc = NA,
    pr_auc = NA
  )
  
  rocs0  <- roc(response = event, predictor = pred, quiet = TRUE, direction = "<")
  au_roc <- as.numeric(auc(rocs0))
  
  rocs <- data.frame(fpr = 1 - rocs0$specificities, tpr = rocs0$sensitivities, model = name)
  
  prs0 <- pr.curve(
    scores.class0 = pred[event != 0],
    scores.class1 = pred[event == 0],
    curve = TRUE
  )
  
  pr_auc <- as.numeric(prs0$auc.integral)
  prs <- as.data.frame(prs0$curve)
  
  colnames(prs) <- c("recall", "precision", "threshold")
  prs$model <- name
  
  list(
    model = name,
    roc = rocs,
    pr  = prs,
    n   = length(event),
    ro_auc = au_roc,
    pr_auc = pr_auc
  )
}

# Run tidymodel
run_model <- function(name, wf, grid, folds, m_set, x_te, x_va) {
  cv <- tune_grid(
    wf,
    resamples = folds,
    grid      = grid,
    metrics   = m_set
  )
  
  best <- select_best(cv, metric = "roc_auc")
  mod  <- finalize_workflow(wf, best) %>% fit(x_te)
  
  # Predict based on best model
  pred <- predict(mod, x_va, type = "prob") %>% bind_cols(event = x_va$event)
  
  roc <- roc_curve(
    pred, truth = event, .pred_1, event_level = "second"
  ) %>% transmute(
    fpr   = 1 - specificity,
    tpr   = sensitivity,
    model = name
  )
  
  pr <- pr_curve(
    pred, truth = event, .pred_1, event_level = "second"
  ) %>% transmute(
    recall    = recall,
    precision = precision,
    threshold = .threshold,
    model     = name
  )
  
  out <- tibble(
    k      = NA,
    n_used = nrow(x_te),
    model  = name,
    roc    = list(roc),
    pr     = list(pr),
    n      = nrow(x_va),
    ro_auc = roc_auc(pred, truth = event, .pred_1, event_level = "second")$.estimate,
    pr_auc = pr_auc(pred, truth = event, .pred_1, event_level = "second")$.estimate
  )
  
  list(cv = cv, best = best, fit = mod, pred = pred, out = out)
}

# 3. Assign Criteria ───────────────────────────────────────────────────────────
coh <- meta$nested %>% select(
  origin, matched_index, group, time, sentrix_id, age, sex, starts_with("smk"), ces_dur, type
) %>% mutate(
  event     = if_else(group != "control", 1L, 0L),
  ever      = if_else(smk_status != "none", 1L, 0L),
  current   = if_else(smk_status == "current", 1L, 0L),
  former    = if_else(smk_status == "quit", 1L, 0L),
  smk_amnt  = if_else(smk_status != "none", smk_amnt, 0),
  smk_dur   = if_else(smk_status != "none", smk_dur, 0),
  packyears = if_else(smk_status != "none", smk_dur * smk_amnt / 20, 0),
  ces_dur   = case_when(
    smk_status == "none"    ~ NA, 
    smk_status == "current" ~ 0, 
    TRUE                    ~ ces_dur
  )
)

# Screening criteria
coh %<>% mutate(
  NLST   = if_else(
    age >= 55 & age <= 74 & packyears >= 30 & (
      current == 1L | (former == 1L & ces_dur <= 15)
    ), 1L, 0L, missing = NA
  ),
  USPSTF = if_else(
    age >= 50 & age <= 80 & packyears >= 20 & (
      current == 1L | (former == 1L & ces_dur <= 15)
    ), 1L, 0L, missing = NA
  ),
  KNLCS  = if_else(
    age >= 55 & age <= 74 & packyears >= 30 & (
      current == 1L | (former == 1L & ces_dur <= 15)
    ), 1L, 0L, missing = NA
  ),
  NELSON = if_else(
    age >= 50 & age <= 75 & (
      current == 1L | (former == 1L & ces_dur <= 10)
    ) & (
      (smk_amnt > 15 & smk_dur > 25) | (smk_amnt > 10 & smk_dur > 30)
    ), 1L, 0L, missing = NA
  ),
  NCCN   = if_else(
    age >= 50 & packyears >= 20 & (
      current == 1L | (former == 1L & ces_dur <= 15)
    ), 1L, 0L, missing = NA
  )
)

# Compute eigengenes
ids  <- colnames(m); cpgs <- names(net$colors)
expr <- m[cpgs, ] %>% impute %>% t

egs <- moduleEigengenes(expr, colors = net$colors)$eigengenes
egs$sentrix_id <- ids

# Redefine eigengenes
egs %<>% as_tibble %>% relocate(sentrix_id, .before = 0)
coh %<>% left_join(egs)

# Split set
set <- coh %>% filter(sentrix_id %in% ids) %>% split(., .$origin)

tr <- set$NCC; tv <- set$KCPS

te <- tr
va <- tv %>% filter(smk_status != "quit")

# 4. Module based Prediction ───────────────────────────────────────────────────
me_vars <- grep("^ME\\d+$", names(te), value = TRUE); me_vars %<>% setdiff(., "ME0")
b <- combinate(vars = me_vars, n_draw = 400, k_min = 1, k_max = 3, seed = 1)

terms <- b %>% append(
  list(
    "smk_status",
    c("age", "sex", "smk_status"),
    c("age", "sex", "smk_status", "packyears")
  )
)

label <- vapply(
  b, paste, collapse = "+", FUN.VALUE = character(1)
) %>% c(
  "reduced", "core", "saturated"
)

specs <- tibble(
  terms = terms,
  label = label,
  k     = lengths(terms),
  formula = map(
    terms, ~ as.formula(paste("event ~", paste(.x, collapse = " + ")))
  )
)

# Fit data with test samples only
specs %<>% mutate(
  mf     = map(
    formula, ~ model.frame(.x, data = te, na.action = na.omit)
  ),
  n_used = map_int(mf, nrow),
  fit    = map(
    formula, ~ safe_glm(.x, te)
  ),
  status = case_when(
    n_used < 20           ~ "few_n",
    map_lgl(fit, is.null) ~ "failed",
    TRUE ~ "ok"
  )
)

# Evaluate CV in test samples
specs %<>% filter(status == "ok")

cands <- specs %>% mutate(
  score  = map(
    formula, ~ cv_score(te, .x, max_k = 5, seed = 1)
  ),
  ro_auc = map_dbl(score, "ro_auc"),
  pr_auc = map_dbl(score, "pr_auc")
) %>% select(
  formula, label, k, ro_auc, pr_auc, n_used
)

# Settle up with top models
set  <- cands %>% arrange(desc(ro_auc)) %>% slice_head(n = 20)
base <- cands %>% filter(!str_detect(label, "^ME"))

set %<>% bind_rows(base) %>% pull(label)

# Evaluate in validation set
eva <- specs %>% filter(
  label %in% set
) %>% mutate(
  .auc = map2(
    fit, label, ~ evaluate(.x, va, .y)
  )
) %>% select(
  k, n_used, .auc
) %>% unnest_wider(
  .auc
)

# 5. Machine-learning based Prediction ─────────────────────────────────────────
tidymodels_prefer()

x_te <- t(m[cpgs, te$sentrix_id]); x_va <- t(m[cpgs, va$sentrix_id])

x_te %<>% as_tibble %>% mutate(event = factor(te$event), .before = 1)
x_va %<>% as_tibble %>% mutate(event = factor(va$event), .before = 1)

# Resampling / metrics
folds <- vfold_cv(x_te, v = 5, strata = event)
m_set <- metric_set(roc_auc, pr_auc)

# Try Elastic Net
rec <- recipe(
  event ~ ., data = x_te
) %>% step_impute_median(
  all_predictors()
) %>% step_zv(
  all_predictors()
) %>% step_normalize(
  all_predictors()
)

spec <- logistic_reg(
  penalty = tune(),
  mixture = tune()
) %>% set_engine(
  "glmnet"
) %>% set_mode(
  "classification"
)

wf <- workflow() %>% add_recipe(rec) %>% add_model(spec)

# Tune grids
grid <- grid_regular(
  mixture(range = c(0, 1)), penalty(range = c(-6, 0)),
  levels = c(mixture = 5, penalty = 20)
)

# Add result
en <- run_model(
  name = "elastic_net",
  wf   = wf,
  grid = grid,
  folds = folds,
  m_set = m_set,
  x_te  = x_te,
  x_va  = x_va
)$out

eva %<>% bind_rows(en)

# Try XGBoost
rec <- recipe(
  event ~ ., data = x_te
) %>% step_impute_median(
  all_predictors()
) %>% step_zv(
  all_predictors()
)

spec <- boost_tree(
  trees = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  mtry  = tune(),
  min_n = tune(),
  sample_size = tune()
) %>% set_engine(
  "xgboost", eval_metric = "auc", tree_method = "hist"
) %>% set_mode(
  "classification"
)

wf <- workflow() %>% add_recipe(rec) %>% add_model(spec)

# Tune grids
grid <- grid_space_filling(
  trees(range = c(300, 1500)),
  tree_depth(range = c(2, 6)),
  learn_rate(range = c(-4, -1)),
  finalize(mtry(), x_te),
  min_n(range = c(1, 30)),
  sample_prop(range = c(0.5, 1.0)),
  size = 30
)

# Add result
xgb <- run_model(
  name = "xgboost",
  wf   = wf,
  grid = grid,
  folds = folds,
  m_set = m_set,
  x_te  = x_te,
  x_va  = x_va
)$out

eva %<>% bind_rows(xgb)

# 6. Compare Everything ────────────────────────────────────────────────────────
me   <- eva %>% filter(str_detect(model, "^ME\\d"))
best <- me %>% arrange(desc(ro_auc)) %>% slice(1) %>% pull(model)

eva %<>% filter(!model %in% me$model | model %in% best)

# Save
write_rds(eva, "object/roc-eval.rds")
eva <- read_rds("object/roc-eval.rds")

# Model labels and colors
legend_order <- c(
  "3-ME signature",
  "Model 1: smoking status",
  "Model 2: Model 1 + age + sex",
  "Model 3: Model 2 + packyears",
  "Elastic Net",
  "XGBoost"
)

roc_cols <- c(
  "Model 1: smoking status" = "#CDD2D4",
  "Model 2: Model 1 + age + sex" = "#BDC3C6",
  "Model 3: Model 2 + packyears" = "#AEB5B9",
  "Elastic Net" = "#247E70",
  "XGBoost" = "#12554F",
  "3-ME signature" = "#2a5a7a"
)

# All model ROC curves
roc_all <- eva %>%
  select(roc) %>%
  unnest(roc) %>%
  filter(
    model %in% c(
      best, "reduced", "core", "saturated",
      "elastic_net", "xgboost"
    ),
    is.finite(fpr),
    is.finite(tpr)
  ) %>%
  mutate(
    label = case_when(
      model == best ~ "3-ME signature",
      model == "reduced" ~ "Model 1: smoking status",
      model == "core" ~ "Model 2: Model 1 + age + sex",
      model == "saturated" ~ "Model 3: Model 2 + packyears",
      model == "elastic_net" ~ "Elastic Net",
      model == "xgboost" ~ "XGBoost",
      TRUE ~ model
    ),
    label = factor(label, levels = legend_order)
  ) %>%
  arrange(model, fpr, tpr)

# Main model ROC
main_roc <- roc_all %>%
  filter(model == best) %>%
  arrange(fpr, tpr)

other_roc <- roc_all %>%
  filter(model != best)

stopifnot(nrow(main_roc) >= 2)

# AUROC from stored ROC coordinates
main_auc <- with(
  main_roc,
  sum(
    diff(fpr) *
      (head(tpr, -1) + tail(tpr, -1)) / 2
  )
)

stopifnot(is.finite(main_auc))

# ACS 2023 screening criteria
te <- te %>%
  mutate(
    ACS2023 = case_when(
      is.na(age) |
        is.na(packyears) |
        is.na(smk_status) ~ NA_integer_,
      TRUE ~ as.integer(
        between(age, 50, 80) &
          packyears >= 20 &
          smk_status %in% c("current", "quit")
      )
    )
  )

# Screening criteria
proto <- c(
  "NLST", "USPSTF", "KNLCS", "NELSON", "ACS2023"
)

stopifnot(all(c("event", proto) %in% names(te)))

safe_div <- function(a, b) {
  if (b == 0) NA_real_ else a / b
}

scrn <- purrr::map_dfr(proto, \(p) {
  data <- tibble(
    event = te$event,
    rule = te[[p]]
  ) %>%
    filter(!is.na(event), !is.na(rule))
  
  tp <- sum(data$rule == 1 & data$event == 1)
  fp <- sum(data$rule == 1 & data$event == 0)
  fn <- sum(data$rule == 0 & data$event == 1)
  tn <- sum(data$rule == 0 & data$event == 0)
  
  tibble(
    system = p,
    fpr = safe_div(fp, fp + tn),
    tpr = safe_div(tp, tp + fn),
    n = nrow(data)
  )
})

# Screening labels: upper-left area
scrn_label <- scrn %>%
  filter(is.finite(fpr), is.finite(tpr)) %>%
  arrange(tpr, fpr) %>%
  mutate(
    label_x = 0.27,
    label_y = seq(0.97, 0.71, length.out = n()),
    display = if_else(
      system == "ACS2023", "ACS", system
    )
  )

# Plot ROC
p_e <- ggplot(
  roc_all,
  aes(x = fpr, y = tpr, color = label, group = model)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "#C7CDD0",
    linetype = "22",
    linewidth = 0.15
  ) +
  
  # Other model ROC curves
  geom_step(
    data = other_roc,
    direction = "hv",
    linewidth = 0.15,
    alpha = 0.5,
    lineend = "round"
  ) +
  
  # Main model ROC curve (drawn on top)
  geom_step(
    data = main_roc,
    direction = "hv",
    linewidth = 0.35,
    alpha = 1,
    lineend = "round"
  ) +
  
  # Screening criterion connectors
  geom_segment(
    data = scrn_label,
    aes(
      x = fpr,
      y = tpr,
      xend = label_x - 0.015,
      yend = label_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.15,
    color = "black",
    show.legend = FALSE
  ) +
  
  # Screening criterion labels
  geom_text(
    data = scrn_label,
    aes(
      x = label_x,
      y = label_y,
      label = display
    ),
    inherit.aes = FALSE,
    size = 1.8,
    family = "Arial",
    color = "black",
    hjust = 0,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  
  # Main model AUROC
  annotate(
    "text",
    x = 0.97,
    y = 0.08,
    label = sprintf("AUROC = %.2f", main_auc),
    hjust = 1,
    vjust = 0.5,
    size = 1.8,
    family = "Arial",
    color = "black"
  ) +
  
  scale_color_manual(
    values = roc_cols,
    breaks = legend_order,
    name = NULL,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keywidth = unit(3, "mm"),
      keyheight = unit(1, "mm"),
      override.aes = list(
        linewidth = 0.30,
        alpha = 1
      )
    )
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(add = 0)
  ) +
  labs(
    x = "False positive rate",
    y = "Sensitivity",
    tag = "e"
  ) +
  basic +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "top",
    legend.box.just = "top",
    legend.key.spacing.x = unit(1, "mm"),
    legend.key.spacing.y = unit(0.5, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

ggsave(
  "result/roc-eval.svg",
  plot = p_e,
  device = svglite::svglite,
  width = 3.2,
  height = 1.5,
  dpi = 2000,
  bg = "transparent"
)

# 6. Assign Class Probability ──────────────────────────────────────────────────
fit <- specs %>% filter(label %in% best) %>% pull(fit) %>% .[[1]]

# Predict by model
coh %<>% mutate(
  .hat  = predict(fit, newdata = ., type = "response"),
  .true = .data[["event"]],
  .pred = if_else(.hat > 0.5, 1L, 0L),
  hit   = .true == .pred,
)

# Save object
write_rds(coh, "object/obs-pred.rds")
coh <- read_rds("object/obs-pred.rds")

# Draw rates
recall <- coh %>% filter(
  !is.na(type), !is.na(hit), origin != "NCC"
) %>% group_by(
  type
) %>% summarise(
  n = n(),
  n_right = sum(hit, na.rm = TRUE),
  n_wrong = sum(!hit, na.rm = TRUE),
  .groups = "drop"
) %>% mutate(
  right_rate = n_right / n,
  wrong_rate = n_wrong / n
) %>% mutate(
  type = fct_reorder(type, right_rate, .desc = TRUE)
) %>% pivot_longer(
  cols = c(n_right, n_wrong),
  names_to  = "status",
  values_to = "count"
) %>% mutate(
  status = recode(status, n_right = "Correct", n_wrong = "Incorrect"),
  status = factor(status, levels = c("Incorrect", "Correct")),
  rate = count / n
) %>% ggplot(
  aes(x = type, y = rate, fill = status)
) + geom_col(
  width = 0.8
) + scale_fill_manual(
  values = c(
    "Correct"   = alpha("#405E73", 0.85),
    "Incorrect" = alpha("#a3a3a3", 0.25)
  )
) + scale_y_continuous(
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  labels = percent_format(accuracy = 1),
  expand = expansion(add = c(0, 0.05))
) + labs(
  x = NULL, y = "Proportion", fill = NULL
) + basic + theme(
  # plot.margin      = margin(l = 10, r = 5, b = 5, t = 10),
  # plot.tag         = element_text(
  #   family = "Arial", face = "bold", size = 16, margin = margin(l = -5, t = -10)
  # ),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
  legend.position  = "top",
  legend.direction = "horizontal",
  legend.box       = "horizontal",
  legend.title = element_blank(),
  legend.text = element_text(size = 6, margin = margin(l = 2)),
  
  legend.key.height = unit(2.3, "mm"),
  legend.key.width = unit(2.3, "mm"),
  
  legend.key.spacing.x = unit(2, "mm"),
  legend.key.spacing.y = unit(0, "mm"),
  legend.margin = margin(0, 0, 0, b = -4)
) + labs(
  tag = "f"
)

# Save
ggsave("result/recall.svg", device = svglite, recall, width = 1.8, height = 1.5, dpi = 800)

# Compare LUSC/SCLC to rest
prog <- coh %>% filter(
  origin != "NCC",
  event  != 0,
  !is.na(type),
  !is.na(time),
  !is.na(.hat)
) %>% mutate(
  histology = if_else(
    type %in% c("LUSC", "SCLC"), "LUSC/SCLC", "Rest"
  ),
  time_bin = case_when(
    time > 10 ~ ">10",
    time > 5  ~ "10–5",
    time > 2  ~ "5–2",
    TRUE      ~ "<2"
  ),
  time_bin = factor(
    # time_bin, levels = c(">10", "10–5", "5–2", "<2")
    time_bin, levels = c("<2", "5–2", "10–5", ">10")
  )
) %>% group_by(
  histology, time_bin
) %>% summarise(
  risk = median(.hat, na.rm = TRUE),
  q25  = quantile(.hat, 0.25, na.rm = TRUE),
  q75  = quantile(.hat, 0.75, na.rm = TRUE),
  n    = n(),
  .groups = "drop"
)

pd <- position_dodge(width = 0.2)

prog <- ggplot(
  prog,
  aes(
    x = time_bin,
    y = risk,
    group = histology,
    color = histology
  )
) + geom_errorbar(
  aes(
    ymin = q25,
    ymax = q75
  ),
  width     = 0.1,
  linewidth = 0.25,
  position  = pd
) + geom_point(
  aes(size = n),
  position = pd
) + geom_line(
  linewidth = 0.25,
  position  = pd
) + scale_color_manual(
  values = c(
    "LUSC/SCLC" = "#eb6d5f",
    "Rest"      = "#9d9d9d"
  ),
  name = NULL
) + scale_y_continuous(
  limits = c(0.25, 0.75),
  breaks = c(0.25, 0.5, 0.75),
  labels = percent_format(accuracy = 1),
  expand = expansion(add = 0)
) + scale_size_continuous(
  range = c(0.2, 1.0),
  guide = "none"
) + labs(
  x = "Time-to-diagnosis (years)",
  y = "Pred. risk"
)

prog <- prog + basic + theme(
  # plot.margin     = margin(l = 12, r = 5, b = 5, t = 10),
  # plot.tag        = element_text(
  #   family = "Arial", face = "bold", size = 16, margin = margin(l = -10, t = -10)
  # ),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.box       = "horizontal",
  legend.text = element_text(size = 6, margin = margin(l = 2)),
  legend.key.height = unit(2, "mm"),
  legend.key.width = unit(2, "mm"),
  
  legend.key.spacing.x = unit(2, "mm"),
  legend.key.spacing.y = unit(0, "mm"),
  legend.margin = margin(0, 0, 0, b = -4)
  
) + labs(
  tag = "g"
)

# Save
ggsave("result/progress.svg", device = svglite, prog, width = 1.8, height = 1, dpi = 800,
  bg = "transparent")
