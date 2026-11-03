# 1. Set Environment ───────────────────────────────────────────────────────────
library(slider)
library(survival)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Load data
hr <- read_rds("object/hazard-risk.rds")

meta %<>% map2(hr, left_join)

# 2. Generate Functions ────────────────────────────────────────────────────────
# Arrange conditional logistic
statistics <- function(fit, type = NULL, analysis = NULL) {
  t <- tidy(fit) %>% filter(term %in% "risk_score")
  
  tibble(
    type     = type, 
    analysis = analysis,
    or       = exp(t$estimate),
    lower    = exp(t$estimate - 1.96 * t$std.error),
    upper    = exp(t$estimate + 1.96 * t$std.error),
    p_value  = t$p.value,
    star     = p_to_star(t$p.value)
  )
}

# 3. Evaluate Hazard Ratio ─────────────────────────────────────────────────────
ncc <- meta$nested %>% filter(origin != "KCPS")

ncc$risk_score <- scale(ncc$risk_score)[, 1]

pre <- meta$nested %>% filter(
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
) %>% select(matched_index, time_bin, type)

pre$type %<>% fct_recode(
  "LUAD"  = "Adenocarcinoma (incl. variants)",
  "LUSC"  = "Squamous cell carcinoma",
  "SCLC"  = "Small cell carcinoma",
  "NSCLC" = "Non-small cell lung cancer",
  "Other" = "Other (rare, unspecific)"
)

pre %<>% select(matched_index, time_bin, type) %>% distinct

# Make delta
d <- ncc %>% select(
  matched_index, group, risk_score
) %>% left_join(
  pre, by = "matched_index"
) %>% group_by(
  matched_index
) %>% summarise(
  time_bin = time_bin[1],
  delta    = risk_score[group != "control"][1] - risk_score[group != "prediagnostic"][1],
  .groups = "drop"
) %>% filter(
  is.finite(delta), !is.na(time_bin)
)

# Compute significance on non-tie comparisons
agg <- d %>% group_by(
  time_bin
) %>% summarise(
  n_pairs = n(),
  n_eff   = sum(delta != 0),
  n_win   = sum(delta > 0),
  p_value = ifelse(n_eff > 0, binom.test(n_win, n_eff, p = 0.5)$p.value, NA),
  star    = p_to_star(p_value),
  .groups = "drop"
)

# Asterisk position
pos <- d %>% group_by(time_bin) %>% summarise(y = max(delta, na.rm = TRUE) + 0.5, .groups = "drop")
agg %<>% left_join(pos) %>% mutate(group1 = time_bin, group2 = time_bin)

trend <- ggplot(
  d, aes(x = time_bin, y = delta, fill = time_bin)
) + geom_hline(
  yintercept  = 0,
  linetype    = "dotted",
  linewidth   = 0.35,
  color       = "#c6c6c6",
  show.legend = FALSE
) + geom_violin(
  trim = TRUE,
  width = 0.7,
  alpha = 0.50,
  color = NA,
  show.legend = FALSE
) + geom_boxplot(
  width = 0.2,
  alpha = 0.9,
  linewidth = 0.25,
  outlier.shape = NA,
  show.legend = FALSE
) + geom_point(
  size  = 0.5,
  alpha = 0.35,
  shape = 16,
  position = position_jitter(width = 0.05, height = 0),
  show.legend = FALSE
) + scale_fill_manual(
  values = setNames(
    sequential_hcl(4, palette = "Teal", rev = TRUE),
    c(">10", "5-10", "2-5", "<2")
  )
) + labs(
  x = "Time-to-diagnosis (years)",
  y = expression(Delta~"Risk score")
) + stat_pvalue_manual(
  agg,
  x     = "time_bin",
  label = "star",
  size  = 2.5,
  y.position = "y",
  tip.length = 0
)

trend <- trend + basic + labs(tag = "e")

# Save
ggsave("result/or-trend.svg", device = svglite, trend, width = 1.65, height = 1.6, dpi = 800, bg = "transparent")

# Data for morphological study
morph <- ncc %>% select(matched_index, group, risk_score)
morph %<>% left_join(pre, by = "matched_index")

# Compute risk regardless of morphology
morph %<>% mutate(event = if_else(group != "control", 1L, 0L))

tot <- clogit(event ~ risk_score + strata(matched_index), data = morph) %>% statistics(type = "All", analysis = "Case vs. matched control") 
tot$n_case <- sum(morph$event)

# Compute risk compared to control by morphology
sat <- map_dfr(
  levels(morph$type), ~ {
    data <- morph %>% filter(type %in% .x)
    
    out <- clogit(event ~ risk_score + strata(matched_index), data = data) %>% statistics(type = .x, analysis = "Case vs. matched control") 
    out$n_case <- sum(data$event)
    
    return(out)
  }
)

# Compare one type to others
red <- map_dfr(
  levels(pre$type), ~ {
    data <- pre %>% mutate(y = type %in% .x); data$y %<>% as.integer
    sc   <- filter(morph, event != 0L)
    
    data %<>% left_join(sc, by = "matched_index")
    
    out <- data %>% glm(y ~ risk_score, data = ., family = binomial) %>% statistics(type = .x, analysis = "In tumors")
    out$n_case <- sum(data$y)
    
    return(out)
  }
)

# Combine all
main <- bind_rows(sat, red)

main$type %<>% fct_reorder(., main$n_case, .desc = FALSE) %>% fct_relevel("Other", after = 0)
main$analysis %<>% fct_relevel("Case vs. matched control")

# Remove non-significants
main$star <- if_else(main$star == "n.s", NA, main$star)

# Morphology upon risk score
dodge <- position_dodge(width = 0.65)

const <- ggplot(
  main, aes(x = or, y = type, color = analysis, shape = analysis)
) + geom_vline(
  xintercept = 1,
  linetype   = "dotted",
  linewidth  = 0.35,
  color      = "#c6c6c6"
) + geom_xsidevline(
  xintercept = 1,
  linetype   = "dotted",
  linewidth  = 0.35,
  color      = "#c6c6c6"
) + geom_point(
  aes(size = n_case),
  position = dodge
) + geom_xsidepoint(
  data = tot,
  aes(size = n_case)
) + geom_segment(
  aes(x = lower, xend = upper),
  position = dodge,
  linewidth = 0.35,
  lineend = "round"
) + geom_xsidesegment(
  data = tot,
  aes(x = lower, xend = upper),
  linewidth = 0.35,
  lineend = "round"
) + geom_text(
  aes(x = or, label = star),
  position = dodge,
  hjust = -0.2,
  vjust = 0.2,
  size  = 2.5,
  show.legend = FALSE
) + geom_xsidetext(
  data = tot,
  aes(x = or, label = star),
  position = dodge,
  hjust = -0.3,
  vjust = 0.2,
  size = 2.5,
  show.legend = FALSE
) + scale_x_log10(
  limits = c(NA, 100),
  oob = squish,
  breaks = c(0.1, 1, 10, 100),
  labels = c("0.1", "1.0", "10", "100")
) + scale_y_discrete(
  expand = expansion(add = 0.5)
) + scale_size(
  range  = c(1, 1.7),
  trans  = "sqrt",
  breaks = c(25, 50, 75, 100),
  guide  = "none"
) + scale_shape_manual(
  breaks = c("Case vs. matched control", "In tumors"),
  values = c("Case vs. matched control" = 15, "In tumors" = 16)
) + scale_color_manual(
  values = c("Case vs. matched control" = "#4d5d7f", "In tumors" = "#eb6d5f")
) + coord_cartesian(
  clip = "off"
) + labs(
  x = "OR per +1 SD",
  y = NULL,
  color = NULL,
  shape = NULL,
  size  = NULL
) + basic + theme(
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.box = "horizontal",
  
  # Plot and legend spacing
  legend.box.spacing = unit(0, "mm"),
  legend.margin = margin(t = 4, r = 0, b = -2, l = 0),
  
  ggside.panel.scale.x = 0.15
) + labs(
  tag = "f"
)

# Save
ggsave("result/or-const.svg", device = svglite, const, width = 1.65, height = 1.6, dpi = 800, bg = "transparent")

