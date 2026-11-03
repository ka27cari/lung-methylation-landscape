# Plot outputs for: 05 longitudinal
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# Section: Paired longitudinal profiles (original 14_longitudinal-pairs.r)
# Figure | Longitudinal prediagnostic TM trend


library(tidyverse)
library(grid)
library(ggside)

plot_font <- "Arial"
panel_tag <- "a"
plot_dpi <- 800
output_dir <- "result"

point_col <- "#93b6c7"
segment_col <- "#000000"
text_col <- "#000000"

point_size <- 1.3
point_alpha <- 0.8
segment_lwd <- 0.3

# 0 근처를 조금 벌리고, 먼 과거는 압축
time_trans <- scales::pseudo_log_trans(sigma = 1.5, base = 10)

# 1. Load

pair_file <- file.path(build_dir, "analysis", "figure_data", "LONGITUDINAL_TERMINAL_MASS_PAIRS.csv")

if (!file.exists(pair_file)) stop("Missing: ", pair_file)

pred_df <- read_csv(pair_file, show_col_types = FALSE) %>%
  transmute(
    matched_index = as.character(matched_index),
    time = as.numeric(time),
    TM = as.numeric(pred_TM)
  ) %>%
  filter(is.finite(time), is.finite(TM), time >= 0, TM >= 0, TM <= 1) %>%
  mutate(
    x = -time,
    time_bin = case_when(
      time > 5 ~ ">5 y",
      time >= 3 & time <= 5 ~ "3-5 y",
      time < 3 ~ "<3 y"
    )
  )

stopifnot(nrow(pred_df) > 0L, !anyDuplicated(pred_df$matched_index))

# 2. Bin summaries

max_time <- max(pred_df$time)

bin_summary <- pred_df %>%
  group_by(time_bin) %>%
  summarise(
    TM_med = median(TM),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    time_bin = factor(
      time_bin,
      levels = c(">5 y", "3-5 y", "<3 y")
    ),
    
    x_start = case_when(
      time_bin == ">5 y"  ~ -max_time,
      time_bin == "3-5 y" ~ -5,
      time_bin == "<3 y"  ~ -3
    ),
    
    x_end = case_when(
      time_bin == ">5 y"  ~ -5,
      time_bin == "3-5 y" ~ -3,
      time_bin == "<3 y"  ~ 0
    ),
    
    label = sprintf("%.2f", TM_med)
  ) %>%
  arrange(time_bin)

print(bin_summary)

# 3. Axis

max_time <- max(pred_df$time)

x_breaks <- c(-20, -15, -10, -5, -3, -2, -1, 0)
x_breaks <- x_breaks[x_breaks >= -max_time & x_breaks <= 0]

# 4. Plot
prob_tail_trans <- scales::new_transform(
  name = "arcsine_sqrt",
  transform = function(x) asin(sqrt(x)),
  inverse   = function(x) sin(x)^2,
  domain = c(0, 1)
)

p_long_tm_pred <- ggplot(pred_df, aes(x = x, y = TM)) +
  geom_point(
    colour = point_col,
    size = point_size,
    alpha = point_alpha,
    shape = 16
  ) +
  geom_ysideboxplot(
    data = pred_df,
    aes(x = "case", y = TM),
    inherit.aes = FALSE,
    orientation = "x",
    width = 0.55,
    fill = "#295A7A",
    colour = segment_col,
    alpha = 0.9,
    linewidth = 0.2,
    outlier.shape = NA
  ) +
  scale_ysidex_discrete(
    name = NULL
  ) +
  geom_vline(
    xintercept = c(-5, -3),
    colour = "#BDBDBD",
    linewidth = 0.18,
    linetype = "22"
  ) +
  geom_segment(
    data = bin_summary,
    aes(
      x = x_start,
      xend = x_end,
      y = TM_med,
      yend = TM_med
    ),
    inherit.aes = FALSE,
    colour = segment_col,
    linewidth = segment_lwd,
    lineend = "butt"
  ) +
  geom_text(
    data = bin_summary,
    aes(
      x = x_end,
      y = TM_med,
      label = label
    ),
    inherit.aes = FALSE,
    family = plot_font,
    size = 2.0,
    colour = text_col,
    hjust = 1.05,
    vjust = -0.55
  ) +
  scale_x_continuous(
    trans = time_trans,
    limits = c(-max_time, 0),
    breaks = x_breaks,
    labels = function(x) format(x, trim = TRUE, scientific = FALSE),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    trans = prob_tail_trans,
    limits = c(0, 1),
    breaks = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1),
    labels = c("0", ".05", ".10", ".25", ".50", ".75", ".90", ".95", "1"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Time-to-diagnosis (years)",
    y = "Terminal prob.",
    tag = panel_tag
  ) +
  basic +
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    ggside.panel.scale = 0.2
  )

print(p_long_tm_pred)

# 5. Save

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "longi-trend.svg"),
  device = svglite,
  p_long_tm_pred,
  width = 2.25,
  height = 1.6,
  units = "in",
  dpi = plot_dpi,
  bg = "transparent"
)

# Figure a2 | Longitudinal paired TM by residual-time group
plot_dpi <- 800
output_dir <- "result"

pred_col <- "#93b6c7"
case_col <- "#295A7A"
line_col <- "#BDBDBD"
median_col <- "#666666"

point_size <- 1.2
point_alpha <- 0.8

pair_lwd <- 0.2
trend_lwd <- 0.3
median_lwd <- 0.55

# 1. Load ----------------------------------------------------------------------

pair_file <- file.path(build_dir, "analysis", "figure_data", "LONGITUDINAL_TERMINAL_MASS_PAIRS.csv")
if (!file.exists(pair_file)) stop("Missing: ", pair_file)

pairs <- read_csv(pair_file, show_col_types = FALSE) %>%
  transmute(matched_index = as.character(matched_index), time = as.numeric(time), pred_TM = as.numeric(pred_TM), case_TM = as.numeric(case_TM)) %>%
  filter(is.finite(time), is.finite(pred_TM), is.finite(case_TM), time >= 0, between(pred_TM, 0, 1), between(case_TM, 0, 1)) %>%
  mutate(
    time_group = case_when(time > 5 ~ ">5 years", time < 3 ~ "<3 years", TRUE ~ NA_character_),
    time_group = factor(time_group, levels = c(">5 years", "<3 years"))
  ) %>%
  filter(!is.na(time_group))

stopifnot(nrow(pairs) > 0L, !anyDuplicated(pairs$matched_index))

# 2. Long format + matched jitter ---------------------------------------------

set.seed(1)

pair_jitter <- pairs %>%
  distinct(matched_index, time_group) %>%
  mutate(jitter = runif(n(), -0.10, 0.10))

x_pred <- 1.05
x_case <- 1.95

plot_df <- pairs %>%
  pivot_longer(c(pred_TM, case_TM), names_to = "sample_type", values_to = "TM") %>%
  mutate(
    sample_type = recode(sample_type, pred_TM = "Prediagnostic", case_TM = "Case"),
    sample_type = factor(sample_type, levels = c("Prediagnostic", "Case"))
  ) %>%
  left_join(pair_jitter, by = c("matched_index", "time_group")) %>%
  mutate(x_base = if_else(sample_type == "Prediagnostic", x_pred, x_case),
         x = x_base + jitter)

pair_lines <- plot_df %>%
  select(matched_index, time_group, sample_type, x, TM) %>%
  pivot_wider(names_from = sample_type, values_from = c(x, TM))

# 3. Medians -------------------------------------------------------------------
median_df <- plot_df %>%
  group_by(time_group, sample_type) %>%
  summarise(
    TM_med = median(TM),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    x = if_else(sample_type == "Prediagnostic", x_pred, x_case),
    xmin = x - 0.16,
    xmax = x + 0.16,
    label = sprintf("%.2f", TM_med),
    
    label_vjust = if_else(
      sample_type == "Prediagnostic",
      1.45,   # 아래
      -0.75   # 위
    )
  )

print(pairs %>% count(time_group))
print(median_df)

# 4. Plot ----------------------------------------------------------------------
stat_df <- pairs %>%
  group_by(time_group) %>%
  group_modify(~ {
    
    test <- wilcox.test(
      .x$case_TM,
      .x$pred_TM,
      paired = TRUE,
      alternative = "greater",
      exact = FALSE
    )
    
    tibble(
      n = nrow(.x),
      pred_median = median(.x$pred_TM),
      case_median = median(.x$case_TM),
      median_change = median(.x$case_TM - .x$pred_TM),
      V = unname(test$statistic),
      p = test$p.value
    )
  }) %>%
  ungroup()

print(stat_df)

bracket_df <- stat_df %>%
  mutate(
    xmin = x_pred,
    xmax = x_case,
    y = 1.075,
    y_tip = 1.010,
    label_y = 1.08,
    label = sprintf("italic(p) == %.3f", p)
  )

trend_df <- median_df %>%
  select(time_group, sample_type, x, TM_med)

p_long_tm_pair <- ggplot() +
  
  # Individual matched pairs
  geom_segment(
    data = pair_lines,
    aes(
      x = x_Prediagnostic,
      xend = x_Case,
      y = TM_Prediagnostic,
      yend = TM_Case,
      group = matched_index
    ),
    colour = line_col,
    linewidth = pair_lwd,
    alpha = 0.55
  ) +
  
  # Prediagnostic points
  geom_point(
    data = plot_df %>% filter(sample_type == "Prediagnostic"),
    aes(x = x, y = TM),
    colour = pred_col,
    size = point_size,
    alpha = point_alpha,
    shape = 16
  ) +
  
  # Case points
  geom_point(
    data = plot_df %>% filter(sample_type == "Case"),
    aes(x = x, y = TM),
    colour = case_col,
    size = point_size,
    alpha = point_alpha,
    shape = 16
  ) +
  
  # Representative trend: median -> median
  geom_line(
    data = trend_df,
    aes(
      x = x,
      y = TM_med,
      group = time_group
    ),
    inherit.aes = FALSE,
    colour = median_col,
    linewidth = trend_lwd,
    lineend = "round"
  ) +
  
  # Median labels
  # geom_text(
  #   data = median_df,
  #   aes(
  #     x = x,
  #     y = TM_med,
  #     label = label,
  #     vjust = label_vjust
  #   ),
  #   inherit.aes = FALSE,
  #   family = plot_font,
  #   size = 2.0
  # ) +
  geom_segment(
    data = bracket_df,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.2
  ) +
  geom_text(
    data = bracket_df,
    aes(
      x = (xmin + xmax) / 2,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    family = plot_font,
    size = 2.0,
    parse = TRUE,
    vjust = 0
  ) +
  facet_wrap(~time_group, nrow = 1) +
  
  scale_x_continuous(
    limits = c(0.65, 2.35),
    breaks = c(x_pred, x_case),
    labels = c("Prediagnostic", "Case"),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.10),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0.075)
  )+
  
  labs(
    x = NULL,
    y = "Terminal probability",
    tag = panel_tag
  ) +
  
  basic +
  
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(10, "pt"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
  )

ggsave(
  file.path(output_dir, "paired-TM.svg"),
  device = svglite,
  p_long_tm_pair,
  width = 2.1, height = 1.6, units = "in", dpi = plot_dpi, bg = "transparent"
)

# Section: Observed route transitions in paired samples (original 15_match-case.r)
# Figure b | Dominant route transition matrix

library(tidyverse)
library(grid)

plot_font <- "Arial"
panel_tag <- "b"
plot_dpi <- 800
output_dir <- "result"
route_levels <- c("T1/T2", "T3", "T4")

# 1. Load ----------------------------------------------------------------------

roles <- read_csv(file.path(build_dir, "analysis", "tables", "FINAL_STATE_ROLE_CLASSIFICATION.csv"), show_col_types = FALSE) %>%
  mutate(state = as.integer(state), role = as.character(macrostate_role))

post <- read_csv(file.path(build_dir, "cache", "longitudinal", "LONGITUDINAL_POSTERIOR.csv"), show_col_types = FALSE) %>%
  transmute(sentrix_id = as.character(sentrix_id), matched_index = as.character(matched_index), group = tolower(as.character(group)))

route <- read_csv(file.path(build_dir, "cache", "longitudinal", "LONGITUDINAL_ROUTE_SCORES.csv"), show_col_types = FALSE) %>%
  mutate(sentrix_id = as.character(sentrix_id))

# 2. Terminal-route column map -------------------------------------------------

terminal_map <- roles %>%
  filter(role == "Terminal") %>%
  arrange(state) %>%
  mutate(terminal = paste0("T", row_number()), route_col = paste0("route_to_T", state))

route_lookup <- setNames(terminal_map$route_col, terminal_map$terminal)

stopifnot(
  all(c("T1", "T2", "T3", "T4") %in% names(route_lookup)),
  all(unname(route_lookup) %in% names(route))
)

# 3. Dominant route per longitudinal sample ------------------------------------

route_df <- route %>%
  transmute(
    sentrix_id,
    T1 = .data[[route_lookup[["T1"]]]],
    T2 = .data[[route_lookup[["T2"]]]],
    T3 = .data[[route_lookup[["T3"]]]],
    T4 = .data[[route_lookup[["T4"]]]]
  ) %>%
  inner_join(post, by = "sentrix_id") %>%
  filter(if_all(c(T1, T2, T3, T4), is.finite)) %>%
  mutate(T12 = T1 + T2, route_sum = T12 + T3 + T4) %>%
  filter(route_sum > 0) %>%
  mutate(
    T12 = T12 / route_sum,
    T3 = T3 / route_sum,
    T4 = T4 / route_sum,
    dominant_route = case_when(
      T12 >= T3 & T12 >= T4 ~ "T1/T2",
      T3 >= T12 & T3 >= T4 ~ "T3",
      TRUE ~ "T4"
    )
  )

# 4. Pair prediagnostic and case -----------------------------------------------

pair_route <- route_df %>%
  select(matched_index, group, dominant_route) %>%
  distinct() %>%
  pivot_wider(names_from = group, values_from = dominant_route) %>%
  filter(!is.na(prediagnostic), !is.na(case))

stopifnot(!anyDuplicated(pair_route$matched_index))

# 5. Transition matrix ----------------------------------------------------------

plot_df <- pair_route %>%
  count(prediagnostic, case, name = "n") %>%
  complete(
    prediagnostic = route_levels,
    case = route_levels,
    fill = list(n = 0)
  ) %>%
  group_by(case) %>%
  mutate(
    case_n = sum(n),
    fraction = if_else(case_n > 0, n / case_n, 0),
    pct = round(100 * fraction)
  ) %>%
  ungroup() %>%
  mutate(
    case = factor(case, levels = rev(route_levels)),
    prediagnostic = factor(prediagnostic, levels = route_levels),
    label = paste0(n, "/", case_n, "\n", pct, "%"),
    concordant = as.character(case) == as.character(prediagnostic)
  )

exact_match <- mean(pair_route$prediagnostic == pair_route$case)

cat(sprintf(
  "Dominant-route concordance: %d/%d = %.1f%%\n",
  sum(pair_route$prediagnostic == pair_route$case),
  nrow(pair_route),
  100 * exact_match
))

print(plot_df)


# 6. Plot ----------------------------------------------------------------------

p_route_transition_heatmap <- ggplot(
  plot_df,
  aes(x = prediagnostic, y = case, fill = fraction)
) +
  geom_tile(
    linewidth = 0.5
  ) +
  geom_tile(
    data = plot_df %>% filter(concordant),
    fill = NA,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(
      label = label,
      colour = fraction >= 0.55
    ),
    family = plot_font,
    size = 1.9,
    lineheight = 0.9
  ) +
  scale_colour_manual(
    values = c(
      `FALSE` = "black",
      `TRUE` = "white"
    ),
    guide = "none"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#EEF3F6",
      "#C8DCE8",
      "#91BAD2",
      "#5796BE",
      "#2F739F"
    ),
    limits = c(0, 1),
    name = "within-case\nfraction",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(13, "mm"),
      barwidth = unit(1.5, "mm"),
      ticks = TRUE,
      frame.colour = NA
    )
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Prediagnostic route",
    y = "Matched case route",
    tag = panel_tag
  ) +
  basic +
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.x.top = element_text(size = 6, margin = margin(b = -1)),
    axis.text.y     = element_text(size = 6, margin = margin(r = -1)),
    axis.title = element_text(size = 7),
    plot.tag.position = c(0.05, 0.95),
    
    legend.title = element_text(
      size = 6,
      hjust = 0.5,
      lineheight = 0.9,
      margin = margin(b = 4)
    ),
    legend.text = element_text(
      size = 5.5
    ),
    legend.margin = margin(l = -4)
  )

p_route_transition_heatmap <- p_route_transition_heatmap +
  theme(legend.position = "none")

print(p_route_transition_heatmap)

# 7. Save ----------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "longitudinal-dominant-route-transition-heatmap.svg"),
  device = svglite,
  p_route_transition_heatmap,
  width = 1.8, height = 1.6, units = "in", dpi = plot_dpi, bg = "transparent"
)

# Figure c | Route-specific precursor validation by terminal branch
# 3 facets: T1/T2, T3, T4
# Each facet: Source | Longitudinal

library(tidyverse)
library(grid)

plot_font <- "Arial"
plot_dpi <- 800
output_dir <- "result"
panel_tag <- "c"
route_levels <- c("T1/T2", "T3", "T4")
source_cut <- 0.60

# 1. Load ----------------------------------------------------------------------

roles <- read_csv(file.path(build_dir, "analysis", "tables", "FINAL_STATE_ROLE_CLASSIFICATION.csv"), show_col_types = FALSE) %>%
  mutate(state = as.integer(state), role = as.character(macrostate_role))

route_aff <- read_csv(file.path(build_dir, "analysis", "tables", "SINKHORN_TERMINAL_ROUTE_AFFINITY.csv"), show_col_types = FALSE) %>%
  mutate(state = as.integer(state), terminal_state = as.integer(terminal_state))

long_post <- read_csv(file.path(build_dir, "cache", "longitudinal", "LONGITUDINAL_POSTERIOR.csv"), show_col_types = FALSE) %>%
  mutate(sentrix_id = as.character(sentrix_id), matched_index = as.character(matched_index), group = tolower(as.character(group)))

long_route <- read_csv(file.path(build_dir, "cache", "longitudinal", "LONGITUDINAL_ROUTE_SCORES.csv"), show_col_types = FALSE) %>%
  mutate(sentrix_id = as.character(sentrix_id))

# 2. State map -----------------------------------------------------------------

state_map <- roles %>%
  filter(n >= 5L | role == "Terminal") %>%
  arrange(factor(role, levels = c("Source", "Intermediate", "Terminal")), state) %>%
  group_by(role) %>%
  mutate(display_index = row_number()) %>%
  ungroup() %>%
  mutate(state_label = case_when(
    role == "Source" ~ paste0("B", display_index),
    role == "Intermediate" ~ paste0("I", display_index),
    role == "Terminal" ~ paste0("T", display_index)
  ))

terminal_map <- state_map %>% filter(role == "Terminal") %>% select(state, terminal = state_label)
bi_map <- state_map %>% filter(role %in% c("Source", "Intermediate")) %>% select(state, state_label, role)

# 3. Source route affinity ------------------------------------------------------

source_affinity <- route_aff %>%
  inner_join(terminal_map, by = c("terminal_state" = "state")) %>%
  mutate(route_branch = case_when(
    terminal %in% c("T1", "T2") ~ "T1/T2",
    terminal == "T3" ~ "T3",
    terminal == "T4" ~ "T4"
  )) %>%
  group_by(state, route_branch) %>%
  summarise(affinity = sum(net95_absorption_affinity, na.rm = TRUE), .groups = "drop") %>%
  complete(state = bi_map$state, route_branch = route_levels, fill = list(affinity = 0)) %>%
  group_by(state) %>%
  mutate(route_sum = sum(affinity), source_affinity = if_else(route_sum > 0, affinity / route_sum, NA_real_)) %>%
  ungroup() %>%
  left_join(bi_map, by = "state") %>%
  filter(route_sum > 0, !is.na(state_label)) %>%
  group_by(state) %>%
  mutate(max_affinity = max(source_affinity, na.rm = TRUE),
         selected = source_affinity >= source_cut & abs(source_affinity - max_affinity) < 1e-12) %>%
  ungroup()

active_states <- source_affinity %>%
  distinct(state, state_label, role) %>%
  arrange(factor(role, levels = c("Source", "Intermediate")), state)

state_order <- active_states$state_label
state_levels_plot <- rev(state_order)
n_intermediate <- sum(active_states$role == "Intermediate")
separator_y <- n_intermediate + 0.5

# 4. Longitudinal matched case dominant route ----------------------------------

case_route <- long_route %>%
  select(sentrix_id, starts_with("route_to_T")) %>%
  pivot_longer(starts_with("route_to_T"), names_to = "route_col", values_to = "route_score") %>%
  mutate(terminal_state = as.integer(str_remove(route_col, "^route_to_T"))) %>%
  inner_join(terminal_map, by = c("terminal_state" = "state")) %>%
  mutate(route_branch = case_when(
    terminal %in% c("T1", "T2") ~ "T1/T2",
    terminal == "T3" ~ "T3",
    terminal == "T4" ~ "T4"
  )) %>%
  group_by(sentrix_id, route_branch) %>%
  summarise(route_score = sum(route_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(sentrix_id) %>%
  mutate(route_sum = sum(route_score), route_score = if_else(route_sum > 0, route_score / route_sum, NA_real_)) %>%
  ungroup() %>%
  inner_join(long_post %>% select(sentrix_id, matched_index, group), by = "sentrix_id") %>%
  filter(group == "case", is.finite(route_score)) %>%
  group_by(matched_index) %>%
  slice_max(route_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(matched_index, case_branch = route_branch)

# 5. Prediagnostic occupancy ---------------------------------------------------

state_cols <- paste0("P_state_", sprintf("%02d", active_states$state))
stopifnot(all(state_cols %in% names(long_post)))

pred_occ <- long_post %>%
  filter(group == "prediagnostic") %>%
  select(matched_index, all_of(state_cols)) %>%
  pivot_longer(all_of(state_cols), names_to = "state_col", values_to = "occupancy") %>%
  mutate(state = as.integer(str_extract(state_col, "\\d+$"))) %>%
  inner_join(active_states, by = "state") %>%
  inner_join(case_route, by = "matched_index")

# 6. Longitudinal one-vs-rest AUC ----------------------------------------------

auc_rank <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  n1 <- length(x); n0 <- length(y)
  r <- rank(c(x, y), ties.method = "average")
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

wilcox_greater <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, alternative = "greater", exact = FALSE)$p.value)
}

long_stats <- pred_occ %>%
  group_by(state, state_label, role) %>%
  group_modify(~ map_dfr(route_levels, function(br) {
    x <- .x$occupancy[.x$case_branch == br]
    y <- .x$occupancy[.x$case_branch != br]
    a <- auc_rank(x, y)
    tibble(route_branch = br, auc = a, effect = a - 0.5, p = wilcox_greater(x, y))
  })) %>%
  ungroup() %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  left_join(source_affinity %>% select(state, state_label, route_branch, source_affinity, selected),
            by = c("state", "state_label", "route_branch")) %>%
  filter(!is.na(state_label)) %>%
  mutate(sig = is.finite(q) & q < 0.05,
         validated = selected & sig & effect > 0)

# 7. Plot data -----------------------------------------------------------------

cell_df <- expand_grid(
  route_branch = route_levels,
  panel = c("Source", "Longitudinal"),
  state_label = state_levels_plot
) %>%
  mutate(
    route_branch = factor(route_branch, levels = route_levels),
    panel = factor(panel, levels = c("Source", "Longitudinal")),
    x_num = if_else(panel == "Source", 1, 2),
    y_num = match(state_label, state_levels_plot),
    cell_xmin = x_num - 0.46,
    cell_xmax = x_num + 0.46,
    cell_ymin = y_num - 0.46,
    cell_ymax = y_num + 0.46
  )

long_no_source_df <- cell_df %>%
  filter(panel == "Longitudinal") %>%
  left_join(
    source_affinity %>% select(route_branch, state_label, source_affinity),
    by = c("route_branch", "state_label")
  ) %>%
  filter(is.na(source_affinity) | source_affinity == 0)

source_plot_df <- source_affinity %>%
  filter(!is.na(state_label)) %>%
  mutate(
    route_branch = factor(route_branch, levels = route_levels),
    x_num = 1,
    y_num = match(state_label, state_levels_plot),
    cell_xmin = x_num - 0.46,
    cell_xmax = x_num + 0.46,
    cell_ymin = y_num - 0.46,
    cell_ymax = y_num + 0.46,
    bar_xmin = cell_xmin,
    bar_xmax = cell_xmin + 0.84 * source_affinity,
    bar_ymin = y_num - 0.18,
    bar_ymax = y_num + 0.18
  )

long_plot_df <- long_stats %>%
  filter(!is.na(state_label)) %>%
  mutate(
    route_branch = factor(route_branch, levels = route_levels),
    x_num = 2,
    y_num = match(state_label, state_levels_plot),
    cell_xmin = x_num - 0.46,
    cell_xmax = x_num + 0.46,
    cell_ymin = y_num - 0.46,
    cell_ymax = y_num + 0.46,
    effect_clip = pmax(pmin(effect, 0.5), -0.5),
    dx = 0.38 * effect_clip / 0.5,
    bar_xmin = if_else(dx < 0, x_num + dx, x_num),
    bar_xmax = if_else(dx > 0, x_num + dx, x_num),
    bar_ymin = y_num - 0.18,
    bar_ymax = y_num + 0.18
  )

# 8. Plot ----------------------------------------------------------------------

p_route_precursor_validation <- ggplot() +
  geom_rect(
    data = cell_df,
    aes(xmin = cell_xmin, xmax = cell_xmax, ymin = cell_ymin, ymax = cell_ymax),
    fill = "#FBFBFB", colour = "#D8D8D8", linewidth = 0.15
  ) +
  
  geom_rect(
    data = source_plot_df %>% filter(selected),
    aes(xmin = cell_xmin, xmax = cell_xmax, ymin = cell_ymin, ymax = cell_ymax),
    fill = "#F3F3F3", colour = NA
  ) +
  
  geom_rect(
    data = source_plot_df,
    aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bar_ymin, ymax = bar_ymax, fill = source_affinity),
    colour = NA
  ) +
  
  scale_fill_gradient(
    low = "#E7EFF5", high = "#3F6E99",
    limits = c(0, 1), name = "Route\naffinity",
    guide = guide_colourbar(
      title.position = "top", title.hjust = 0,
      barheight = unit(9, "mm"), barwidth = unit(1.3, "mm")
    )
  ) +
  
  geom_rect(
    data = source_plot_df %>% filter(selected),
    aes(xmin = cell_xmin, xmax = cell_xmax, ymin = cell_ymin, ymax = cell_ymax),
    fill = NA, colour = "black", linewidth = 0.15
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # Source affinity가 있는 경우에만 significant 배경 강조
  geom_rect(
    data = long_plot_df %>% filter(sig, selected),
    aes(xmin = cell_xmin, xmax = cell_xmax, ymin = cell_ymin, ymax = cell_ymax),
    fill = "#F3F3F3", colour = NA
  ) +
  
  geom_rect(
    data = long_plot_df,
    aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bar_ymin, ymax = bar_ymax, fill = auc),
    colour = NA
  ) +
  
  geom_rect(
    data = long_plot_df %>% filter(sig),
    aes(xmin = cell_xmin, xmax = cell_xmax, ymin = cell_ymin, ymax = cell_ymax),
    fill = NA, colour = "black", linewidth = 0.15
  ) +
  geom_segment(
    data = cell_df %>% filter(panel == "Longitudinal"),
    aes(x = x_num, xend = x_num, y = y_num - 0.28, yend = y_num + 0.28),
    colour = "black", linewidth = 0.15
  ) +
  facet_wrap(~route_branch, nrow = 1) +
  
  scale_fill_gradient2(
    low = "#5B80B2", mid = "#F7F7F7", high = "#D55C5C",
    midpoint = 0.5, limits = c(0.1, 0.9), oob = scales::squish,
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
    name = "AUC",
    guide = guide_colourbar(
      title.position = "top", title.hjust = 0,
      barheight = unit(9, "mm"), barwidth = unit(1.3, "mm")
    )
  ) +
  
  scale_x_continuous(
    limits = c(0.5, 2.5),
    breaks = c(1, 2),
    labels = c("Nested", "Longit."),
    position = "bottom",
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0.5, length(state_levels_plot) + 0.5),
    breaks = seq_along(state_levels_plot),
    labels = state_levels_plot,
    expand = c(0, 0)
  ) +
  
  labs(x = NULL, y = NULL, tag = panel_tag) +
  basic +
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6, margin = margin(t = 2)),
    axis.text.y = element_text(size = 6),
    strip.background = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.box.spacing = unit(0.5, "mm"),
    legend.margin = margin(l = 3.5),
    legend.title = element_text(size = 6, margin = margin(b = 3)),
    legend.text = element_text(size = 5.2, margin = margin(l = 2))
  )

print(p_route_precursor_validation)

# 9. Save ----------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "route-precursor-validation.svg"),
  device = svglite,
  p_route_precursor_validation,
  width = 3, height = 1.6, units = "in", dpi = plot_dpi, bg = "transparent"
)

# Section: Program and blood-cell context by route (original 16_src-term-prog-cell.r)
# Figure d | Terminal blood-composition + program displacement dot map
# Single-panel version


library(tidyverse)
library(grid)
library(ggtext)

plot_font <- "Arial"
plot_dpi <- 800
output_dir <- "result"
panel_tag <- "d"

# 1. Settings -----------------------------------------------------------------

cell_types <- c("Bas","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Eos","Mono","Neu","NK","Treg")

cell_labels <- c(
  Neu = "Neutrophil", Mono = "Monocyte", Eos = "Eosinophil", Bas = "Basophil",
  Bnv = "B naive", Bmem = "B memory", CD4nv = "CD4 naive", CD4mem = "CD4 memory",
  CD8nv = "CD8 naive", CD8mem = "CD8 memory", NK = "NK cell", Treg = "Treg"
)

cell_order <- c(
  "Neu","Mono","Eos","Bas",
  "Bnv","Bmem","CD4nv","CD4mem","CD8nv","CD8mem","NK","Treg"
)

program_order <- c("P1","P2","P6","P3","P9","P7","P5","P8","P4")
route_levels <- c("T1","T2","T3","T4","Pooled")

# Actual physical spacing in one coordinate system
x_step <- 0.72
myeloid_lymphoid_gap <- 0.3
cell_program_gap <- 1.1

y_step <- 0.78
pooled_gap <- 0.14

dot_size <- 2.5
outline_size <- 2.73
outline_stroke <- 0.23

low_col <- "#2c3173"
mid_col <- "#FAFAFA"
high_col <- "#de1c10"

program_family_cols <- c(
  ReprPC = "#292473",
  Bivalent = "#C2185B"
)

# 2. Display coordinates -------------------------------------------------------

cell_map <- tibble(
  cell_type = cell_order,
  display_order = seq_along(cell_order)
) %>%
  mutate(
    lymphoid = !(cell_type %in% c("Neu","Mono","Eos","Bas")),
    plot_x = (display_order - 1) * x_step +
      myeloid_lymphoid_gap * as.numeric(lymphoid)
  )

program_start <- max(cell_map$plot_x) + cell_program_gap

program_map <- tibble(
  program = program_order,
  plot_x = program_start + (seq_along(program_order) - 1) * x_step
)

route_map <- tibble(
  branch = route_levels,
  plot_y = c(
    4 * y_step + pooled_gap,
    3 * y_step + pooled_gap,
    2 * y_step + pooled_gap,
    1 * y_step + pooled_gap,
    0
  )
)

program_tick_labels <- setNames(program_order, program_order)

program_tick_labels[c("P3","P9")] <- paste0(
  "<span style='color:", program_family_cols["ReprPC"], "'>",
  c("P3","P9"), "</span>"
)

program_tick_labels[c("P5","P8")] <- paste0(
  "<span style='color:", program_family_cols["Bivalent"], "'>",
  c("P5","P8"), "</span>"
)

# 3. Load ---------------------------------------------------------------------

meta <- read_csv(
  file.path(build_dir, "data", "nested.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    sentrix_id = as.character(sentrix_id),
    across(all_of(cell_types), as.numeric)
  )

post <- read_csv(
  file.path(build_dir, "cache", "model", "SOURCE_POSTERIOR.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    sentrix_id = as.character(sentrix_id),
    group = tolower(as.character(group)),
    hard_state = as.integer(hard_state)
  )

roles <- read_csv(
  file.path(build_dir, "analysis", "tables", "FINAL_STATE_ROLE_CLASSIFICATION.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    state = as.integer(state),
    role = as.character(macrostate_role),
    n = as.numeric(n)
  )

centroids <- read_csv(
  file.path(build_dir, "cache", "model", "STATE_CENTROIDS.csv"),
  show_col_types = FALSE
) %>%
  mutate(state = as.integer(state))

D <- post %>%
  inner_join(meta, by = "sentrix_id")

# 4. Blood composition: CLR + robust scaling ----------------------------------

A <- as.matrix(D[, cell_types])
storage.mode(A) <- "double"

col_med <- apply(A, 2, median, na.rm = TRUE)

for (j in seq_along(cell_types)) {
  bad <- !is.finite(A[, j])
  A[bad, j] <- col_med[j]
}

positive_values <- A[is.finite(A) & A > 0]

if (!length(positive_values)) {
  stop("No positive blood-composition values available.")
}

pseudo_count <- min(positive_values) / 2
A[A <= 0] <- pseudo_count

B <- log(A)
B <- B - rowMeans(B)

pred_mask <- D$group == "prediagnostic"

robust_center <- apply(B[pred_mask, , drop = FALSE], 2, median)
q10 <- apply(B[pred_mask, , drop = FALSE], 2, quantile, probs = 0.10)
q90 <- apply(B[pred_mask, , drop = FALSE], 2, quantile, probs = 0.90)

robust_scale <- q90 - q10
robust_scale[!is.finite(robust_scale) | robust_scale < 1e-12] <- 1

Bs <- sweep(B, 2, robust_center, "-")
Bs <- sweep(Bs, 2, robust_scale, "/")

# 5. State definitions ---------------------------------------------------------

source_states <- roles %>%
  filter(role == "Source") %>%
  pull(state)

terminal_states <- roles %>%
  filter(role == "Terminal") %>%
  arrange(state) %>%
  pull(state)

if (length(terminal_states) != 4L) {
  stop("Expected four terminal states; found ", length(terminal_states))
}

terminal_map <- tibble(
  state = terminal_states,
  branch = paste0("T", seq_along(terminal_states))
)

# 6. Blood-composition displacement -------------------------------------------

base_mask <- pred_mask & D$hard_state %in% source_states

if (!any(base_mask)) {
  stop("No prediagnostic Source-state samples.")
}

base_vec <- apply(
  Bs[base_mask, , drop = FALSE],
  2,
  median
)

get_cell_shift <- function(mask, label) {
  if (sum(mask) == 0L) stop("No samples for ", label)
  
  term_vec <- apply(
    Bs[mask, , drop = FALSE],
    2,
    median
  )
  
  tibble(
    branch = label,
    feature = cell_types,
    value = as.numeric(term_vec - base_vec),
    feature_type = "Cell"
  )
}

cell_shift <- bind_rows(
  get_cell_shift(pred_mask & D$hard_state %in% terminal_states, "Pooled"),
  get_cell_shift(pred_mask & D$hard_state == terminal_states[1], "T1"),
  get_cell_shift(pred_mask & D$hard_state == terminal_states[2], "T2"),
  get_cell_shift(pred_mask & D$hard_state == terminal_states[3], "T3"),
  get_cell_shift(pred_mask & D$hard_state == terminal_states[4], "T4")
) %>%
  left_join(
    cell_map %>% select(cell_type, plot_x),
    by = c("feature" = "cell_type")
  )

# 7. Program displacement ------------------------------------------------------

program_long <- centroids %>%
  select(state, starts_with("program_")) %>%
  pivot_longer(
    cols = starts_with("program_"),
    names_to = "program_raw",
    values_to = "score"
  ) %>%
  mutate(
    program_index = as.integer(str_extract(program_raw, "[0-9]+")),
    program = paste0("P", program_index)
  ) %>%
  filter(program %in% program_order) %>%
  left_join(
    roles %>% select(state, role, n),
    by = "state"
  )

stopifnot(
  all(is.finite(program_long$score)),
  all(is.finite(program_long$n)),
  all(program_long$n > 0)
)

source_program <- program_long %>%
  filter(role == "Source") %>%
  group_by(program) %>%
  summarise(
    base = weighted.mean(score, n),
    .groups = "drop"
  )

terminal_program <- program_long %>%
  filter(state %in% terminal_states) %>%
  left_join(
    terminal_map,
    by = "state"
  ) %>%
  left_join(
    source_program,
    by = "program"
  )

program_terminal_shift <- terminal_program %>%
  transmute(
    branch,
    feature = program,
    value = score - base,
    feature_type = "Program"
  )

program_pooled_shift <- terminal_program %>%
  group_by(program) %>%
  summarise(
    term = weighted.mean(score, n),
    base = first(base),
    .groups = "drop"
  ) %>%
  transmute(
    branch = "Pooled",
    feature = program,
    value = term - base,
    feature_type = "Program"
  )

program_shift <- bind_rows(
  program_pooled_shift,
  program_terminal_shift
) %>%
  left_join(
    program_map,
    by = c("feature" = "program")
  )

# 8. Combine ------------------------------------------------------------------

plot_df <- bind_rows(
  cell_shift,
  program_shift
) %>%
  left_join(
    route_map,
    by = "branch"
  )

stopifnot(
  nrow(plot_df) ==
    length(route_levels) * (length(cell_order) + length(program_order)),
  !anyNA(plot_df$plot_x),
  !anyNA(plot_df$plot_y),
  all(is.finite(plot_df$value))
)

clip_lim <- quantile(abs(plot_df$value), 0.95, na.rm = TRUE)
clip_lim <- max(clip_lim, 0.25)

plot_df <- plot_df %>%
  mutate(
    value_clip = pmax(
      pmin(value, clip_lim),
      -clip_lim
    )
  )

# 9. Data-driven outlines ------------------------------------------------------

terminal_n <- D %>%
  filter(
    group == "prediagnostic",
    hard_state %in% terminal_states
  ) %>%
  count(
    hard_state,
    name = "n_branch"
  ) %>%
  left_join(
    terminal_map,
    by = c("hard_state" = "state")
  )

small_n_cut <- median(terminal_n$n_branch)

outline_pairs <- plot_df %>%
  filter(branch != "Pooled") %>%
  left_join(
    terminal_n %>% select(branch, n_branch),
    by = "branch"
  ) %>%
  group_by(feature_type, feature) %>%
  mutate(
    other_mean = (sum(value) - value) / (n() - 1),
    distinct_score = abs(value - other_mean)
  ) %>%
  arrange(
    desc(distinct_score),
    n_branch,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  filter(n_branch <= small_n_cut) %>%
  select(
    branch,
    feature_type,
    feature
  )

plot_df <- plot_df %>%
  left_join(
    outline_pairs %>% mutate(outline = TRUE),
    by = c("branch","feature_type","feature")
  ) %>%
  mutate(
    outline = replace_na(outline, FALSE)
  )

# 10. Axis labels --------------------------------------------------------------

x_axis <- bind_rows(
  cell_map %>%
    transmute(
      plot_x,
      label = unname(cell_labels[cell_type])
    ),
  
  program_map %>%
    transmute(
      plot_x,
      label = unname(program_tick_labels[program])
    )
) %>%
  arrange(plot_x)

# Extra space around outer points, without changing their internal spacing
x_min <- min(x_axis$plot_x) - 0.35
x_max <- max(x_axis$plot_x) + 0.35

y_axis <- min(route_map$plot_y) - 0.38
y_max <- max(route_map$plot_y) + 0.35

# Separate x-axis line segments
cell_axis_start <- min(cell_map$plot_x) - 0.20
cell_axis_end <- max(cell_map$plot_x) + 0.20

program_axis_start <- min(program_map$plot_x) - 0.20
program_axis_end <- max(program_map$plot_x) + 0.20

# 11. Plot --------------------------------------------------------------------

p_dotmap <- ggplot(
  plot_df,
  aes(x = plot_x, y = plot_y)
) +
  
  geom_point(
    aes(colour = value_clip),
    shape = 16,
    size = dot_size
  ) +
  
  geom_point(
    data = plot_df %>% filter(outline),
    shape = 1,
    size = outline_size,
    colour = "black",
    stroke = outline_stroke,
    show.legend = FALSE
  ) +
  
  scale_colour_gradient2(
    low = low_col,
    mid = mid_col,
    high = high_col,
    midpoint = 0,
    limits = c(-clip_lim, clip_lim),
    oob = scales::squish,
    name = "difference\n(terminal-base)",
    guide = guide_colorbar(
      direction = "vertical",
      barheight = unit(13, "mm"),
      barwidth = unit(1.3, "mm"),
      ticks = TRUE,
      label.position = "right",
      title.position = "right",
      title.hjust = 0.5,
      theme = theme(
        legend.text = element_text(
          family = plot_font,
          size = 5.4,
          margin = margin(l = 2)
        ),
        
        legend.title = element_text(
          family = plot_font,
          size = 6,
          angle = 270,
          hjust = 0.5,
          vjust = 0.5,
          margin = margin(l = 2)
        )
      )
    )
  ) +
  
  scale_x_continuous(
    breaks = x_axis$plot_x,
    labels = x_axis$label,
    limits = c(x_min, x_max),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    breaks = route_map$plot_y,
    labels = route_map$branch,
    limits = c(y_axis, y_max),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    tag = panel_tag
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  basic +
  
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    
    axis.text.x = ggtext::element_markdown(
      family = plot_font,
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      margin = margin(t = 2)
    ),
    
    axis.text.y = element_text(
      family = plot_font,
      size = 6,
      margin = margin(r = 2)
    ),
    
    legend.position = "right",
    legend.box.spacing = unit(0, "mm"),
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 5,
      unit = "pt"
    )
  )

print(p_dotmap)

# 12. Save --------------------------------------------------------------------

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  file.path(
    output_dir,
    "blood-program-route-dotmap.svg"
  ),
  device = svglite,
  p_dotmap,
  width = 4,
  height = 1.55,
  units = "in",
  dpi = plot_dpi,
  bg = "transparent"
)


# Section: Hematopoietic displacement (original 17_hemotopoietic-displace.r)
# Figure e | Pooled vs terminal-aware longitudinal coupling

library(tidyverse)
library(ggtext)
library(grid)

plot_font <- "Arial"
panel_tag <- "e"
plot_dpi <- 800
output_dir <- "result"

point_col <- "#7A8794"
line_col <- "#914453"
n_boot <- 5000
n_perm <- 9999

# 1. Load ---------------------------------------------------------------------

coupling <- read_csv(
  file.path(build_dir, "analysis", "figure_data", "LONGITUDINAL_TERMINAL_AWARE_COUPLING.csv"),
  show_col_types = FALSE
)

stopifnot(all(c(
  "matched_index",
  "pooled_fate_displacement", "pooled_hemo_displacement",
  "terminal_aware_fate_displacement", "terminal_aware_hemo_displacement"
) %in% names(coupling)))

# 2. Long format ---------------------------------------------------------------

plot_df <- bind_rows(
  coupling %>% transmute(
    matched_index,
    model = "Pooled",
    fate = pooled_fate_displacement,
    hemo = pooled_hemo_displacement
  ),
  coupling %>% transmute(
    matched_index,
    model = "Terminal-conditioned",
    fate = terminal_aware_fate_displacement,
    hemo = terminal_aware_hemo_displacement
  )
) %>%
  filter(is.finite(fate), is.finite(hemo)) %>%
  mutate(model = factor(model, levels = c("Pooled", "Terminal-conditioned")))

# 3. Bootstrap CI + permutation p ---------------------------------------------

calc_stats <- function(d, n_boot = 5000, n_perm = 9999, seed = 20261005) {
  x <- d$fate
  y <- d$hemo
  n <- length(x)
  
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  
  set.seed(seed)
  boot_rho <- replicate(n_boot, {
    ii <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[ii], y[ii], method = "spearman"))
  })
  
  ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)
  
  set.seed(seed + 1)
  perm_rho <- replicate(n_perm, {
    suppressWarnings(cor(x, sample(y), method = "spearman"))
  })
  
  perm_p <- (1 + sum(perm_rho >= rho, na.rm = TRUE)) / (n_perm + 1)
  
  tibble(
    rho = rho,
    ci_lo = unname(ci[1]),
    ci_hi = unname(ci[2]),
    perm_p = perm_p,
    n = n
  )
}

stat_df <- plot_df %>%
  group_by(model) %>%
  group_modify(~ calc_stats(.x, n_boot = n_boot, n_perm = n_perm)) %>%
  ungroup()

print(stat_df)

# 4. Statistical labels --------------------------------------------------------

label_pos <- plot_df %>%
  group_by(model) %>%
  summarise(
    x = min(fate) + 0.03 * diff(range(fate)),
    y = max(hemo) - 0.03 * diff(range(hemo)),
    .groups = "drop"
  )

stat_df <- stat_df %>%
  mutate(
    p_text = if_else(
      perm_p < 0.001,
      "<i>p</i> &lt; 0.001",
      paste0("<i>p</i> = ", sprintf("%.3f", perm_p))
    ),
    label = paste0(
      "\u03c1 = ", sprintf("%.2f", rho), " [", sprintf("%.2f", ci_lo), ", ", sprintf("%.2f", ci_hi), "]",
      # "<br>95% bootstrap CI [", sprintf("%.2f", ci_lo), ", ", sprintf("%.2f", ci_hi), "]",
      "<br>permutation ", p_text
    )
  ) %>%
  left_join(label_pos, by = "model")

# 5. Plot ---------------------------------------------------------------------

p_coupling <- ggplot(plot_df, aes(x = fate, y = hemo)) +
  geom_point(
    colour = point_col,
    size = 0.8,
    alpha = 0.55,
    shape = 16
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = line_col,
    linewidth = 0.35,
    lineend = "round"
  ) +
  ggtext::geom_richtext(
    data = stat_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    family = plot_font,
    size = 1.9,
    hjust = 0,
    vjust = 1,
    lineheight = 1.0,
    fill = NA,
    label.color = NA,
    label.padding = unit(0, "pt")
  ) +
  facet_wrap(
    ~model,
    nrow = 1,
    scales = "free"
  ) +
  labs(
    x = "Program displacement",
    y = "Cell displacement",
    tag = panel_tag
  ) +
  basic +
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    panel.spacing.x = unit(3, "mm"),
    axis.text = element_text(size = 6.6),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(
      family = plot_font,
      size = 7,
      face = "plain",
      margin = margin(b = 2)
    ),
    plot.margin = margin(l = 12, r = 5, b = 5, t = 5),
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 16,
      margin = margin(l = -10)
    ),
    plot.tag.location = "plot",
    plot.tag.position = c(0, 0.95)
  )

print(p_coupling)

# 6. Save ---------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "coupling-pooled-terminal-aware.svg"),
  device = svglite,
  p_coupling,
  width = 2.95,
  height = 1.55,
  units = "in",
  dpi = plot_dpi,
  bg = "transparent"
)
