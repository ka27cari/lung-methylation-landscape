# Plot outputs for: 04 associations
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# Section: Terminal-conditioned methylation eigengenes (original 11_me-bayes.r)
# Figure 4h | Terminal-conditioned Bayesian ME proximity effects
# - no overall title / subtitle / caption
# - significant MEs are shown in bold and marked with †
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)
library(grid)

plot_font <- "Arial"
panel_tag <- "i"
sig_mark <- "\u2020"

me <- read_csv(
  file.path(build_dir, "analysis", "tables", "TERMINAL_CONDITIONED_ME_HORIZON.csv"),
  show_col_types = FALSE
) %>%
  filter(
    config == "C_ASSA",
    predictor == "proximity",
    condition %in% c("Pooled", "T1", "T2", "T3", "T4"),
    feature %in% paste0("ME", 0:6)
  ) %>%
  mutate(
    condition = factor(condition, levels = c("Pooled", "T1", "T2", "T3", "T4")),
    feature = factor(feature, levels = rev(paste0("ME", 0:6))),
    interval_supported = (`CrI_2.5` > 0) | (`CrI_97.5` < 0),
    supported = as.logical(supported)
  )

# stopifnot(nrow(me) == 28L)
stopifnot(all(me$supported == me$interval_supported))

condition_cols <- c(
  Pooled = "#3F3F3F",
  T1 = "#5E49B8",
  T2 = "#008C8C",
  T3 = "#CF9000"
)

x_lim <- max(abs(c(me$`CrI_2.5`, me$`CrI_97.5`)), na.rm = TRUE)
x_lim <- ceiling(x_lim * 10) / 10
label_x <- -x_lim * 1.10
left_lim <- label_x - x_lim * 0.05

plot_df <- me %>%
  mutate(
    feature_chr = as.character(feature),
    feature_display = if_else(supported, paste0(feature_chr, " ", sig_mark), feature_chr),
    feature_face = if_else(supported, "bold", "plain")
  )

p_terminal_conditioned_me <- ggplot(plot_df, aes(y = feature)) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.15,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = `CrI_2.5`, xend = `CrI_97.5`, yend = feature, colour = condition),
    linewidth = 0.60,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = posterior_median, colour = condition),
    size = 1.8,
    shape = 16,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, label = feature_display, fontface = feature_face),
    hjust = 1,
    family = plot_font,
    size = 2.2,
    colour = "#222222",
    show.legend = FALSE
  ) +
  facet_wrap(~ condition, nrow = 1) +
  scale_colour_manual(values = condition_cols) +
  scale_x_continuous(
    limits = c(left_lim, x_lim),
    breaks = scales::breaks_pretty(n = 4)
  ) +
  labs(tag = panel_tag, x = NULL, y = NULL) +
  basic +
  coord_cartesian(clip = "off") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.title = element_blank(),
    plot.tag.position = c(-0.01, 0.95),
    # panel.spacing.x = unit(5, "mm")
  )

print(p_terminal_conditioned_me)

dir.create("result", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "result/me-bayes.svg",
  p_terminal_conditioned_me,
  device = svglite,
  width = 5.3,
  height = 1.4,
  units = "in",
  bg = "transparent"
)


# Section: Out-of-fold program discrimination (original 12_prog-auc.r)
# Conditional OOF AUC comparison
# Original horizontal bar format
#
# DATA ARE READ DIRECTLY FROM:
# INCREMENTAL_INFORMATION_OOF_SUMMARY.csv
#
# Bar = mean AUC above chance (0.5)
# Error bar = SD across repeated OOF runs
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)
library(grid)

plot_font <- "Arial"

# 1. Locate REAL result table --------------------------------------------------

auc_filename <- "INCREMENTAL_INFORMATION_OOF_SUMMARY.csv"

auc_candidates <- c(
  file.path(build_dir, "analysis", "tables", auc_filename),
  file.path(dirname(build_dir), "build", "analysis", "tables", auc_filename)
)

auc_file <- auc_candidates[file.exists(auc_candidates)][1]

if (length(auc_file) == 0 || is.na(auc_file)) {
  stop(
    "Cannot find ", auc_filename,
    "\nChecked:\n",
    paste(auc_candidates, collapse = "\n")
  )
}

cat("AUC table:", auc_file, "\n")

auc_all <- read_csv(auc_file, show_col_types = FALSE)

stopifnot(all(c("model", "AUC_mean", "AUC_sd") %in% names(auc_all)))

# 2. Main models to display ----------------------------------------------------

model_order <- c(
  "Epidemiology",
  "Epidemiology+Blood",
  "Epidemiology+Blood+ME",
  "Epidemiology+Blood+ME+TM",
  "Epidemiology+Blood+ME+Programs",
  "Programs",
  "TM"
)

model_labels <- c(
  "Epidemiology" = "Epidemiology",
  "Epidemiology+Blood" = "Epidemiology + Blood",
  "Epidemiology+Blood+ME" = "Epidemiology + Blood + ME",
  "Epidemiology+Blood+ME+TM" = "Epidemiology + Blood + ME + TM",
  "Epidemiology+Blood+ME+Programs" = "Epidemiology + Blood + ME + Programs",
  "Programs" = "Programs only",
  "TM" = "TM only"
)

auc_df <- auc_all %>%
  filter(model %in% model_order) %>%
  mutate(
    label = unname(model_labels[model]),
    auc_above_chance = AUC_mean - 0.5,
    error_low = AUC_mean - AUC_sd - 0.5,
    error_high = AUC_mean + AUC_sd - 0.5,
    feature_group = case_when(
      str_detect(model, "Programs") ~ "Programs",
      str_detect(model, "TM") ~ "TM",
      TRUE ~ "Reference"
    )
  ) %>%
  mutate(label = fct_reorder(label, AUC_mean))

print(
  auc_df %>%
    arrange(desc(AUC_mean)) %>%
    select(model, AUC_mean, AUC_sd)
)

# 3. Colours -------------------------------------------------------------------

bar_cols <- c(
  Reference = "#BFC3C7",
  TM = "#527A7A",
  Programs = "#7568A0"
)

# 4. Plot ----------------------------------------------------------------------

p_auc <- ggplot(
  auc_df,
  aes(x = label, y = auc_above_chance, fill = feature_group)
) +
  geom_errorbar(
    aes(ymin = error_low, ymax = error_high),
    width = 0.17,
    linewidth = 0.2,
    colour = "#333333"
  ) + geom_col(
    width = 0.68,
    colour = NA
  ) +
  geom_text(
    aes(label = sprintf("%.3f", AUC_mean)),
    hjust = -0.3,
    family = plot_font,
    size = 1.85,
    # colour = "#222222"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = bar_cols,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 0.35),
    breaks = seq(0, 0.35, by = 0.05),
    labels = function(x) sprintf("%.2f", x + 0.5),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Conditional OOF AUC",
    tag = "h"
  ) +
  basic +
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    plot.tag.position  = c(-0.02, 0.95),
    plot.margin = margin(l = 15, r = 5, b = 5, t = 5),
  )

ggsave(
  "result/prog-aucs.svg",
  device = svglite,
  p_auc,
  width = 3.2,
  height = 1.3,
  units = "in",
  bg = "transparent"
)

