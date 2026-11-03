# Plot outputs for: 06 external
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# External PBMC validation and route-level comparison.

library(tidyverse)

read_build <- function(...) read.csv(file.path(build_dir, ...), check.names = FALSE)

# Load the frozen, source-control-aligned projection produced by the pipeline.
prog <- sprintf("program_%02d", 1:9)
terminal <- c(T1 = 2L, T2 = 6L, T3 = 7L, T4 = 8L)
S <- read_build("cache", "model", "SOURCE_PROGRAM_SCORES.csv")
P <- read_build("cache", "model", "SOURCE_POSTERIOR.csv")
R <- read_build("cache", "model", "SOURCE_ROUTE_SCORES.csv")
E <- read_build("cache", "external", "EXTERNAL_PROGRAM_SCORES.csv")
Epost <- read_build("cache", "external", "EXTERNAL_POSTERIOR.csv")
A <- read_build("analysis", "tables", "STATE_TO_TERMINAL_ABSORPTION.csv")

P <- P[match(S$sentrix_id, P$sentrix_id), ]
R <- R[match(S$sentrix_id, R$sentrix_id), ]
Epost <- Epost[match(E$sample, Epost$sample), ]
A <- A[order(A$state), paste0("T", unname(terminal))]
Xs <- as.matrix(S[, prog]); Xa <- as.matrix(E[, prog])
source_control <- P$group == "control"
external_control <- E$group == "Control"
Q <- as.matrix(Epost[, sprintf("P_state_%02d", 0:15)])
E$TM <- Epost$TM
E$route <- names(terminal)[max.col(Q %*% as.matrix(A))]
S$route <- names(terminal)[max.col(as.matrix(R[, paste0("route_to_T", unname(terminal))]))]
source_route <- S$route
ext_control <- external_control
auroc <- read_build("analysis", "tables", "EXTERNAL_VALIDATION_SUMMARY.csv")$AUC[[1]]

plot_font <- "Arial"
E$group <- factor(
  E$group,
  levels = c("Control", "LUAD")
)
stat_lab <- sprintf("AUC = %.2f", auroc)

gse_id <- "GSE275371"

p <- ggplot(E, aes(x = group, y = TM)) +
  geom_violin(
    aes(fill = group),
    width = 0.3,
    scale = "width",
    trim = TRUE,
    colour = NA,
    alpha = 0.45
  ) +
  geom_point(
    aes(colour = group),
    position = position_jitter(
      width = 0.10,
      height = 0,
      seed = 1
    ),
    size = 0.85,
    alpha = 0.75,
    shape = 16
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.20,
    linewidth = 0.35,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.10,
    label = stat_lab,
    family = plot_font,
    size = 1.9
  ) +
  scale_fill_manual(
    values = c(
      Control = "#9A9A9A",
      LUAD = "#db6d5f"
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      Control = "#666666",
      LUAD = "#db6d5f"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1.15),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = gse_id,
    y = "Terminal mass",
    tag = "f"
  ) +
  basic +
  theme(
    text = element_text(family = plot_font),
    axis.text = element_text(size = 6.5),
    axis.title.y = element_text(size = 7)
  )

# 4. Save ----------------------------------------------------------------------

ggsave(
  file.path("result/external-terminal-mass.svg"),
  device = svglite,
  p,
  width = 1.7,
  height = 1.6,
  units = "in",
  bg = "transparent"
)

families <- list(
  "ReprPC" = c("P3", "P9"),
  "Bivalent" = c("P5", "P8")
)

family_colors <- c(
  "ReprPC" = "#292473",
  "Bivalent" = "#C2185B"
)

# Source/external control-relative directions, faceted by preferred route.
source_delta <- colMeans(Xs[source_control, , drop = FALSE])
external_delta <- colMeans(Xa[external_control, , drop = FALSE])
directions <- bind_rows(lapply(names(terminal), function(route) {
  dx <- colMeans(Xs[P$group == "prediagnostic" & S$route == route, , drop = FALSE]) - source_delta
  dy <- colMeans(Xa[E$group == "LUAD" & E$route == route, , drop = FALSE]) - external_delta
  data.frame(route = route, program = paste0("P", 1:9),
             source = dx, external = dy, agreement = sign(dx) == sign(dy))
}))

directions <- directions %>%
  mutate(
    family = case_when(
      program %in% families$ReprPC  ~ "ReprPC",
      program %in% families$Bivalent ~ "Bivalent",
      TRUE ~ "Other"
    ),
    bold_label = (
      (route == "T1" & program %in% c("P3", "P8")) |
        (route == "T2" & program %in% c("P5", "P8")) |
        (route == "T3" & program %in% paste0("P", c(2, 3, 6, 7, 9))) |
        (route == "T4" & program %in% c("P3", "P9"))
    ),
    label = if_else(!agreement | bold_label, program, NA_character_),
    fontface = if_else(bold_label, "bold", "plain")
  )

n_total <- sum(!is.na(directions$agreement))
n_same <- sum(directions$agreement, na.rm = TRUE)

concordance_label <- sprintf(
  "Directional concordance: %d/%d (%.1f%%)",
  n_same, n_total, 100 * n_same / n_total
)

p <- ggplot(directions, aes(source, external, color = family)) +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = 0, ymax = Inf,
    fill = "#DCE9E3", alpha = 0.35
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = 0,
    ymin = -Inf, ymax = 0,
    fill = "#DCE9E3", alpha = 0.35
  ) +
  annotate(
    "text",
    x = 0, y = Inf,
    label = concordance_label,
    hjust = 0.57,
    vjust = 1.5,
    size = 1.8
  ) +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = "22", color = "#c6c6c6") +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "22", color = "#c6c6c6") +
  geom_abline(slope = 1, linewidth = 0.2, intercept = 0, color = "#c6c6c6") +
  geom_point(size = 0.3) +
  geom_text_repel(
    data = subset(directions, !agreement),
    aes(label = paste0(route, ":", program)),
    size = 1.8,
    fontface = "plain",
    box.padding = 0.1,
    point.padding = 0,
    force = 2,
    force_pull = 0.5,
    max.overlaps = Inf,
    max.time = 3,
    segment.size = 0.15,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      family_colors,
      "Other" = "#000000"
    ),
    breaks = c("ReprPC", "Bivalent", "Other")
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = 0.2)
  ) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = 0.2)
  ) +
  labs(
    x = "Source prediagnostic − control",
    y = "External LUAD − control",
    color = NULL,
    tag = "g"
  ) +
  basic +
  theme(
    plot.margin = margin(l = 15, r = 5, b = 5, t = 5),
    plot.tag.position = c(-0.03, 0.95),
    legend.justification = "top",
    legend.box.spacing = unit(1, "mm")
  ); ggsave(
  file.path("result/prog-external.svg"),
  device = svglite,
  p,
  width = 3.2,
  height = 1.6,
  units = "in",
  bg = "transparent"
)

########################################################################
# Preferred-route proportions using the pipeline's frozen external posterior.
route_counts <- bind_rows(
  data.frame(stage = "Source", group = ifelse(source_control, "Control", "PreDx/Case"),
             route = source_route),
  data.frame(stage = "External", group = ifelse(ext_control, "Control", "PreDx/Case"),
             route = E$route)) %>%
  count(stage, group, route) %>% group_by(stage, group) %>%
  mutate(proportion = n / sum(n)) %>% ungroup()
route_counts$stage <- factor(route_counts$stage,
                             levels = c("Source", "External"))

condition_cols <- c(
  T4 = "#3F3F3F",
  T1 = "#5E49B8",
  T2 = "#008C8C",
  T3 = "#CF9000"
)

p <- ggplot(route_counts, aes(stage, proportion, fill = route)) +
  geom_col(alpha = 0.8) +
  facet_wrap(~group, ncol = 2) +
  scale_fill_manual(
    values = condition_cols,
    breaks = c("T1", "T2", "T3", "T4")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::label_percent()
  ) +
  labs(
    x = NULL,
    y = "Proportions",
    fill = "Preferred\nroute",
    tag = "h"
  ) +
  basic +
  theme(
    strip.text = element_text(size = 6.5),
    legend.title = element_text(
      margin = margin(b = 2, r = 4)
    ),
    legend.justification = "top",
    legend.box.spacing = unit(0, "mm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# 4. Save ----------------------------------------------------------------------

ggsave(
  file.path("result/external-prop.svg"),
  device = svglite,
  p,
  width = 2.2,
  height = 1.6,
  units = "in",
  bg = "transparent"
)
