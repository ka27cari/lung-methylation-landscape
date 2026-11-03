# Plot outputs for: 01 programs
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# Section: Recurrent methylation programs (original 01_recur.r)
source(file.path(repo, "r", "figures", "figure_environment.R"))
selected_file <- file.path(
  build_dir,
  "cache/model/SELECTED_PROGRAMS.csv"
)

families_file <- file.path(
  build_dir,
  "cache/discovery/RECURRENT_FAMILIES_ALL.csv"
)

audit_file <- file.path(
  build_dir,
  "cache/model/PROGRAM_SELECTION_AUDIT.csv"
)

selected <- read_csv(
  selected_file,
  show_col_types = FALSE
)

families <- read_csv(
  families_file,
  show_col_types = FALSE
)

audit <- read_csv(
  audit_file,
  show_col_types = FALSE
)

xlsx_file <- file.path(manifest_dir, "41467_2025_59333_MOESM34_ESM.xlsx")
source_file <- file.path(build_dir, "cache/discovery/SOURCE_CPG_ORIENTATION.csv")

stopifnot(file.exists(xlsx_file), file.exists(source_file))

# Only strip a validated EPICv2 suffix; non-cg records cannot join Kirmani cg IDs.
base_id <- function(x) {
  z <- ifelse(grepl("^cg[0-9]{8}(_.*)?$", x), sub("_.*$", "", x), NA_character_)
  z
}

src <- readr::read_csv(source_file, show_col_types = FALSE) %>%
  mutate(base = base_id(cpg))
stopifnot(all(c("cpg", "direction", "base") %in% names(src)))
src_unique <- src %>% filter(!is.na(base)) %>%
  group_by(base) %>% summarise(n_signs = n_distinct(sign(direction)),
                               direction = if (n_distinct(sign(direction)) == 1) first(sign(direction)) else NA_real_,
                               .groups = "drop")
if (any(src_unique$n_signs > 1)) warning("Ambiguous direction for duplicate base IDs; excluded from +/- only")

# Reproduce the background exactly as in code/06_reproduction.r:
# EPICv2 probes also annotated as legacy 450K; unique CpG base IDs.
library(missMethyl)
anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
stopifnot(all(c("Name", "Methyl450_Loci") %in% names(anno)))
bg <- anno %>% filter(!is.na(Methyl450_Loci), Methyl450_Loci != "") %>%
  transmute(base = base_id(Name)) %>% filter(!is.na(base)) %>% pull(base) %>% unique()

esm <- readxl::read_xlsx(xlsx_file)
stopifnot(all(c("chip_group", "CpG") %in% names(esm)))
info <- esm %>% filter(chip_group %in% c("any", "ASXL1", "DNMT3A", "TET2")) %>%
  transmute(regulator = as.character(chip_group), base = as.character(CpG)) %>%
  filter(base %in% bg) %>% distinct()
chip_sets <- split(info$base, info$regulator)

# Same one-sided Fisher test and Haldane odds-ratio estimator as original code.
fisher_enrich <- function(query, bg, chip) {
  query <- intersect(unique(query), unique(bg))
  out <- setdiff(unique(bg), query)
  a <- sum(query %in% chip); b <- length(query) - a
  c <- sum(out %in% chip);   d <- length(out) - c
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  p <- stats::fisher.test(mat, alternative = "greater")$p.value
  or <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  tibble(a = a, b = b, c = c, d = d, OR_Haldane = or, p = p)
}

conflicted::conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(dplyr::first)

sets <- list(All_7763 = src_unique$base,
             Prediagnostic_up = src_unique$base[src_unique$direction > 0],
             Prediagnostic_down = src_unique$base[src_unique$direction < 0])
results <- purrr::imap_dfr(sets, function(query, name) {
  purrr::imap_dfr(chip_sets, function(chip, regulator) {
    fisher_enrich(query, bg, chip) %>% mutate(set = name, regulator = regulator, .before = 1)
  })
}) %>% group_by(set) %>% mutate(q_per_set = p.adjust(p, method = "BH")) %>% ungroup() %>%
  mutate(q_all_tests = p.adjust(p, method = "BH"), n_bg = a + b + c + d,
         n_query = a + b, n_chip_bg = a + c,
         observed_percent = ifelse(n_query > 0, 100 * a / n_query, NA_real_),
         expected_percent = 100 * n_chip_bg / n_bg) %>%
  arrange(factor(set, levels = names(sets)), factor(regulator, levels = names(chip_sets)))

chip_genes <- c("ASXL1", "DNMT3A", "TET2")
pal <- c(
  "ASXL1"  = "#7cb8ba",
  "DNMT3A" = "#5b7e8a",
  "TET2"   = "#db6d5f"
)
tr <- pseudo_log_trans(sigma = 0.05)

ovr <- results %>%
  filter(set == "All_7763", regulator %in% chip_genes) %>%
  mutate(
    regulator = factor(regulator, levels = chip_genes),
    x = as.numeric(regulator),
    p_obs = observed_percent / 100,
    p_exp = expected_percent / 100,
    star = p_to_star(q_all_tests)
  )

bar <- ovr %>%
  transmute(regulator, x, prop = p_obs, key = as.character(regulator))

exp <- ovr %>%
  transmute(regulator, x, prop = p_exp)

sig <- ovr %>%
  filter(!is.na(q_all_tests), q_all_tests < 0.05, a > 2) %>%
  mutate(
    base = tr$transform(pmax(p_exp, p_obs)),
    y = tr$inverse(base + 0.12)
  )

enr <- ggplot(bar, aes(x = x, y = prop, fill = key)) +
  geom_col(width = 0.45) +
  geom_segment(
    data = exp,
    aes(x = x - 0.27, xend = x + 0.27, y = prop, yend = prop, color = "Expected"),
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
  scale_y_continuous(
    trans = tr,
    limits = c(0, 0.02),
    breaks = c(0, 0.005, 0.01, 0.02),
    labels = percent_format(accuracy = 0.1),
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
    name = NULL,
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
    guide = guide_legend(order = 1, override.aes = list(linewidth = 0.45))
  ) +
  labs(x = NULL, y = "CHIP-associated", tag = "c") +
  basic +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.spacing = unit(2, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6, margin = margin(l = 2)),
    legend.key.height = unit(2.3, "mm"),
    legend.key.width = unit(2.3, "mm"),
    legend.key.spacing.x = unit(0, "mm"),
    legend.key.spacing.y = unit(0.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(2, "mm"),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

# Direction concordance: same / opposite
unique_direction <- function(v) {
  s <- unique(sign(v[is.finite(v) & v != 0]))
  if (length(s) == 1L) s[[1L]] else NA_real_
}

chip_direction <- esm %>%
  filter(chip_group %in% chip_genes) %>%
  transmute(regulator = as.character(chip_group),
            base = base_id(as.character(CpG)),
            effect = as.numeric(beta)) %>%
  filter(!is.na(base), base %in% bg) %>%
  group_by(regulator, base) %>%
  summarise(chip_sign = unique_direction(effect), .groups = "drop")

direction_summary <- src_unique %>%
  filter(base %in% bg) %>%
  select(base, source_sign = direction) %>%
  inner_join(chip_direction, by = "base") %>%
  mutate(
    comparable = !is.na(source_sign) & !is.na(chip_sign),
    same = comparable & source_sign == chip_sign
  ) %>%
  group_by(regulator) %>%
  summarise(n_comparable = sum(comparable),
            n_same = sum(same, na.rm = TRUE), .groups = "drop")

table_numbers <- ovr %>%
  mutate(regulator = as.character(regulator)) %>%
  left_join(direction_summary, by = "regulator") %>%
  mutate(
    n_comparable = replace_na(n_comparable, 0L),
    n_same = replace_na(n_same, 0L),
    same_pct = if_else(n_comparable > 0, round(100 * n_same / pmax(n_comparable, 1)), NA_real_),
    so_label = if_else(n_comparable > 0,
                       sprintf("%d/%d", as.integer(same_pct), as.integer(100 - same_pct)), "–"),
    so_label = if_else(so_label == "100/0", "All", so_label)
  )

stopifnot(all(table_numbers$n_comparable <= table_numbers$a),
          n_distinct(table_numbers$n_query) == 1)

set_total <- table_numbers %>%
  summarise(Total = first(n_query))

# Original text-only table
chip_table_plot <- ggplot() +
  geom_hline(yintercept = c(2.5, 1.5, 0.5),
             linewidth = 0.18, colour = "black") +
  geom_text(data = table_numbers,
            aes(x = x, y = 3, label = a),
            size = 2, family = "Arial") +
  geom_text(data = set_total,
            aes(x = 2, y = 2, label = Total),
            size = 2, family = "Arial") +
  geom_text(data = table_numbers,
            aes(x = x, y = 1, label = so_label),
            size = 2, family = "Arial") +
  scale_x_continuous(limits = c(0.5, 3.5), expand = expansion(add = 0)) +
  scale_y_continuous(
    limits = c(0.5, 3.45),
    breaks = c(1, 2, 3),
    labels = c("S/O (%)", "Total", "Overlap"),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(x = NULL, y = NULL) +
  basic +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(
      size = 6, hjust = 1, margin = margin(r = 2)
    ),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

chip_enrichment_table <- (enr / chip_table_plot) +
  patchwork::plot_layout(heights = c(1.0, 0.5)) + plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = NA,
        colour = NA
      )
    )
  ) &
  theme(
    plot.margin = margin(l = 7, t = 5, b = 1, r = 2),
    plot.tag.location = "panel",
    plot.tag.position = c(-0.75, 1),
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14
    )
  )

ggsave(
  "result/chip-cpgs.svg",
  plot = chip_enrichment_table,
  device = svglite::svglite,
  width = 1.9,
  height = 1.5,
  bg = "transparent"
)

# Effective dimensions ────────────────────────────────────────────────────────

file <- file.path(build_dir, "anchor/selection/anchor/ANCHOR_PR10_PR60_RECOMPUTED.csv")
if (!file.exists(file)) stop("Anchor PCA sensitivity CSV not found.")

x <- readr::read_csv(file, show_col_types = FALSE)
stopifnot(all(c("anchor", "pr10", "pr20", "pr30", "pr40") %in% names(x)))

pca_df <- x %>%
  tidyr::pivot_longer(cols = c(pr10, pr20, pr30, pr40),
                      names_to = "pc", values_to = "pr") %>%
  dplyr::mutate(anchor_pct = anchor * 100,
                pc = factor(pc,
                            levels = c("pr10", "pr20", "pr30", "pr40"),
                            labels = c("PC10", "PC20", "PC30", "PC40")))

min_df <- pca_df %>%
  dplyr::group_by(pc) %>%
  dplyr::slice_min(pr, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

stopifnot(nrow(min_df) == 4, all(abs(min_df$anchor_pct - 2) < 1e-8))
anchor_n <- x$n_anchor[abs(x$anchor - 0.02) < 1e-8]

stopifnot(nrow(min_df) == 4, all(abs(min_df$anchor_pct - 2) < 1e-8))

p_anchor <- ggplot(pca_df, aes(x = anchor_pct, y = pr, color = pc, group = pc)) +
  geom_vline(xintercept = 2, linetype = "22", linewidth = 0.25, color = "grey30") +
  geom_line(linewidth = 0.5) +
  annotate("text", x = 2.2, y = Inf,
           label = sprintf("2%% selected (%d CpGs)", anchor_n),
           hjust = 0, vjust = 1.2, size = 1.8,
           family = "Arial") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10),
                labels = c("0.5", "1", "2", "5", "10"),
                limits = c(0.45, 15)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.6),
                     breaks = c(2, 5, 10, 15, 20),
                     expand = expansion(mult = c(0.04, 0.1))) +
  scale_color_manual(
    values = c(
      "PC10" = "#284B63",
      "PC20" = "#56839A",
      "PC30" = "#659C98",
      "PC40" = "#9889A8"
    ),
    name = "Leading PCs"
  ) +
  labs(x = "Reference anchor (%)",
       y = "Effective dimensions",
       tag = "d") +
  basic + theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = -4, l = 0),
    legend.key.spacing.x = unit(2, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6, margin = margin(l = 1))
  )

# Save
ggsave("result/anchor.svg", device = svglite, p_anchor, width = 1.9, height = 1.6, dpi = 800)

# Restrict to the actual program-selection universe ────────────────────────────
depth_df <- families %>%
  filter(abs(cosine_threshold - 0.80) < 1e-8,
         strategy == "kmeans_local",
         representation == "counterfactual_a0.02_k10_rank") %>%
  mutate(stability_group = case_when(
    n_ranks < 5 ~ "Partial",
    seed_stability_mean >= 0.95 ~ "High",
    TRUE ~ "Other"
  )) %>%
  count(n_ranks, stability_group, name = "n_family") %>%
  complete(n_ranks = 2:5, stability_group = c("High", "Other", "Partial"), fill = list(n_family = 0)) %>%
  mutate(stability_group = factor(stability_group, levels = c("High", "Other", "Partial")))

label_df <- depth_df %>%
  group_by(n_ranks) %>%
  summarise(n_total = sum(n_family),
            n_high = sum(n_family[stability_group == "High"]),
            n_low = sum(n_family[stability_group == "Other"]), .groups = "drop")

# Dashed connectors from the high-stability portion
high_bar <- label_df %>% filter(n_ranks == 5)
stopifnot(nrow(high_bar) == 1)

y_max <- max(label_df$n_total) * 1.18

connectors <- tibble(
  x = 5.35,
  y = c(high_bar$n_total, high_bar$n_low),
  xend = 6,
  yend = c(y_max * 0.75, y_max * 0.05)
)

p_depth <- ggplot(depth_df, aes(n_ranks, n_family, fill = stability_group)) +
  geom_col(width = 0.7) +
  geom_text(data = label_df, aes(x = n_ranks, y = n_total, label = n_total),
            inherit.aes = FALSE, vjust = -0.4, size = 2) +
  geom_text(data = high_bar,
            aes(x = n_ranks, y = n_low + n_high / 2, label = n_high),
            inherit.aes = FALSE, color = "white", size = 2) +
  annotate("text", x = 3.05, y = Inf,
           label = "Ranks: 12, 14, 16, 20, 24",
           vjust = 1.5, size = 1.9) +
  geom_segment(data = connectors,
               aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, linetype = "22",
               linewidth = 0.2, color = "black") +
  scale_fill_manual(
    values = c("Partial" = "grey75", "Other" = "#9db3c2", "High" = "#2a5a7a"),
    breaks = c("Other", "High"),
    labels = c("<0.95", "≥0.95"),
    name = "Seed stability",
    guide = guide_legend(nrow = 1, title.position = "left")
  ) +
  scale_x_continuous(breaks = 2:5, labels = paste0(2:5, "/5")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(xlim = c(1.55, 5.45),
                  ylim = c(0, max(label_df$n_total) * 1.18),
                  expand = FALSE, clip = "off") +
  labs(x = "Recurrence across NMF ranks", y = "Candidate programs (n)", tag = "e") +
  basic + theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = -4, l = 0),
    legend.key.spacing.x = unit(2, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.title = element_text(size = 6, margin = margin(r = 5)),
    legend.text = element_text(size = 6, margin = margin(l = 2)),
    plot.margin = margin(l = 10, r = 60, b = 5, t = 5)
  )

# Save
ggsave("result/fam-depth.svg", device = svglite, p_depth, width = 2.5, height = 1.6, dpi = 800)

root <- base
outdir <- file.path(root, "06_audit")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidates <- file.path(build_dir, "analysis/tables/PROGRAM_DUAL_RESOLUTION_COVERAGE_CURVE.csv")
file <- candidates[file.exists(candidates)][1]
if (is.na(file)) stop("Coverage curve CSV not found.")

x <- readr::read_csv(file, show_col_types = FALSE) %>%
  arrange(n_programs)

stopifnot(all(c("n_programs", "uncovered_n") %in% names(x)))
stopifnot(all(is.finite(x$n_programs)), all(is.finite(x$uncovered_n)))
stopifnot(!anyDuplicated(x$n_programs))

p8 <- filter(x, n_programs == 8)
p9 <- filter(x, n_programs == 9)
stopifnot(nrow(p8) == 1, nrow(p9) == 1)

# Visualization
p_basis <- ggplot(x, aes(x = n_programs, y = uncovered_n)) +
  geom_hline(yintercept = 20, linetype = "22",
             linewidth = 0.2, color = "gray30") +
  geom_line(linewidth = 0.45, color = "#8199AE") +
  geom_point(shape = 16, size = 0, color = "#8199AE") +
  geom_point(data = p9, shape = 16, size = 2, color = "#2A5A7A") +
  geom_text(
    data = p9, aes(label = "9 selected"),
    nudge_x = 0.35, hjust = 0, vjust = -0.7,
    size = 2, family = "Arial", color = "#2A5A7A"
  ) +
  annotate(
    "text", x = 0, y = 20, label = "Cutoff: 20",
    hjust = 0, vjust = -0.65,
    size = 1.9, family = "Arial"
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 5),
    breaks = c(0, 7, 20, 100, 600),
    expand = expansion(mult = c(0.04, 0.14)),
    position = "right"
  ) +
  labs(
    x = "Programs",
    y = "Uncovered samples"
  ) +
  basic +
  theme(
    legend.position = "none",
    axis.line.y.left = element_blank(),
    axis.line.y.right = element_line(linewidth = 0.35)
  )

# Save
ggsave("result/basis_selection.svg", device = svglite, p_basis, width = 1.6, height = 1.2, dpi = 800)

# Stability prioritization among all 5/5-rank families ─────────────────────────
full_rank <- fam_use %>%
  filter(n_ranks == length(ranks_use)) %>%
  arrange(desc(seed_stability_mean), desc(context_n), family_id) %>%
  mutate(candidate_rank = row_number()) %>%
  left_join(
    selected %>% select(family_id, program_index, selection_rank),
    by = "family_id"
  ) %>%
  left_join(
    audit %>% select(family_id, decision, max_cosine_to_selected),
    by = "family_id"
  ) %>%
  mutate(
    status = case_when(
      !is.na(program_index) ~ "Selected",
      decision == "collapsed_near_duplicate" ~ "Near-duplicate",
      TRUE ~ "Other full-rank"
    )
  )

# position at which selection stopped
last_examined_rank <- full_rank %>%
  filter(status %in% c("Selected", "Near-duplicate")) %>%
  summarise(x = max(candidate_rank)) %>%
  pull(x)

p_stability <- ggplot(full_rank, aes(x = candidate_rank, y = seed_stability_mean)) +
  annotate(
    "rect",
    xmin = 0.5, xmax = last_examined_rank + 0.5,
    ymin = -Inf, ymax = Inf,
    fill = "#8da8ba", alpha = 0.25, color = NA
  ) +
  geom_line(color = "grey80", linewidth = 0.15) +
  geom_line(
    data = full_rank %>% filter(candidate_rank <= last_examined_rank),
    aes(x = candidate_rank, y = seed_stability_mean),
    inherit.aes = FALSE,
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = last_examined_rank + 0.5,
    linetype = "dashed", linewidth = 0.2
  ) +
  geom_point(
    data = full_rank %>% filter(status == "Other full-rank"),
    shape = 16, size = 0.45, color = "grey80"
  ) +
  geom_point(
    data = full_rank %>% filter(status == "Selected"),
    aes(shape = k, color = k),
    size = 0.85
  ) +
  annotate(
    "text",
    x = sqrt(last_examined_rank),
    y = 0.902,
    label = "10 candidates",
    hjust = 0.5, vjust = 0,
    size = 1.9,
    color = "grey25"
  ) +
  scale_shape_manual(
    values = c("k4" = 16, "k6" = 17),
    breaks = context_breaks,
    labels = context_labels,
    name = "Local context"
  ) +
  scale_color_manual(
    values = cols_k,
    breaks = context_breaks,
    labels = context_labels,
    name = "Local context"
  ) +
  scale_x_continuous(
    breaks = c(1, 10, 25, 50, 75, 100),
    labels = c("first", "10", "25", "50", "75", "100"),
    limits = c(0.5, 75),
    expand = 0,
    trans = "pseudo_log"
  ) +
  scale_y_continuous(
    limits = c(0.9, 1.005),
    breaks = c(0.9, 1.0)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = c(16, 17), size = 1.1)
    )
  ) +
  labs(
    x = "Candidate rank",
    y = "Seed stability",
    tag = "d"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    
    plot.margin = margin(l = 10, r = 5, b = 5, t = 5),
    plot.tag = element_text(
      family = "Arial", face = "bold", size = 16,
      margin = margin(l = -5)
    ),
    plot.tag.location = "plot",
    plot.tag.position = c(0, 0.95),
    
    legend.title = element_text(
      size = 5.5, hjust = 0.5,
      margin = margin(b = 1.5)
    ),
    legend.text = element_text(
      size = 5.3, hjust = 0.5,
      margin = margin(l = 0)
    ),
    
    legend.key.width = unit(2.5, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.spacing.y = unit(0.3, "mm"),
    legend.spacing.x = unit(1, "mm"),
    legend.box.spacing = unit(0, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.just = "center",
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.95),
    legend.justification = c(1, 1),
    legend.direction = "vertical"
  )

# Save
ggsave("result/stability.png", p_stability, width = 1.85, height = 1.5, dpi = 800)

# Sequential redundancy pruning ────────────────────────────────────────────────
audit_df <- audit %>%
  filter(
    final_decision %in% c(
      "selected_geometry_basis",
      "collapsed_near_duplicate"
    )
  ) %>%
  left_join(
    fam_use %>% select(family_id, context_n, k),
    by = "family_id"
  ) %>%
  mutate(
    cosine_plot = if_else(
      max_cosine_to_selected < 0,
      0,
      max_cosine_to_selected
    ),
    status = case_when(
      final_decision == "selected_geometry_basis" ~ "Selected",
      final_decision == "collapsed_near_duplicate" ~ "Near-duplicate",
      TRUE ~ NA_character_
    ),
    program_index = if_else(
      final_decision == "selected_geometry_basis",
      selected_program_index,
      NA_real_
    )
  )

redundant <- ggplot(audit_df, aes(x = cosine_plot, y = seed_stability)) +
  geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.2) +
  geom_point(
    data = audit_df %>% filter(status == "Selected"),
    aes(shape = k, color = k),
    size = 2.3
  ) +
  geom_text(
    data = audit_df %>% filter(status == "Selected"),
    aes(label = program_index),
    size = 1.85, color = "white"
  ) +
  geom_point(
    data = audit_df %>% filter(status == "Near-duplicate"),
    shape = 4, size = 2, stroke = 0.35, color = "red"
  ) +
  annotate(
    "text",
    x = 1.0, y = 0.996,
    label = "near-duplicate",
    hjust = 1.12, vjust = 0.5,
    size = 1.7, color = "red"
  ) +
  scale_shape_manual(
    values = c("k4" = 16, "k6" = 17),
    guide = "none"
  ) +
  scale_color_manual(
    values = cols_k,
    breaks = c("k4", "k6"),
    labels = parse(text = c("italic(k) == 4", "italic(k) == 6")),
    # name = "Local context",
    name = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(16, 17),
        size = 1.1
      )
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 0.7, 0.99),
    labels = c("first", "0.5", "0.7", "0.99"),
    limits = c(-0.03, 1.05)
  ) +
  scale_y_continuous(
    limits = c(0.96, 1.001),
    breaks = c(0.96, 0.98, 1.00)
  ) +
  labs(
    x = "Max. cosine to others",
    y = "Seed stability",
    tag = "d"
  ) +
  basic + theme(
    legend.position = "top",
    legend.text = element_text(
      size = 5.5,
      margin = margin(l = 2)
    ),
    legend.box.margin = margin(b = -5),
    legend.key.height = unit(2, "mm"),
    legend.key.width  = unit(1.4, "mm"),
    legend.spacing.x  = unit(1, "mm")
  )

# Save
ggsave("result/redundancy.png", redundant, width = 1.65, height = 1.5, dpi = 800)

# Section: Program ontology and annotation (original 02_ontology.r)
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)
library(ape)
library(ggplotify)
library(ggtext)
library(patchwork)

# Paths ------------------------------------------------------------------------
landscape_dir <- dirname(build_dir)
ontology_dir  <- file.path(landscape_dir, "ontology")

cosine_file <- file.path(
  ontology_dir,
  "results",
  "hierarchy_cosine_recomputed_in_R.csv"
)

mm_dir_500 <- file.path(ontology_dir, "results", "top500")
mm_dir_1000 <- file.path(ontology_dir, "results", "top1000")

kycg_file <- file.path(
  ontology_dir,
  "kycg",
  "results",
  "KYCG_RECURRENT_FEATURES.csv"
)

# Settings ---------------------------------------------------------------------
display_order <- c("P1", "P2", "P6", "P3", "P9", "P7", "P5", "P8", "P4")

families <- list(
  "ReprPC" = c("P3", "P9"),
  "Bivalent" = c("P5", "P8")
)

family_colors <- c(
  "ReprPC" = "#292473",
  "Bivalent" = "#C2185B"
)

# Read cosine matrix -----------------------------------------------------------
read_cosine_matrix <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  p <- paste0("P", 1:9)
  
  if (all(p %in% names(x))) {
    other_cols <- setdiff(names(x), p)
    
    for (row_col in other_cols) {
      rn <- as.character(x[[row_col]])
      
      if (length(rn) == 9 && all(rn %in% p)) {
        m <- as.matrix(x[, p])
        storage.mode(m) <- "numeric"
        rownames(m) <- rn
        m <- m[p, p, drop = FALSE]
        m <- (m + t(m)) / 2
        diag(m) <- 1
        return(m)
      }
    }
    
    if (nrow(x) == 9) {
      m <- as.matrix(x[, p])
      storage.mode(m) <- "numeric"
      rownames(m) <- p
      m <- (m + t(m)) / 2
      diag(m) <- 1
      return(m)
    }
  }
  
  stop("Cannot identify cosine matrix format.")
}

cosine <- read_cosine_matrix(cosine_file)
cosine <- cosine[display_order, display_order, drop = FALSE]
n_program <- length(display_order)

# Heatmap data -----------------------------------------------------------------
heat_df <- expand.grid(
  row_program = display_order,
  col_program = display_order,
  stringsAsFactors = FALSE
) %>%
  mutate(
    row_i = match(row_program, display_order),
    col_i = match(col_program, display_order),
    similarity = mapply(function(r, c) cosine[r, c], row_program, col_program),
    x = col_i,
    y = n_program + 1 - row_i
  ) %>%
  filter(row_i <= col_i)

# Text colour for values
sim_min <- min(heat_df$similarity, na.rm = TRUE)

heat_df <- heat_df %>%
  mutate(
    similarity_label = sprintf("%.2f", similarity),
    text_color = if_else(similarity >= 0.75, "white", "black")
  )

# Custom palette
heat_colors <- c(
  "#FFFFFF",
  "#F7FBFF",
  "#b5dcf5",
  "#90bfde",
  "#5488ab",
  "#2a5a7a"
)

low_sim_df <- heat_df %>%
  filter(row_i < col_i, similarity <= 0.5)

# Family outlines --------------------------------------------------------------
make_pair_outline <- function(family) {
  idx <- sort(match(families[[family]], display_order))
  if (length(idx) != 2 || diff(idx) != 1) stop(family, " members are not adjacent.")
  
  a <- idx[1]
  b <- idx[2]
  ya <- n_program + 1 - a
  yb <- n_program + 1 - b
  
  tibble(
    family = family,
    x = c(a - .5, b + .5, b + .5, b - .5, b - .5, a - .5),
    y = c(ya + .5, ya + .5, yb - .5, yb - .5, ya - .5, ya - .5)
  )
}

outline_df <- bind_rows(lapply(names(families), make_pair_outline))

# Family labels: beside upper block --------------------------------------------
family_label_df <- bind_rows(lapply(names(families), function(family) {
  idx <- sort(match(families[[family]], display_order))
  top_idx <- idx[1]
  
  tibble(
    family = family,
    x = top_idx - 0.60,
    y = n_program + 1 - top_idx,
    label = family
  )
}))

# Plot -------------------------------------------------------------------------
p_correl <- ggplot(heat_df, aes(x, y, fill = similarity)) +
  geom_tile(color = "black", linewidth = 0.2) +
  geom_polygon(
    data = outline_df,
    aes(x, y, group = family, color = family),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.75,
    linejoin = "mitre"
  ) +
  geom_text(
    data = family_label_df,
    aes(x, y, label = label, color = family),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 0.5,
    size = 2.2
  ) +
  geom_point(
    data = low_sim_df,
    aes(x, y),
    inherit.aes = FALSE,
    shape = 4,
    color = "black",
    size = 0.75,
    stroke = 0.4
  ) +
  scale_fill_gradientn(
    name = "Cosine similarity",
    colours = heat_colors,
    limits = c(sim_min, 1),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      title.hjust = 0,
      barwidth = grid::unit(13, "mm"),
      barheight = grid::unit(1.5, "mm"),
      ticks = FALSE,
      theme = theme(
        legend.title = element_text(size = 6, margin = margin(b = 2)),
        legend.text = element_text(size = 5.5, margin = margin(t = 1))
      )
    )
  ) +
  scale_color_manual(values = family_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq_len(n_program), labels = display_order,
    limits = c(.5, n_program + .5), expand = c(0, 0), position = "top"
  ) +
  scale_y_continuous(
    breaks = seq_len(n_program), labels = rev(display_order),
    limits = c(.5, n_program + .5), expand = c(0, 0), position = "right"
  ) +
  labs(x = NULL, y = NULL, tag = "f") +
  basic + 
  theme(
    axis.line = element_blank(),
    plot.margin = margin(l = 20, r = 5, b = 5, t = 5),
    axis.ticks = element_blank(),
    plot.tag.location  = "plot",
    plot.tag.position  = c(-0.05, 0.95),
    legend.position = "inside",
    legend.position.inside = c(-0.05, -0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = NA, color = NA),
    legend.margin = margin(0, 0, 0, 0)
  )

ggsave("result/program-correl.svg", device = svglite, p_correl, width = 1.7, height = 1.5, dpi = 800)

# Additional paths/settings ----------------------------------------------------
mm_dir_500 <- file.path(ontology_dir, "results", "top500")
mm_dir_1000 <- file.path(ontology_dir, "results", "top1000")
kycg_file <- file.path(ontology_dir, "kycg", "results", "KYCG_RECURRENT_FEATURES.csv")

family_source_target <- c("ReprPC" = "T2", "Bivalent PC" = "T3")
mm_fdr_cut <- 0.05
kycg_fdr_cut <- 0.05
kycg_min_runs <- 3
max_ontology_terms <- 16

forced_mm_terms <- c(
  "homophilic cell adhesion via plasma membrane adhesion molecules",
  "Cadherin signaling"
)

# Unrooted program tree =====================================================
distance <- 1 - cosine
diag(distance) <- 0
distance[distance < 0] <- 0

hc <- hclust(as.dist(distance), method = "average")
phy <- as.phylo(hc)

get_descendant_tips <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) return(node)
  children <- tree$edge[tree$edge[, 1] == node, 2]
  unlist(lapply(children, function(z) get_descendant_tips(tree, z)), use.names = FALSE)
}

edge_family <- rep(NA_character_, nrow(phy$edge))

for (i in seq_len(nrow(phy$edge))) {
  child <- phy$edge[i, 2]
  labs <- phy$tip.label[unique(get_descendant_tips(phy, child))]
  
  for (family in names(families)) {
    if (all(labs %in% families[[family]])) {
      edge_family[i] <- family
      break
    }
  }
}

edge_color <- ifelse(is.na(edge_family), "black", family_colors[edge_family])
edge_width <- ifelse(is.na(edge_family), 0.65, 1.2)

tip_family <- rep(NA_character_, length(phy$tip.label))
names(tip_family) <- phy$tip.label
for (family in names(families)) tip_family[phy$tip.label %in% families[[family]]] <- family

tip_color <- ifelse(is.na(tip_family), "black", family_colors[tip_family])

p_tree <- ggplotify::as.ggplot(~{
  par(mar = c(1.2, 1.8, 1.2, 1.8), xpd = NA)
  
  ape::plot.phylo(
    phy,
    type = "unrooted",
    use.edge.length = FALSE,
    show.tip.label = FALSE,
    edge.color = edge_color,
    edge.width = edge_width,
    no.margin = FALSE
  )
  
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tip <- length(phy$tip.label)
  
  tree_span <- max(diff(range(pp$xx)), diff(range(pp$yy)))
  offset <- tree_span * 0.055
  
  for (i in seq_len(n_tip)) {
    parent <- phy$edge[match(i, phy$edge[, 2]), 1]
    
    dx <- pp$xx[i] - pp$xx[parent]
    dy <- pp$yy[i] - pp$yy[parent]
    d <- sqrt(dx^2 + dy^2)
    
    ux <- dx / d
    uy <- dy / d
    
    lx <- pp$xx[i] + offset * ux
    ly <- pp$yy[i] + offset * uy
    
    if (ux > 0.20) {
      adj <- c(0, 0.5)
    } else if (ux < -0.20) {
      adj <- c(1, 0.5)
    } else if (uy > 0) {
      adj <- c(0.5, 0)
    } else {
      adj <- c(0.5, 1)
    }
    
    text(
      lx, ly,
      labels = phy$tip.label[i],
      col = tip_color[i],
      font = 1,
      cex = 0.55,
      adj = adj,
      xpd = NA
    )
  }
}) +
  labs(tag = "a") +
  theme(
    plot.margin = margin(l = 10, r = 5, b = 5, t = 5),
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      margin = margin(l = -5)
    ),
    plot.tag.location = "plot",
    plot.tag.position = c(0, 0.95)
  )

ggsave("result/program-tree.svg", device = svglite, p_tree, width = 1, height = 1.4, dpi = 800)

# Common annotation order ======================================================
program_order <- c("P1", "P2", "P6", "P3", "P9", "P4", "P7", "P5", "P8")
target_order <- c("P1", "P2", "P6", "P3", "P9", "T2", "P4", "P7", "P5", "P8", "T3")

program_role_colors <- c(
  P1 = "#D88900",
  P6 = "#7A4E93",
  P4 = "#5C6B7A"
)

target_labels <- c(
  P1 = "P1", P2 = "P2", P6 = "P6", P3 = "P3", P9 = "P9",
  T2 = "P3+P9", P7 = "P7", P5 = "P5", P8 = "P8",
  T3 = "P5+P8", P4 = "P4"
)

mm_fdr_cut <- 0.05
max_ontology_terms <- 12

forced_mm_terms <- c(
  "homophilic cell adhesion via plasma membrane adhesion molecules",
  "Cadherin signaling"
)

target_labels["P6"] <- sprintf(
  "P6<br><span style='color:%s'>ReprPC</span>",
  program_role_colors["P6"]
)

target_labels["T2"] <- sprintf(
  "<span style='color:%s'>P3+P9</span>",
  family_colors["T2"]
)

target_labels["P4"] <- sprintf(
  "P4<br><span style='color:%s'>Heterochromatin</span>",
  program_role_colors["P4"]
)

target_labels["P9"] <- sprintf(
  "<span style='color:%s'>P9</span>",
  family_colors["T2"]
)

target_labels["P7"] <- "P7"

ontology_low  <- "#BCE5E0"
ontology_high <- "#007E87"

get_mm_files <- function(directory, top_n) {
  list.files(
    directory,
    pattern = paste0("__top", top_n, "__(GO|KEGG)__full_EPICv2\\.csv$"),
    full.names = TRUE
  )
}

read_mm_file <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  
  if ("TERM" %in% names(x)) {
    x <- x %>% mutate(source = "GO", term = as.character(TERM))
  } else if ("Description" %in% names(x)) {
    x <- x %>% mutate(source = "KEGG", term = as.character(Description))
  } else {
    return(NULL)
  }
  
  x
}

mm_files <- c(
  get_mm_files(mm_dir_500, 500),
  get_mm_files(mm_dir_1000, 1000)
)

if (!length(mm_files)) stop("No missMethyl result files found.")

mm <- map_dfr(mm_files, read_mm_file) %>%
  mutate(
    FDR = as.numeric(FDR),
    DE = as.numeric(DE),
    top_n_per_program = as.numeric(top_n_per_program)
  ) %>%
  filter(target %in% target_order)

# GO Biological Process only
if ("ONTOLOGY" %in% names(mm)) {
  mm <- mm %>% filter(source != "GO" | ONTOLOGY == "BP")
}

# Remove broader duplicate adhesion term
mm <- mm %>%
  filter(
    !str_detect(
      str_to_lower(term),
      "^cell-cell adhesion via plasma.?membrane adhesion molecules$"
    )
  )

# Significant terms ------------------------------------------------------------
mm_sig <- mm %>%
  filter(
    is.finite(FDR),
    FDR < mm_fdr_cut,
    top_n_per_program %in% c(500, 1000)
  )

# Collapse top500 / top1000.
# Threshold information is retained only for term selection.
mm_sum <- mm_sig %>%
  group_by(target, source, term) %>%
  summarise(
    n_thresholds = n_distinct(top_n_per_program),
    best_FDR = min(FDR),
    DE = max(DE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(score = -log10(pmax(best_FDR, 1e-300)))

# Representative ontology terms -----------------------------------------------
term_stats <- mm_sum %>%
  group_by(source, term) %>%
  summarise(
    n_targets = n_distinct(target),
    max_thresholds = max(n_thresholds),
    total_thresholds = sum(n_thresholds),
    best_FDR = min(best_FDR),
    max_DE = max(DE, na.rm = TRUE),
    .groups = "drop"
  )

robust_terms <- term_stats %>%
  filter(max_thresholds == 2 | n_targets >= 2)

best_per_target <- mm_sum %>%
  group_by(target) %>%
  arrange(desc(n_thresholds), best_FDR, desc(DE), .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  distinct(source, term)

forced_terms <- mm_sum %>%
  filter(term %in% forced_mm_terms) %>%
  distinct(source, term)

selected_terms <- bind_rows(
  forced_terms,
  robust_terms %>% select(source, term),
  best_per_target
) %>%
  distinct(source, term) %>%
  left_join(term_stats, by = c("source", "term")) %>%
  arrange(
    desc(n_targets),
    desc(max_thresholds),
    desc(total_thresholds),
    best_FDR,
    desc(max_DE)
  ) %>%
  slice_head(n = max_ontology_terms) %>%
  select(source, term)

# Plot data --------------------------------------------------------------------
mm_plot <- mm_sum %>%
  semi_join(selected_terms, by = c("source", "term")) %>%
  mutate(
    plot_term = case_when(
      term == "homophilic cell adhesion via plasma membrane adhesion molecules" ~
        "GO  Homophilic cell adhesion",
      TRUE ~ paste0(if_else(source == "GO", "GO:  ", "KEGG:  "), term)
    )
  )

# Remove targets with no displayed ontology term
plot_targets <- target_order[target_order %in% unique(mm_plot$target)]

mm_plot <- mm_plot %>%
  filter(target %in% plot_targets) %>%
  mutate(target = factor(target, levels = plot_targets))

# Order terms: more broadly recurrent terms toward top
term_order <- mm_plot %>%
  group_by(plot_term) %>%
  summarise(
    n_targets = n_distinct(target),
    best_FDR = min(best_FDR),
    max_DE = max(DE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(n_targets, desc(best_FDR), max_DE) %>%
  pull(plot_term)

mm_plot <- mm_plot %>%
  mutate(plot_term = factor(plot_term, levels = unique(term_order)))

# Labels only for targets that remain
plot_target_labels <- target_labels[plot_targets]

# Plot -------------------------------------------------------------------------
p_mm <- ggplot(mm_plot, aes(x = target, y = plot_term)) +
  geom_point(
    aes(size = DE, color = score),
    shape = 16
  ) +
  scale_size_area(
    max_size = 4,
    breaks = c(10, 20, 30),
    name = "DE",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      label.position = "right",
      nrow = 1
    )
  ) +
  scale_color_gradient(
    low = ontology_low,
    high = ontology_high,
    name = expression(-log[10](FDR)),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(13, "mm"),
      barheight = unit(1.5, "mm"),
      ticks = FALSE
    )
  ) +
  scale_x_discrete(
    limits = plot_targets,
    labels = plot_target_labels,
    drop = TRUE
  ) +
  scale_y_discrete(
    labels = function(x) str_wrap(x, width = 42)
  ) +
  labs(x = NULL, y = NULL, tag = "b") +
  basic +
  theme(
    axis.line = element_blank(),
    axis.text.x = ggtext::element_markdown(
      size = 5.5, angle = 0, hjust = 0.5, vjust = 1
    ),
    axis.text.y = element_text(size = 5.5),
    plot.margin = margin(l = 15, r = 5, b = 5, t = 5),
    plot.tag.position = c(-0.03, 0.95),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 5.5, margin = margin(b = 1)),
    legend.text = element_text(size = 5, margin = margin(t = 2)),
    legend.spacing.y = unit(0, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.4
    ),
    legend.spacing.x = unit(5, "mm"),
    legend.key.spacing.x = unit(2, "mm"),
    legend.box.spacing =  unit(1, "mm")
  )

ggsave("result/miss-methyl.svg", device = svglite,  p_mm, width = 3, height = 1.4, dpi = 800)

# knowYourCG regulatory context =======================================================
program_order <- c("P1", "P2", "P7", "P6", "P3", "P9", "P4", "P5", "P8")
target_order <- c("P1", "P3", "P9", "T2", "P2", "P4", "P6", "P5", "P8", "T3")

program_role_colors <- c(
  P1 = "#D88900",
  P6 = "#7A4E93",
  P4 = "#5C6B7A"
)

family_colors <- c(
  T2 = "#292473",
  T3 = "#C2185B"
)

kycg_fdr_cut <- 0.05
kycg_min_runs <- 3

# kycg_low <- "#F4C4D8"
# kycg_high <- "#9A1550"
kycg_low  <- "#bacedb"
kycg_high <- "#49708a"

# X-axis labels ----------------------------------------------------------------
target_labels <- setNames(target_order, target_order)

target_labels["P1"] <- sprintf(
  "P1<br><span style='color:%s'>Promoter / CGI</span>",
  program_role_colors["P1"]
)

target_labels["P2"] <- "P2"

target_labels["P3"] <- "P3"
target_labels["P9"] <- "P9"

target_labels["T2"] <- sprintf(
  "P3+P9<br><span style='color:%s'>ReprPC</span>",
  family_colors["T2"]
)

target_labels["P4"] <- sprintf(
  "P4<br><span style='color:%s'>Heterochromatin</span>",
  program_role_colors["P4"]
)

target_labels["P6"] <- "P6"

target_labels["P5"] <- "P5"
target_labels["P8"] <- "P8"

target_labels["T3"] <- sprintf(
  "P5+P8<br><span style='color:%s'>Bivalent</span>",
  family_colors["T3"]
)

# Read KYCG --------------------------------------------------------------------
kycg <- read_csv(kycg_file, show_col_types = FALSE) %>%
  mutate(
    n_significant_runs = as.numeric(n_significant_runs),
    best_global_FDR = as.numeric(best_global_FDR),
    max_overlap = as.numeric(max_overlap)
  ) %>%
  filter(
    target %in% target_order,
    source_class == "biological_regulatory"
  )

# Regulatory categories --------------------------------------------------------
kycg_def <- tribble(
  ~category, ~pattern,
  "CpG island", "^CGI;Island$",
  "CGI shore", "^CGI;.*Shore$",
  "Active TSS", "^(ChromHMM;TssA|REMCChromHMM;[0-9]+_TssA)$",
  "TSS flank", "(TssFlnk|TssAFlnk)",
  "H3K4me3", "^HM;H3K4me3",
  "ReprPC", "^(ChromHMM;ReprPC|REMCChromHMM;[0-9]+_ReprPC)$",
  "Weak ReprPC", "ReprPCWk$",
  "Bivalent TSS", "TssBiv$",
  "Heterochromatin", "(ChromHMM;Het$|REMCChromHMM;[0-9]+_Het$)",
  "ZNF / repeats", "ZNF/Rpts",
  "PMD", "^PMD;",
  "H3K27me3", "^HM;H3K27me3",
  "H2AK119", "^HM;H2AK119",
  "JARID2", "^TFBSrm;JARID2$",
  "CTCF", "^CTCFbind;"
)

kycg_candidates <- map_dfr(seq_len(nrow(kycg_def)), function(i) {
  kycg %>%
    filter(grepl(kycg_def$pattern[i], dbname, perl = TRUE)) %>%
    mutate(category = kycg_def$category[i])
})

# Significant + robust results -------------------------------------------------
# n_significant_runs is only a filtering criterion, not a visual encoding.
kycg_plot <- kycg_candidates %>%
  filter(
    is.finite(best_global_FDR),
    best_global_FDR < kycg_fdr_cut,
    n_significant_runs >= kycg_min_runs
  ) %>%
  group_by(target, category) %>%
  arrange(
    desc(n_significant_runs),
    best_global_FDR,
    desc(max_overlap),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    target = factor(target, levels = target_order),
    score = -log10(pmax(best_global_FDR, 1e-300))
  )

# Category order ---------------------------------------------------------------
kycg_order <- c(
  "CpG island",
  "Active TSS",
  "H3K4me3",
  "CGI shore",
  "TSS flank",
  "ReprPC",
  "Weak ReprPC",
  "Bivalent TSS",
  "Heterochromatin",
  "ZNF / repeats",
  "PMD",
  "H3K27me3",
  "H2AK119",
  "JARID2",
  "CTCF"
)

# Keep only categories that actually occur, but keep ALL X-axis targets
present_categories <- kycg_order[kycg_order %in% unique(kycg_plot$category)]

kycg_plot <- kycg_plot %>%
  mutate(
    category = factor(category, levels = rev(present_categories)),
    fill_state = factor(
      if_else(best_global_FDR < 1e-15, "filled", "open"),
      levels = c("filled", "open")
    )
  )

highlight_df <- bind_rows(
  tibble(
    group = "Promoter",
    category = c("CpG island", "Active TSS", "H3K4me3"),
    xmin = 0.62, xmax = 1.38
  ),
  tibble(
    group = "ReprPC",
    category = c("CGI shore", "ReprPC"),
    xmin = 1.55, xmax = 4.45
  ),
  tibble(
    group = "Bivalent",
    category = "Bivalent TSS",
    xmin = 7.55, xmax = 10.45
  )
) %>%
  filter(category %in% as.character(kycg_plot$category)) %>%
  mutate(y = match(category, levels(kycg_plot$category)))

# Plot -------------------------------------------------------------------------
p1_band <- highlight_df %>%
  filter(group == "Promoter") %>%
  summarise(
    xmin = min(xmin),
    xmax = max(xmax),
    ymin = min(y) - 0.43,
    ymax = max(y) + 0.43
  )

p_kycg <- ggplot(kycg_plot, aes(x = target, y = category)) +
  geom_point(
    aes(size = max_overlap, color = score, shape = fill_state),
    stroke = 0.4
  ) +
  geom_rect(
    data = p1_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#D88900", alpha = 0.13, colour = NA
  ) +
  geom_rect(
    data = filter(highlight_df, group != "Promoter"),
    aes(xmin = xmin, xmax = xmax, ymin = y - 0.43,
        ymax = y + 0.43, fill = group),
    inherit.aes = FALSE,
    alpha = 0.12, colour = NA
  ) +
  scale_fill_manual(
    values = c("ReprPC" = "#292473", "Bivalent" = "#C2185B"),
    guide = "none"
  ) +
  scale_size_area(
    max_size = 3.7,
    name = "Overlap",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      nrow = 1,
      order = 1
    )
  ) +
  scale_shape_manual(
    name = "FDR",
    values = c(filled = 16, open = 1),
    labels = c("< 1e-15", "< 0.05"),
    guide = guide_legend(
      order = 3,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      nrow = 1,
      override.aes = list(size = 2.5, color = "grey35")
    )
  ) +
  scale_color_gradient(
    low = kycg_low,
    high = kycg_high,
    name = expression(-log[10]("FDR")),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(13, "mm"),
      barheight = unit(1.5, "mm"),
      ticks = FALSE,
      order = 2
    )
  ) +
  scale_x_discrete(
    limits = target_order,
    labels = target_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(present_categories),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    tag = "g"
  ) +
  basic +
  theme(
    axis.line = element_blank(),
    # axis.ticks = element_blank(),
    axis.text.x = ggtext::element_markdown(
      size = 5.5,
      lineheight = 0.85,
      angle = 0,
      hjust = 0.5,
      vjust = 1
    ),
    axis.text.y = element_text(size = 5.5),
    
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.4
    ),
    plot.margin = margin(l = 10, r = 5, b = 2, t = 5),
    plot.tag.position = c(0.00, 0.95),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "right",
    legend.title = element_text(size = 5.5, margin = margin(b = 1)),
    legend.text = element_text(size = 5, margin = margin(t = 2)),
    legend.spacing.x = unit(5, "mm"),
    legend.spacing.y = unit(0, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.spacing.x = unit(2, "mm"),
    legend.box.spacing =  unit(1, "mm")
  )

ggsave("result/kycg.svg", device = svglite, p_kycg, width = 3.8, height = 2, dpi = 800)

# Section: Cell-type context of programs (original 05_type.r)
source(file.path(repo, "r", "figures", "figure_environment.R"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(scales)
})

# PATH + LOAD

if (!exists("find_path")) {
  find_path <- function(target, max_up = 3L) {
    p <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    for (i in 0:max_up) {
      if (basename(p) == target) return(p)
      candidate <- file.path(p, target)
      if (dir.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
      p2 <- dirname(p)
      if (identical(p2, p)) break
      p <- p2
    }
    stop("Could not find directory: ", target, call. = FALSE)
  }
}

base <- Sys.getenv("LUNG_METHYLATION_PROJECT_ROOT", unset = "")
if (!nzchar(base)) base <- find_path("Project", 3)
setwd(file.path(base, "05_ncc"))

project_dir <- file.path(base, "05_ncc")
cache_dir <- file.path(project_dir, "object", "cache")

cloud_df <- read_csv(file.path(cache_dir, "PROGRAM_CPG_CLOUD.csv.gz"), show_col_types = FALSE)
thresholds <- read_csv(file.path(cache_dir, "CPG_BONFERRONI_THRESHOLDS.csv"), show_col_types = FALSE)
program_modes <- read_csv(file.path(cache_dir, "PROGRAM_CLOUD_MODES.csv"), show_col_types = FALSE)
sample_pc <- read_csv(file.path(cache_dir, "PROGRAM_BIPLOT_SAMPLES.csv"), show_col_types = FALSE)
centroids <- read_csv(file.path(cache_dir, "PROGRAM_BIPLOT_CENTROIDS.csv"), show_col_types = FALSE)
loadings <- read_csv(file.path(cache_dir, "PROGRAM_BIPLOT_LOADINGS.csv"), show_col_types = FALSE)
pca_variance <- read_csv(file.path(cache_dir, "PROGRAM_BIPLOT_VARIANCE.csv"), show_col_types = FALSE)

fine_order <- c("Neu", "Mono", "Eos", "Bas", "Bnv", "Bmem", "CD4nv", "CD4mem", "CD8nv", "CD8mem", "Treg", "NK")

fine_labels <- expression(
  Neu    = "Neutrophil",
  Mono   = "Monocyte",
  Eos    = "Eosinophil",
  Bas    = "Basophil",
  Bnv    = "B naive",
  Bmem   = "B memory",
  CD4nv  = "CD4 naive",
  CD4mem = "CD4 memory",
  CD8nv  = "CD8 naive",
  CD8mem = "CD8 memory",
  Treg   = T[reg],
  NK     = "NK cell"
)

myeloid_types <- c("Neu", "Mono", "Eos", "Bas")
naive_nk_types <- c("Bnv", "CD4nv", "CD8nv", "NK")
memory_treg_types <- c("Bmem", "CD4mem", "CD8mem", "Treg")

# Default black; only the two program pairs are highlighted.
program_color <- function(x) {
  case_when(
    x %in% c("P3", "P6") ~ "#292473",
    x %in% c("P5", "P8") ~ "#C2185B",
    TRUE ~ "#000000"
  )
}

# 1. OVERALL CELL-TYPE SENSITIVITY CLOUD

cloud_trans <- pseudo_log_trans(sigma = 1, base = 10)

all_cpg <- cloud_df %>%
  select(base_cpg, broad, fine) %>%
  distinct() %>%
  filter(is.finite(broad), is.finite(fine))

bonf_x <- thresholds %>% filter(gse == "GSE110554") %>% pull(bonf_threshold)
bonf_y <- thresholds %>% filter(gse == "GSE167998") %>% pull(bonf_threshold)

x_t <- cloud_trans$transform(bonf_x)
y_t <- cloud_trans$transform(bonf_y)
x_rng_t <- cloud_trans$transform(c(0, 60))
y_rng_t <- cloud_trans$transform(c(0, 60))

dx <- diff(x_rng_t) * 0.05
dy <- diff(y_rng_t) * 0.040

quad_labels <- tibble(
  label = c("Q1", "Q2", "Q3", "Q4"),
  x = cloud_trans$inverse(c(x_t + dx, x_t - dx, x_t - dx, x_t + dx)),
  y = cloud_trans$inverse(c(y_t + dy, y_t + dy, y_t - dy, y_t - dy))
)

label_layout <- tribble(
  ~program_label, ~dx,    ~dy,     ~hjust,
  "P1",            0.020,  0.000,   0,
  "P2",           -0.020, -0.006,   1,
  "P3",            0.000, -0.050,   0.5,   # point 바로 아래
  "P4",           -0.016,  0.016,   1,
  "P5",           -0.016, -0.018,   1,
  "P6",            0.005,  0.050,   0.5,   # point 바로 위
  "P7",            0.018,  0.006,   0,
  "P8",           -0.018, -0.014,   1,
  "P9",            0.018,  0.006,   0
)

program_plot <- program_modes %>%
  left_join(label_layout, by = "program_label") %>%
  mutate(
    dx = replace_na(dx, 0.01),
    dy = replace_na(dy, 0),
    hjust = replace_na(hjust, 0),
    point_color = case_when(
      program_label %in% c("P3", "P6") ~ "#292473",
      program_label %in% c("P5", "P8") ~ "#C2185B",
      TRUE ~ "#000000"
    ),
    x_plot = cloud_trans$transform(mode_broad),
    y_plot = cloud_trans$transform(mode_fine),
    label_x = cloud_trans$inverse(x_plot + dx * diff(x_rng_t)),
    label_y = cloud_trans$inverse(y_plot + dy * diff(y_rng_t))
  )

ann_x_t <- x_t + diff(x_rng_t) * 0.15
ann_y_t <- y_rng_t[1] + diff(y_rng_t) * 0.045
ann_x <- cloud_trans$inverse(ann_x_t)
ann_y <- cloud_trans$inverse(ann_y_t)

p_cloud <- ggplot(all_cpg, aes(broad, fine)) +
  geom_vline(
    xintercept = bonf_x,
    linewidth = 0.2,
    color = "gray10"
  ) +
  geom_hline(
    yintercept = bonf_y,
    linewidth = 0.2,
    color = "gray10"
  ) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE,
    n = 250,
    alpha = 0.9,
    interpolate = TRUE
  ) +
  scale_x_continuous(
    trans = cloud_trans,
    breaks = c(0, 2, 5, 10, 20, 50),
    limits = c(0, 60),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans = cloud_trans,
    breaks = c(0, 2, 5, 10, 20, 50),
    limits = c(0, 60),
    expand = c(0, 0)
  ) +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"),
    name = "CpG density",
    breaks = function(x) range(x, na.rm = TRUE),
    labels = c("Low", "High"),
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "left",
      label.position = "right",
      barwidth = unit(0.07, "in"),
      barheight = unit(0.55, "in"),
      ticks = FALSE,
      frame.colour = "black"
    )
  ) +

  # geom_text(
  #   data = quad_labels,
  #   aes(x = x, y = y, label = label),
  #   inherit.aes = FALSE,
  #   size = 2.2,
  #   fontface = "bold",
  #   color = "black"
  # ) +
  geom_point(
    data = program_plot,
    aes(x = mode_broad, y = mode_fine, color = point_color),
    inherit.aes = FALSE,
    size = 1,
    shape = 20,
    stroke = 0.3
  ) +
  geom_text(
    data = program_plot,
    aes(x = label_x, y = label_y, label = program_label, color = point_color, hjust = hjust),
    inherit.aes = FALSE,
    size = 2,
    vjust = 0.5
  ) +
  scale_color_identity() +
  annotate(
    "text",
    x = ann_x,
    y = ann_y,
    label = paste0(
      "Q1  Both levels\n",
      "Q2  Subtype-specific\n",
      "Q3  Weak / none\n",
      "Q4  Major-type specific"
    ),
    hjust = 0,
    vjust = 0,
    size = 1.9,
    lineheight = 1.25,
    color = "black"
  ) +
  labs(
    x = "Major cell-type association",
    y = "Leukocyte-subtype association",
    tag = "c"
  ) +
  basic +
  theme(
    plot.tag.position = c(0, 0.972),
    legend.position = "right",
    legend.justification = c(0.5, 1),
    legend.box.just = "top",
    legend.direction = "vertical",
    
    legend.title = element_text(
      size = 6,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(r = 4, unit = "pt")
    ),
    
    legend.text = element_text(size = 5.5, margin = margin(l = 2, unit = "pt")),
    legend.margin = margin(t = 5, r = 0, b = 0, l = -10, unit = "pt")
  )

ggsave(
  "result/cloud.svg",
  p_cloud,
  device = svglite,
  width = 2.8,
  height = 2.3,
  units = "in",
  dpi = 800,
  bg = "white"
)

# 2. RELATIVE PROGRAM-PROFILE PCA BIPLOT

cell_type_colors <- c(
  Neu    = "#8AC0C1",
  Mono   = "#76B2B3",
  Eos    = "#5EA2A4",
  Bas    = "#489295",
  
  Bnv    = "#8EA5AD",
  CD4nv  = "#7C97A1",
  CD8nv  = "#6C8B96",
  
  Bmem   = "#DF7C6F",
  CD4mem = "#D86C5F",
  CD8mem = "#C95B4E",
  
  Treg   = "#CFB686",
  NK     = "#A69AB4"
)

sample_pc <- sample_pc %>%
  mutate(
    type = factor(type, levels = fine_order),
    marginal_x = if_else(as.character(type) %in% myeloid_types, "Myeloid", "Others"),
    marginal_y = if_else(as.character(type) %in% c("Bmem", "CD4mem", "CD8mem"), "Memory", "Others")
  )

# Direct side-curve labels; no marginal-density legend is needed.
x_density_labels <- sample_pc %>%
  group_by(marginal_x) %>%
  summarise(PC2 = median(PC2, na.rm = TRUE), .groups = "drop") %>%
  mutate(density_y = 0.2)

y_density_labels <- sample_pc %>%
  group_by(marginal_y) %>%
  summarise(PC1 = median(PC1, na.rm = TRUE), .groups = "drop") %>%
  mutate(density_x = 0.2)

# Density accents are softer than the point palette.
density_myeloid <- "#76B2B3"  # 청록
density_memory <- "#D86C5F"  # 블루그레이
density_others  <- "#62676B"  # 회색 점선 유지

# Hull에는 NK와 Treg 제외
hull_pc <- sample_pc %>%
  mutate(hull_group = case_when(
    type %in% myeloid_types ~ "Myeloid",
    type %in% c("Bnv", "CD4nv", "CD8nv") ~ "Naive",
    type %in% c("Bmem", "CD4mem", "CD8mem") ~ "Memory",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(hull_group))

# Hull label positions after swapping PC axes
pc1_span <- diff(range(sample_pc$PC1, na.rm = TRUE))
hull_labels <- hull_pc %>%
  group_by(hull_group) %>%
  summarise(
    x = median(PC2),
    y = case_when(
      first(hull_group) == "Myeloid" ~ max(PC1),
      first(hull_group) == "Naive" ~ min(PC1),
      first(hull_group) == "Memory" ~ max(PC1)
    ),
    .groups = "drop"
  )

# Loading arrow 길이
arrow_scale <- 0.62
label_gap <- 0.0025

loadings <- loadings %>%
  mutate(
    PC1_arrow = PC1_plot * arrow_scale,
    PC2_arrow = PC2_plot * arrow_scale,
    PC2_arrow = case_when(
      program_label == "P2" ~ PC2_plot * arrow_scale + 0.0001,
      TRUE ~ PC2_plot * arrow_scale
    ),
    arrow_r = sqrt(PC1_arrow^2 + PC2_arrow^2),
    PC1_label = PC1_arrow + label_gap * PC1_arrow / arrow_r,
    PC2_label = PC2_arrow + label_gap * PC2_arrow / arrow_r,
    PC1_label = case_when(
      program_label == "P9" ~ PC1_label - 0.001,
      program_label == "P3" ~ PC1_label + 0.0001,
      TRUE ~ PC1_label
    ),
    PC2_label = case_when(
      program_label == "P2" ~ PC2_label + 0.0005,
      TRUE ~ PC2_label
    )
  )

load_black <- loadings %>% filter(!program_label %in% c("P3", "P5", "P6", "P8"))
load_36 <- loadings %>% filter(program_label %in% c("P3", "P6"))
load_58 <- loadings %>% filter(program_label %in% c("P5", "P8"))

x_lab <- sprintf("PC2 (%.1f%%)", pca_variance$PC2[1])
y_lab <- bquote(PC1~"("*bold(.(sprintf("%.1f%%", pca_variance$PC1[1])))*")")

# 0 근처는 벌리고 큰 절댓값은 압축
pc_trans <- scales::pseudo_log_trans(sigma = 0.01, base = 10)

biplot <- ggplot(sample_pc, aes(PC2, PC1)) +
  geom_hline(yintercept = 0, linewidth = 0.1, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.1, color = "gray70") +
  
  
  geom_point(
    data = sample_pc %>% filter(type %in% myeloid_types),
    aes(color = type),
    size = 1, alpha = 0.85, shape = 16
  ) +
  scale_color_manual(
    values = cell_type_colors[myeloid_types],
    breaks = myeloid_types,
    labels = fine_labels[myeloid_types],
    name = "Myeloid",
    guide = guide_legend(
      order = 1,
      override.aes = list(size = 1)
    )
  ) +
  
  ggnewscale::new_scale_color() +
  
  geom_point(
    data = sample_pc %>% filter(type %in% naive_nk_types),
    aes(color = type),
    size = 1, alpha = 0.85, shape = 16
  ) +
  scale_color_manual(
    values = cell_type_colors[naive_nk_types],
    breaks = naive_nk_types,
    labels = fine_labels[naive_nk_types],
    name = "Naive + NK",
    guide = guide_legend(
      order = 2,
      override.aes = list(size = 1)
    )
  ) +
  
  ggnewscale::new_scale_color() +
  
  geom_point(
    data = sample_pc %>% filter(type %in% memory_treg_types),
    aes(color = type),
    size = 1, alpha = 0.85, shape = 16
  ) +
  scale_color_manual(
    values = cell_type_colors[memory_treg_types],
    breaks = memory_treg_types,
    labels = fine_labels[memory_treg_types],
    name = expression("Memory + " * T[reg]),
    guide = guide_legend(
      order = 3,
      override.aes = list(size = 1)
    )
  ) +
  
  geom_segment(
    data = load_black,
    aes(x = 0, y = 0, xend = PC2_arrow, yend = PC1_arrow),
    inherit.aes = FALSE,
    color = "black",
    # color = "gray40", 
    alpha = 0.85,
    linewidth = 0.2,
    arrow = arrow(length = grid::unit(0.065, "cm"))
  ) +
  geom_text(
    data = load_black,
    aes(PC2_label, PC1_label, label = program_label),
    inherit.aes = FALSE,
    # color = "gray40",
    alpha = 0.85,
    color = "black",
    # fontface = "bold", 
    size = 1.5
  ) +
  
  geom_segment(
    data = load_36,
    aes(x = 0, y = 0, xend = PC2_arrow, yend = PC1_arrow),
    inherit.aes = FALSE,
    color = "black", 
    alpha = 0.85,
    # color = "#292473", 
    linewidth = 0.2,
    arrow = arrow(length = grid::unit(0.065, "cm"))
  ) +
  geom_text(
    data = load_36,
    aes(PC2_label, PC1_label, label = program_label),
    inherit.aes = FALSE,
    color = "black", 
    alpha = 0.85,
    # color = "#292473", 
    # fontface = "bold",
    size = 1.5
  ) +
  
  geom_segment(
    data = load_58,
    aes(x = 0, y = 0, xend = PC2_arrow, yend = PC1_arrow),
    inherit.aes = FALSE,
    color = "black", 
    alpha = 0.85,
    # color = "#C2185B",
    linewidth = 0.2,
    arrow = arrow(length = grid::unit(0.065, "cm"))
  ) +
  geom_text(
    data = load_58,
    aes(PC2_label, PC1_label, label = program_label),
    inherit.aes = FALSE,
    color = "black", 
    alpha = 0.85,
    # color = "#C2185B",
    # fontface = "bold",
    size = 1.5
  ) +
  
  # Top margin (PC2): Myeloid filled/solid vs Others unfilled/dashed.
  # Use group-specific fixed aesthetics to avoid adding a density legend.
  ggside::geom_xsidedensity(
    data = filter(sample_pc, marginal_x == "Myeloid"),
    aes(x = PC2, y = after_stat(scaled)),
    inherit.aes = FALSE, adjust = 1.35, position = "identity",
    fill = "#8AC0C1", colour = density_myeloid,
    alpha = 0.55, linewidth = 0.25, linetype = "solid",
    outline.type = "upper", show.legend = FALSE
  ) +
  ggside::geom_xsidedensity(
    data = filter(sample_pc, marginal_x == "Others"),
    aes(x = PC2, y = after_stat(scaled)),
    inherit.aes = FALSE, adjust = 1.35, position = "identity",
    fill = NA, colour = density_others,
    linewidth = 0.25, linetype = "dashed",
    outline.type = "upper", show.legend = FALSE
  ) +
  # Right margin (PC1): Memory filled/solid vs Others unfilled/dashed.
  ggside::geom_ysidedensity(
    data = filter(sample_pc, marginal_y == "Memory"),
    aes(x = after_stat(scaled), y = PC1),
    inherit.aes = FALSE, adjust = 1.35, position = "identity",
    fill = "#DF7C6F", colour = density_memory,
    alpha = 0.55, linewidth = 0.25, linetype = "solid",
    outline.type = "upper", show.legend = FALSE
  ) +
  ggside::geom_ysidedensity(
    data = filter(sample_pc, marginal_y == "Others"),
    aes(x = after_stat(scaled), y = PC1),
    inherit.aes = FALSE, adjust = 1.35, position = "identity",
    fill = NA, colour = density_others,
    linewidth = 0.25, linetype = "dashed",
    outline.type = "upper", show.legend = FALSE
  ) +
  # Direct labels on the side panels, matching visualization.R's original style.
  ggside::geom_xsidetext(
    data = filter(x_density_labels, marginal_x == "Myeloid"),
    aes(x = PC2, y = density_y, label = marginal_x),
    inherit.aes = FALSE, colour = density_myeloid,
    size = 1.6, 
    fontface = "bold", 
    show.legend = FALSE
  ) +
  ggside::geom_xsidetext(
    data = filter(x_density_labels, marginal_x == "Others"),
    aes(x = PC2, y = density_y, label = marginal_x),
    inherit.aes = FALSE, colour = density_others,
    size = 1.6,
    fontface = "bold", 
    show.legend = FALSE
  ) +
  ggside::geom_ysidetext(
    data = filter(y_density_labels, marginal_y == "Memory"),
    aes(x = density_x, y = PC1, label = marginal_y),
    inherit.aes = FALSE, colour = density_memory,
    size = 1.6,
    fontface = "bold",
    angle = 270,
    show.legend = FALSE
  ) +
  ggside::geom_ysidetext(
    data = filter(y_density_labels, marginal_y == "Others"),
    aes(x = density_x, y = PC1, label = marginal_y),
    inherit.aes = FALSE, colour = density_others,
    size = 1.6, 
    fontface = "bold", 
    angle = 270, show.legend = FALSE
  ) +
  ggside::ggside(x.pos = "top", y.pos = "right") +
  ggside::scale_xsidey_continuous(
    breaks = 1, labels = "max", limits = c(0, 1.07),
    minor_breaks = NULL, expand = expansion(mult = c(0.03, 0.03))
  ) +
  ggside::scale_ysidex_continuous(
    breaks = 1, labels = "max", limits = c(0, 1.07),
    minor_breaks = NULL, expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_x_continuous(
    trans = pc_trans,
    breaks = c(-0.06, -0.03, -0.01, 0, 0.01, 0.03, 0.06),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    trans = pc_trans,
    breaks = c(-0.06, -0.03, -0.01, 0, 0.01, 0.03, 0.06),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  
  labs(x = x_lab, y = y_lab, tag = "h") +
  basic +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = unit(3, "mm"),
    legend.spacing.y = unit(2, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 6, margin = margin(b = 0.5)),
    ggside.panel.scale.x = 0.2,
    ggside.panel.scale.y = 0.2,
    ggside.panel.spacing = unit(1, "mm"),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_text(size = 6),
    ggside.axis.line = element_line(linewidth = 0.35),
    ggside.axis.ticks = element_line(linewidth = 0.35)
  )

ggsave(
  "result/biplot.svg",
  device = svglite,
  biplot,
  width = 3.1, height = 2,
  units = "in", dpi = 800, bg = "transparent"
)

# Section: Program embedding and biplots (original 06_embed.r)
# UMAP + main panel + 9 program panels
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(uwot)
library(ggnewscale)
library(patchwork)
library(tidyverse)
library(ggrepel)
library(scales)

# ── 1. Load data ──────────────────────────────────────────────────────────────

score <- read_csv(
  file.path(build_dir, "cache", "model", "SOURCE_PROGRAM_SCORES.csv"),
  show_col_types = FALSE
)

post <- read_csv(
  file.path(build_dir, "cache", "model", "SOURCE_POSTERIOR.csv"),
  show_col_types = FALSE
)

dat <- score %>%
  left_join(
    post %>%
      select(sentrix_id, group, time, type, hard_state, TM),
    by = "sentrix_id"
  )

program_cols <- grep("^program_", names(dat), value = TRUE)

# program_1, ..., program_9 순서 보장
program_cols <- program_cols[
  order(as.integer(stringr::str_extract(program_cols, "\\d+")))
]

X <- dat %>%
  select(all_of(program_cols)) %>%
  as.matrix()

stopifnot(
  ncol(X) == 9,
  all(is.finite(X))
)


# ── 2. UMAP ───────────────────────────────────────────────────────────────────

um <- uwot::umap(
  X, seed = 42,
  n_neighbors = 9,
  min_dist = 0.18,
  metric = "euclidean",
  n_components = 5,
  init = "spectral",
  spread = 1,
  repulsion_strength = 1,
  learning_rate = 1,
  n_epochs = 1000,
  ret_model = FALSE,
  verbose = TRUE,
  n_threads = 1
)

saveRDS(um, "object/umap.rds")

# ── 3. Main panel: control vs prediagnostic time gradient ────────────────────
#
# prediagnostic = pale red-ish (far) -> vivid red (near)

emb <- dat %>%
  mutate(
    UMAP1 = um[, 1],
    UMAP2 = um[, 2],
    time_to_dx = if_else(group == "prediagnostic", as.numeric(time), NA_real_),
    time_plot = pmin(time_to_dx, 10)
  )

size <- 0.75

p_ttd <- ggplot() +
  geom_point(
    data = emb %>% filter(group == "control"),
    aes(UMAP1, UMAP2),
    size = size, shape = 16, color = "#D3D3D3"
  ) +
  geom_point(
    data = emb %>% filter(group == "prediagnostic") %>% arrange(desc(time_plot)),
    aes(UMAP1, UMAP2, color = time_plot),
    size = size, shape = 16
  ) +
  scale_color_gradientn(
    colours = c("#B2182B", "#D73027", "#E86A61", "#F3AAA8", "#F6D5D6"),
    values = scales::rescale(c(0, 0.75, 2, 5, 10)),
    limits = c(0, 10),
    breaks = c(0, 2, 5, 10),
    labels = c("0", "2", "5", "\u226510"),
    oob = scales::squish,
    name = "Time-to-diagnosis",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      barwidth = unit(0.57, "in"),
      barheight = unit(0.045, "in"),
      ticks = FALSE,
      frame.colour = NA
    )
  ) +
  labs(tag = "d") +
  theme_void(base_size = 9) +
  theme(
    plot.margin = margin(l = 5, r = 5, b = 10, t = 5),
    plot.tag = element_text(
      family = "Arial", face = "bold", size = 14,
      margin = margin(l = -5)
    ),
    plot.tag.location = "plot",
    plot.tag.position = c(0.05, 0.95),
    
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.05),
    legend.justification = c(1, 0),
    legend.direction = "horizontal",
    
    legend.title = element_text(
      size = 6, hjust = 0.5,
      margin = margin(b = 2.5, unit = "pt")
    ),
    legend.text = element_text(
      size = 5.5,
      margin = margin(t = 2, unit = "pt")
    ),
    legend.margin = margin(b = -13, unit = "pt")
  )

ggsave(
  "result/umap-ttd.svg",
  p_ttd,
  device = svglite,
  width = 1.35,
  height = 1.2,
  units = "in",
  dpi = 800,
  bg = "white"
)

# ── Histology UMAP ────────────────────────────────────────────────────────────

emb_type <- emb %>%
  mutate(
    histology = case_when(
      group == "control" ~ "Normal",
      group == "prediagnostic" & !is.na(type) & type != "nan" ~ type,
      TRUE ~ "Unknown"
    ),
    histology = factor(
      histology,
      levels = c("Normal", "LUAD", "LUSC", "NSCLC", "SCLC", "Other", "Unknown")
    ),
    proximity_alpha = case_when(
      group == "prediagnostic" ~
        0.05 + 0.92 * (1 - pmin(time_to_dx, 10) / 10)^1.5,
      TRUE ~ NA_real_
    )
  )

p_type <- ggplot() +
  
  # Normal: background
  geom_point(
    data = emb_type %>% filter(group == "control"),
    aes(UMAP1, UMAP2, color = histology),
    size = size,
    shape = 16,
    alpha = 0.25
  ) +
  
  # Prediagnostic: far first, near last
  geom_point(
    data = emb_type %>%
      filter(group == "prediagnostic") %>%
      arrange(desc(time_to_dx)),
    aes(
      UMAP1,
      UMAP2,
      color = histology,
      alpha = proximity_alpha
    ),
    size = size,
    shape = 16
  ) +
  
  scale_color_manual(
    values = c(
      "Normal"  = "#D9DCE2",
      "LUAD"    = "#eb6f5d",
      "LUSC"    = "#71c9dd",
      "NSCLC"   = "#6376a0",
      "SCLC"    = "#33b39f",
      "Other"   = "#9da7c3",
      "Unknown" = "#a7dace"
    ),
    drop = TRUE,
    name = NULL
  ) +
  
  scale_alpha_identity(
    guide = "none"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 1, alpha = 1)
    )
  ) +
  
  theme_void(base_size = 9) +
  theme(
    plot.margin = margin(l = 5, r = 5, b = 10, t = 5),
    
    legend.position = "right",
    legend.justification = c(0, 0),
    legend.box.just = "bottom",
    legend.title.align = 0,
    legend.spacing = unit(2, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.text = element_text(
      size = 5.5,
      margin = margin(l = 1)
    ),
    legend.box.spacing = unit(2, "mm"),
    legend.margin = margin(l = -18, 0, b = -9, 0)
  )

# ggsave(
#   "result/umap-type.png",
#   p_type,
#   width = 1.6,
#   height = 1.2,
#   units = "in",
#   dpi = 800,
#   bg = "white"
# )

# ── Program data ──────────────────────────────────────────────────────────────

prog_long <- emb %>%
  select(sentrix_id, group, time_to_dx, UMAP1, UMAP2, all_of(program_cols)) %>%
  pivot_longer(
    cols = all_of(program_cols),
    names_to = "program",
    values_to = "value"
  ) %>%
  mutate(
    program_num = as.integer(str_extract(program, "\\d+")),
    program_lab = paste0("P", program_num)
  ) %>%
  group_by(program_lab) %>%
  mutate(
    value01 = scales::rescale(value, to = c(0, 1))
  ) %>%
  ungroup()

# ── Program-specific color ────────────────────────────────────────────────────

program_color <- function(x) {
  case_when(
    x %in% c("P3", "P6") ~ "#292473",
    x %in% c("P5", "P8") ~ "#C2185B",
    TRUE ~ "#000000"
  )
}

# ── Program panel ─────────────────────────────────────────────────────────────

make_program_panel <- function(
    pg_label,
    point_size = 0.25,
    low_col = "#F2F2F2",
    high_col = "#4B4F8C"
) {
  
  d <- prog_long %>%
    filter(program_lab == pg_label) %>%
    arrange(value01)
  
  ggplot(d, aes(UMAP1, UMAP2, color = value01)) +
    geom_point(
      shape = 16,
      size = point_size
    ) +
    scale_color_gradient(
      low = low_col,
      high = high_col,
      limits = c(0, 1)
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = pg_label,
      hjust = -0.15,
      vjust = 1.15,
      family = "Arial",
      fontface = "bold",
      size = 2.3,
      colour = program_color(pg_label)
    ) +
    theme_void(base_size = 9) +
    theme(
      legend.position = "none",
      plot.margin = margin(1, 1, 1, 1, unit = "pt")
    )
}


# ── Program panels ────────────────────────────────────────────────────────────

p_list <- setNames(
  lapply(paste0("P", 1:9), make_program_panel),
  paste0("P", 1:9)
)

# ── Navigation panel ──────────────────────────────────────────────────────────

x0 <- 0.22
y0 <- 0.22
ov <- 0.035

p_nav <- ggplot() +
  annotate(
    "segment",
    x = x0 - ov, y = y0,
    xend = 0.88, yend = y0,
    linewidth = 0.25,
    lineend = "square",
    arrow = arrow(
      type = "closed",
      angle = 20,
      length = unit(0.045, "in")
    )
  ) +
  annotate(
    "segment",
    x = x0, y = y0 - ov,
    xend = x0, yend = 0.88,
    linewidth = 0.25,
    lineend = "square",
    arrow = arrow(
      type = "closed",
      angle = 20,
      length = unit(0.045, "in")
    )
  ) +
  annotate(
    "text",
    x = 0.57, y = 0.10,
    label = "UMAP 1",
    family = "Arial",
    size = 2.1
  ) +
  annotate(
    "text",
    x = 0.10, y = 0.57,
    label = "UMAP 2",
    family = "Arial",
    angle = 90,
    size = 2.1
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    plot.margin = margin(1, 1, 1, 1, unit = "pt")
  )

# ── Layout: P1-P5 top, P6-P9 bottom, left aligned ─────────────────────────────

design_txt <- "
ABCDE
NFGHI
"

p_face <- wrap_plots(
  A = p_list$P1,
  B = p_list$P2,
  C = p_list$P3,
  D = p_list$P4,
  E = p_list$P5,
  N = p_nav,
  F = p_list$P6,
  G = p_list$P7,
  H = p_list$P8,
  I = p_list$P9,
  design = design_txt
)

ggsave(
  "result/umap-prog.svg",
  device = svglite,
  p_face,
  width = 2.85,
  height = 1.2,
  units = "in",
  dpi = 800,
  bg = "white"
)

# PROGRAM TOP100 CpG OVERLAP: UpSet plot

programs <- paste0("P", 1:9)

loading_file <- file.path("object", "cache", "PROGRAM_LOADING_BASE.csv.gz")
fallback_file <- file.path("result", "PROGRAM_TOP100_FROM_REBUILD.csv")

if (file.exists(loading_file)) {
  loading <- read_csv(loading_file, show_col_types = FALSE)
} else if (file.exists(fallback_file)) {
  loading <- read_csv(fallback_file, show_col_types = FALSE)
} else {
  stop("PROGRAM_LOADING_BASE.csv.gz 파일을 찾을 수 없습니다.")
}

stopifnot(all(c("program_label", "base_cpg", "loading") %in% names(loading)))

# 2. Select top100 CpGs --------------------------------------------------------

top100 <- loading %>%
  filter(
    program_label %in% programs,
    !is.na(base_cpg),
    is.finite(loading)
  ) %>%
  group_by(program_label, base_cpg) %>%
  summarise(loading = sum(loading), .groups = "drop") %>%
  group_by(program_label) %>%
  slice_max(loading, n = 100, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(program_label = factor(program_label, levels = programs))

stopifnot(
  n_distinct(top100$program_label) == 9,
  all(as.integer(table(top100$program_label)) == 100)
)

# 3. Exact intersections -------------------------------------------------------

membership <- top100 %>%
  arrange(base_cpg, program_label) %>%
  group_by(base_cpg) %>%
  summarise(
    combination = paste(as.character(program_label), collapse = "+"),
    n_sets = n_distinct(program_label),
    .groups = "drop"
  )

exact <- membership %>%
  count(combination, n_sets, name = "n_cpg") %>%
  mutate(
    pair_39 = str_detect(combination, "(^|\\+)P3(\\+|$)") &
      str_detect(combination, "(^|\\+)P9(\\+|$)"),
    pair_58 = str_detect(combination, "(^|\\+)P5(\\+|$)") &
      str_detect(combination, "(^|\\+)P8(\\+|$)"),
    category = case_when(
      pair_39 ~ "P3 + P9",
      pair_58 ~ "P5 + P8",
      TRUE ~ "Other"
    )
  )

# 4. Order ALL columns by descending intersection size -------------------------

# 공유 조합은 최소 2개 CpG 이상
# 단독 프로그램은 P1~P9 모두 포함
min_exact_cpg <- 2

plot_df <- exact %>%
  filter(n_sets == 1 | (n_sets >= 2 & n_cpg >= min_exact_cpg)) %>%
  arrange(desc(n_cpg), n_sets, combination)

comb_order <- plot_df$combination

plot_df <- plot_df %>%
  mutate(combination = factor(combination, levels = comb_order))

# 5. Membership matrix ---------------------------------------------------------

active <- top100 %>%
  inner_join(
    membership %>% select(base_cpg, combination),
    by = "base_cpg"
  ) %>%
  distinct(combination, program_label) %>%
  transmute(
    combination,
    program = as.character(program_label),
    active = TRUE
  )

mat <- expand_grid(
  combination = comb_order,
  program = programs
) %>%
  left_join(active, by = c("combination", "program")) %>%
  left_join(
    plot_df %>%
      mutate(combination = as.character(combination)) %>%
      select(combination, category),
    by = "combination"
  ) %>%
  mutate(
    active = replace_na(active, FALSE),
    y = match(program, rev(programs)),
    combination = factor(combination, levels = comb_order)
  )

segs <- mat %>%
  filter(active) %>%
  group_by(combination, category) %>%
  summarise(
    ymin = min(y),
    ymax = max(y),
    .groups = "drop"
  ) %>%
  filter(ymin != ymax)

# 6. Colors --------------------------------------------------------------------

pal <- c(
  "P3 + P9" = "#292473",
  "P5 + P8" = "#C2185B",
  "Other"   = "#8CA9B9"
)

# 7. Upper panel ---------------------------------------------------------------

p_bar <- ggplot(plot_df, aes(x = combination, y = n_cpg)) +
  geom_col(
    aes(fill = category),
    width = 0.72
  ) +
  geom_text(
    aes(label = n_cpg),
    vjust = -0.5,
    size = 2
  ) +
  scale_fill_manual(
    values = pal,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = comb_order,
    drop = FALSE,
    expand = expansion(add = 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, max(plot_df$n_cpg) * 1.20),
    breaks = pretty(c(0, max(plot_df$n_cpg)), n = 3),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Allocated CpGs", tag = "e"
  ) +
  basic +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# 8. Lower panel ---------------------------------------------------------------

p_matrix <- ggplot(mat, aes(x = combination, y = y)) +
  geom_point(
    color = "#DFE3E6",
    size = 1
  ) +
  geom_segment(
    data = segs,
    aes(
      x = combination,
      xend = combination,
      y = ymin,
      yend = ymax,
      color = category
    ),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_point(
    data = filter(mat, active),
    aes(color = category),
    size = 1.3
  ) +
  scale_color_manual(
    values = pal,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = comb_order,
    drop = FALSE,
    expand = expansion(add = 0.5)
  ) +
  scale_y_continuous(
    breaks = 1:9,
    labels = rev(programs),
    limits = c(0.5, 9.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  basic +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.position = "none"
  )

# 9. Combine and save

p_upset <- (p_bar / p_matrix) +
  patchwork::plot_layout(heights = c(1.4, 1.2)) & 
  theme(
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5)
  )

ggsave("result/upset.svg", device = svglite, p_upset, width = 4.2, height = 2, dpi = 800)

# ── Load matched program scores ───────────────────────────────────────────────

score_file <- file.path(
  build_dir, "cache/model/SOURCE_PROGRAM_SCORES.csv"
)

meta_paths <- c(
  file.path(build_dir, "data/nested.csv"),
  file.path(dirname(build_dir), "data/nested.csv")
)

meta_path <- meta_paths[file.exists(meta_paths)][1]

stopifnot(file.exists(score_file), !is.na(meta_path))

score <- read_csv(score_file, show_col_types = FALSE)
nst   <- read_csv(meta_path, show_col_types = FALSE)

pcols <- grep(
  "^program_[0-9]+$",
  names(score),
  value = TRUE
)

pcols <- pcols[
  order(as.integer(str_extract(pcols, "[0-9]+")))
]

stopifnot(length(pcols) == 9)


# ── Paired case-control differences ───────────────────────────────────────────

pair_diff <- score %>%
  inner_join(
    nst %>%
      select(
        sentrix_id,
        origin,
        matched_index,
        group,
        time
      ),
    by = "sentrix_id"
  ) %>%
  filter(
    group %in% c("control", "prediagnostic")
  ) %>%
  mutate(
    pair = paste(origin, matched_index, sep = "_")
  ) %>%
  pivot_longer(
    all_of(pcols),
    names_to = "program",
    values_to = "score"
  ) %>%
  group_by(pair, program) %>%
  filter(
    n() == 2,
    n_distinct(group) == 2
  ) %>%
  summarise(
    time = as.numeric(
      time[group == "prediagnostic"][1]
    ),
    delta = score[group == "prediagnostic"][1] -
      score[group == "control"][1],
    .groups = "drop"
  ) %>%
  filter(
    is.finite(time),
    is.finite(delta)
  ) %>%
  mutate(
    program_lab = factor(
      paste0(
        "P",
        as.integer(str_extract(program, "[0-9]+"))
      ),
      levels = paste0("P", 1:9)
    )
  )


# ── Program colors ────────────────────────────────────────────────────────────

pal <- setNames(
  rep("#7396A5", 9),
  paste0("P", 1:9)
)

pal[c("P3", "P9")] <- "#292473"
pal[c("P5", "P8")] <- "#C2185B"


# ── Base LOESS plot ────────────────────────────────────────────────────────────

x_trans <- scales::pseudo_log_trans(sigma = 3)

p_time_3x3 <- ggplot(
  pair_diff,
  aes(
    time,
    delta,
    color = program_lab,
    fill = program_lab
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.2
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    span = 0.60,
    se = TRUE,
    level = 0.95,
    n = 250,
    alpha = 0.17,
    linewidth = 0.3
  ) +
  facet_wrap(
    ~ program_lab,
    ncol = 3,
    scales = "fixed",
    drop = FALSE
  ) +
  scale_color_manual(
    values = pal,
    guide = "none"
  ) +
  scale_fill_manual(
    values = pal,
    guide = "none"
  ) +
  scale_x_continuous(
    trans = x_trans,
    breaks = c(0, 1, 2, 5, 10, 15),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    # Extra space above for departure annotations
    expand = expansion(mult = c(0.10, 0.22))
  ) +
  labs(
    x = "Time-to-diagnosis",
    y = expression(
      Delta*S[p]~"(case - matched control)"
    ),
    tag = "f"
  ) +
  basic +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 7,
      face = "plain",
      hjust = 0.5
    ),
    panel.spacing = unit(0, "mm"),
    panel.grid = element_blank(),
    legend.position = "none"
  )


# ── Extract LOESS coordinates ─────────────────────────────────────────────────

pb <- ggplot_build(p_time_3x3)

# Layer 1 = geom_hline
# Layer 2 = geom_smooth

sm <- as_tibble(pb$data[[2]]) %>%
  mutate(
    PANEL = as.integer(PANEL)
  ) %>%
  left_join(
    as_tibble(pb$layout$layout) %>%
      transmute(
        PANEL = as.integer(PANEL),
        program_lab
      ),
    by = "PANEL"
  ) %>%
  mutate(
    year = x_trans$inverse(x)
  ) %>%
  filter(
    is.finite(year),
    is.finite(y)
  )


# ── Detect departure time ─────────────────────────────────────────────────────

# Departure criterion:
# Absolute deviation >= 0.10 from distant-time baseline

delta_cut <- 0.10

depart <- sm %>%
  group_by(program_lab) %>%
  group_modify(\(d, key) {
    
    d <- d %>%
      arrange(year)
    
    # Baseline from distant time-to-diagnosis
    far <- d$y[d$year >= 10]
    
    if (length(far) < 5) {
      far <- tail(
        d$y,
        max(5, ceiling(nrow(d) * 0.2))
      )
    }
    
    baseline <- median(far, na.rm = TRUE)
    
    deviation <- abs(d$y - baseline)
    
    # Find first transition from substantial
    # deviation toward distant-time baseline
    idx <- which(
      head(deviation, -1) >= delta_cut &
        tail(deviation, -1) < delta_cut
    )
    
    if (length(idx) == 0) {
      return(
        tibble(
          year = numeric(0),
          curve_y = numeric(0),
          baseline = numeric(0)
        )
      )
    }
    
    i <- idx[1]
    
    # Interpolation on transformed X coordinates
    x0 <- d$x[i]
    x1 <- d$x[i + 1]
    
    v0 <- deviation[i]
    v1 <- deviation[i + 1]
    
    w <- (delta_cut - v0) / (v1 - v0)
    
    cross_x <- x0 + w * (x1 - x0)
    
    # Exact height on the LOESS polyline
    cross_y <- d$y[i] +
      w * (d$y[i + 1] - d$y[i])
    
    tibble(
      year = x_trans$inverse(cross_x),
      curve_y = cross_y,
      baseline = baseline
    )
    
  }) %>%
  ungroup()


# ── Position triangles and labels ─────────────────────────────────────────────

# Shared Y-axis span because facets use fixed scales
yr <- pb$layout$panel_params[[1]]$y.range
dy <- diff(yr)

depart <- depart %>%
  mutate(
    label = sprintf("%.1f yr", year),
    
    # Downward triangle just above LOESS curve
    marker_y = curve_y + 0.13 * dy,
    
    # Text above the triangle
    label_y = curve_y + 0.26 * dy
  )


# ── Add departure annotations ─────────────────────────────────────────────────

p_time_3x3_depart <- p_time_3x3 +
  
  geom_point(
    data = depart,
    aes(
      x = year,
      y = marker_y
    ),
    inherit.aes = FALSE,
    shape = 25,
    size = 1.2,
    stroke = 0,
    fill = "#333333",
    color = "#333333",
    show.legend = FALSE
  ) +
  
  geom_text(
    data = depart,
    aes(
      x = year,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 1.7,
    family = "Arial",
    color = "#333333",
    hjust = 0.5,
    vjust = 0.5,
    show.legend = FALSE
  )


# ── Inspect departure times ───────────────────────────────────────────────────

print(
  depart %>%
    select(
      program_lab,
      year,
      baseline,
      curve_y
    )
)


# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(
  "result/tile-departure.svg",
  plot = p_time_3x3_depart,
  device = svglite::svglite,
  width = 2.8,
  height = 2.4,
  dpi = 800,
  bg = "transparent"
)

