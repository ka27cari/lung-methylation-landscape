# Plot outputs for: 02 states
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# Section: Microstate structure and labels (original 03_state-structure.r)
# Figure 4a/b | Sparse-state display correction | 
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)
library(patchwork)
library(grid)

# 1. Settings

affinity_k <- 10L
sparse_n <- 5L
output_dir <- "result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Muted scientific palette
macro_cols <- c(
  Source       = "#6F8792",
  Intermediate = "#74799A",
  Terminal     = "#8A687D"
)

# BEGIN SHARED DISPLAY MAP 
# Display policy shared byte-for-byte by 03_states, 04_ot and 04c.
# n_source = hard-state count in the complete source dataset, NOT prediagnostic n.
# Filter sparse non-terminal states FIRST, then number the displayed B/I/T states.
# Never use the display index to index a posterior column or transport matrix.
make_state_display_map <- function(roles, sparse_n = 5L) {
  required <- c("state", "n", "macrostate_role")
  if (!all(required %in% names(roles))) stop("Incomplete state-role table.")
  if (length(sparse_n) != 1L || !is.finite(sparse_n) || sparse_n < 0) {
    stop("sparse_n must be one non-negative source-count threshold.")
  }
  role_order <- c("Source", "Intermediate", "Terminal")
  prefix <- c(Source = "B", Intermediate = "I", Terminal = "T")
  raw_ids <- suppressWarnings(as.numeric(as.character(roles$state)))
  counts <- suppressWarnings(as.numeric(as.character(roles$n)))
  role <- as.character(roles$macrostate_role)
  stopifnot(
    all(is.finite(raw_ids)), all(raw_ids == as.integer(raw_ids)),
    !anyDuplicated(raw_ids), all(is.finite(counts)),
    all(counts >= 0), all(counts == as.integer(counts)),
    !anyNA(role), all(role %in% role_order)
  )
  all_states <- data.frame(
    state = as.integer(raw_ids), role = role, n_source = counts,
    stringsAsFactors = FALSE
  )
  all_states <- all_states[
    order(match(all_states$role, role_order), all_states$state), , drop = FALSE
  ]
  rownames(all_states) <- NULL
  all_states$included <- all_states$n_source >= sparse_n |
    all_states$role == "Terminal"
  
  # This is the ONLY map used for labels, ordering and plotted data joins.
  shown <- all_states[all_states$included, , drop = FALSE]
  if (nrow(shown) == 0L) stop("No states remain after the display filter.")
  shown$display_index <- integer(nrow(shown))
  for (r in role_order) {
    ii <- which(shown$role == r)
    shown$display_index[ii] <- seq_along(ii)
  }
  shown$display_label <- paste0(unname(prefix[shown$role]), shown$display_index)
  shown$display_order <- seq_len(nrow(shown))
  shown$plot_y <- nrow(shown) + 1L - shown$display_order
  stopifnot(!anyDuplicated(shown$display_label))
  
  # Legacy labels are retained ONLY in the crosswalk, never for plotting.
  old_index <- integer(nrow(all_states))
  for (r in role_order) {
    ii <- which(all_states$role == r)
    old_index[ii] <- seq_along(ii)
  }
  all_states$previous_label <- paste0(unname(prefix[all_states$role]), old_index)
  all_states$previous_03_label <- ifelse(
    all_states$role == "Source",
    sub("^B", "S", all_states$previous_label), all_states$previous_label
  )
  ii <- match(all_states$state, shown$state)
  all_states$display_label <- shown$display_label[ii]
  all_states$display_order <- shown$display_order[ii]
  all_states$plot_y <- shown$plot_y[ii]
  all_states$sparse_n <- sparse_n
  all_states$exclusion_reason <- ifelse(
    all_states$included, NA_character_, paste0("source n < ", sparse_n)
  )
  list(
    shown = tibble::as_tibble(shown),
    audit = tibble::as_tibble(all_states),
    keep_states = sort(shown$state)
  )
}

check_source_state_counts <- function(posterior, roles) {
  ids <- suppressWarnings(as.integer(as.character(posterior$hard_state)))
  role_ids <- as.integer(roles$state)
  stopifnot(!anyNA(ids), all(ids %in% role_ids))
  observed <- vapply(role_ids, function(s) sum(ids == s), numeric(1))
  if (!isTRUE(all.equal(observed, as.numeric(roles$n), check.attributes = FALSE))) {
    stop("roles$n does not match hard-state counts in the complete source posterior.")
  }
  invisible(TRUE)
}
# END SHARED DISPLAY MAP 

score_file <- file.path(build_dir, "cache/model/SOURCE_PROGRAM_SCORES.csv")
posterior_file <- file.path(build_dir, "cache/model/SOURCE_POSTERIOR.csv")
role_file <- file.path(build_dir, "analysis/tables/FINAL_STATE_ROLE_CLASSIFICATION.csv")
ot_file <- file.path(build_dir, "cache/model/SINKHORN_OT_PLAN.csv")

scores <- read_csv(score_file, show_col_types = FALSE)
posterior <- read_csv(posterior_file, show_col_types = FALSE)
roles <- read_csv(role_file, show_col_types = FALSE)
ot_raw <- read_csv(ot_file, show_col_types = FALSE)
cat("State roles: frozen LP-derived classification; plotted transport: Sinkhorn OT.\n")

# 2. Sanity checks

program_cols <- grep("^program_", names(scores), value = TRUE)
stopifnot(length(program_cols) == 9L)

stopifnot(
  !anyNA(scores$sentrix_id), !anyDuplicated(scores$sentrix_id),
  !anyNA(posterior$sentrix_id), !anyDuplicated(posterior$sentrix_id),
  setequal(scores$sentrix_id, posterior$sentrix_id)
)

roles <- roles %>%
  mutate(
    state = as.integer(state),
    macrostate_role = as.character(macrostate_role)
  )

posterior <- posterior %>%
  mutate(hard_state = as.integer(hard_state)) %>%
  slice(match(scores$sentrix_id, sentrix_id))

stopifnot(
  identical(as.character(scores$sentrix_id), as.character(posterior$sentrix_id))
)

stopifnot(
  setequal(sort(unique(posterior$hard_state)), sort(roles$state))
)

# 3. Filter first, then number the retained states

check_source_state_counts(posterior, roles)
K_model <- nrow(roles)                 # model dimension remains unchanged
stopifnot(K_model == 16L)
state_display <- make_state_display_map(roles, sparse_n)
state_display_map <- state_display$shown
keep_states <- state_display$keep_states

source_states <- state_display_map$state[state_display_map$role == "Source"]
intermediate_states <- state_display_map$state[state_display_map$role == "Intermediate"]
terminal_states <- state_display_map$state[state_display_map$role == "Terminal"]
state_order <- state_display_map$state
K <- length(state_order)              # displayed dimension, not model K

cat("Model states:", K_model, "; displayed states:", K, "\n")
cat("Excluded raw state IDs:", paste(setdiff(roles$state, keep_states), collapse = ", "), "\n")

# 4. Display labels: shared with 04_ot and 04c

display_map <- state_display_map %>%
  transmute(state, macrostate = role, display_label, order_id = display_order)
label_lookup <- setNames(display_map$display_label, as.character(display_map$state))
state_labels <- unname(label_lookup[as.character(state_order)])

# 5. Macrostate block positions

n_source <- length(source_states)
n_intermediate <- length(intermediate_states)
n_terminal <- length(terminal_states)

macro_block <- tibble(
  macrostate = c("Source", "Intermediate", "Terminal"),
  start = c(1L, n_source + 1L, n_source + n_intermediate + 1L),
  end = c(n_source, n_source + n_intermediate, K)
) %>%
  mutate(
    colour = unname(macro_cols[macrostate]),
    macro_label = case_when(
      macrostate == "Source"       ~ "Baseline",
      macrostate == "Intermediate" ~ "Interm.",
      macrostate == "Terminal"     ~ "Terminal"
    ),
    xmin = start - 0.5,
    xmax = end + 0.5,
    xmid = (start + end) / 2,
    ymin = K + 0.5 - end,
    ymax = K + 1.5 - start
  )

macro_boundaries <- c(
  n_source + 0.5,
  n_source + n_intermediate + 0.5
)

# 6. Fixed 9D rNMF sample coordinates

sample_df <- scores %>%
  left_join(
    posterior %>% select(sentrix_id, group, hard_state),
    by = "sentrix_id"
  )

stopifnot(!anyNA(sample_df$hard_state))

X <- as.matrix(sample_df[, program_cols])
storage.mode(X) <- "double"

# 7. Sample affinity

D <- as.matrix(dist(X, method = "euclidean"))

# Keep the original complete-source kernel. Hiding rows must not refit it.
n_sample_model <- nrow(D)
k_eff <- min(affinity_k, n_sample_model - 1L)

sigma <- apply(D, 1, function(z) {
  sort(z, decreasing = FALSE)[k_eff + 1L]
})

sigma[sigma <= 0] <- min(sigma[sigma > 0], na.rm = TRUE)

denom <- outer(sigma, sigma, FUN = "*")

affinity <- exp(
  -(D^2) / pmax(denom, .Machine$double.eps)
)

diag(affinity) <- 1

# 8. Sample ordering

sample_order <- unlist(
  lapply(state_order, function(s) {
    idx <- which(sample_df$hard_state == s)
    
    if (length(idx) <= 2L) {
      return(idx)
    }
    
    D_sub <- D[idx, idx, drop = FALSE]
    hc <- hclust(as.dist(D_sub), method = "average")
    idx[hc$order]
  }),
  use.names = FALSE
)

stopifnot(
  length(sample_order) == sum(sample_df$hard_state %in% keep_states),
  !anyDuplicated(sample_order),
  all(sample_df$hard_state[sample_order] %in% keep_states)
)
n_sample <- length(sample_order)      # only retained samples are drawn
cat("Affinity: calculated on", n_sample_model, "; displayed", n_sample, "samples.\n")

A_plot <- affinity[sample_order, sample_order, drop = FALSE]
# Explicit plotting indices avoid accidental use of original row names as positions.
dimnames(A_plot) <- list(as.character(seq_len(n_sample)), as.character(seq_len(n_sample)))

sample_affinity_long <- as.data.frame(
  as.table(A_plot),
  responseName = "affinity"
) %>%
  transmute(
    row = as.integer(Var1),
    col = as.integer(Var2),
    x = col,
    y = n_sample + 1L - row,
    affinity = affinity
  )

# 9. Sample block boundaries

state_counts <- vapply(
  state_order,
  function(s) sum(sample_df$hard_state == s),
  numeric(1)
)

state_cumulative <- cumsum(state_counts)

state_boundaries_x <- head(state_cumulative, -1) + 0.5
state_boundaries_y <- n_sample - head(state_cumulative, -1) + 0.5

macro_sample_counts <- c(
  Source = sum(sample_df$hard_state %in% source_states),
  Intermediate = sum(sample_df$hard_state %in% intermediate_states),
  Terminal = sum(sample_df$hard_state %in% terminal_states)
)

macro_sample_cumulative <- cumsum(macro_sample_counts)

macro_sample_boundaries_x <- head(macro_sample_cumulative, -1) + 0.5
macro_sample_boundaries_y <- n_sample - head(macro_sample_cumulative, -1) + 0.5

# 10. Sample affinity panel

p_affinity <- ggplot(
  sample_affinity_long,
  aes(x = x, y = y, fill = affinity)
) +
  geom_raster() +
  geom_vline(
    xintercept = state_boundaries_x,
    colour = scales::alpha("white", 0.25),
    linewidth = 0.12
  ) +
  geom_hline(
    yintercept = state_boundaries_y,
    colour = scales::alpha("white", 0.25),
    linewidth = 0.12
  ) +
  geom_vline(
    xintercept = macro_sample_boundaries_x,
    colour = scales::alpha("white", 0.25),
    linewidth = 0.12
  ) +
  geom_hline(
    yintercept = macro_sample_boundaries_y,
    colour = scales::alpha("white", 0.25),
    linewidth = 0.12
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1),
    trans = "sqrt",
    begin = 0.05,
    end = 0.95,
    guide = "none"
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(0.5, n_sample + 0.5),
    ylim = c(0.5, n_sample + 0.5)
  ) +
  labs(
    x = "Samples",
    y = "Samples",
    tag = "a"
  ) +
  basic +
  theme(
    axis.line = element_blank(),
    plot.tag.position  = c(0, 0.88),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    # plot.margin = margin(t = 7, r = 3, b = 4, l = 5)
  )

# 11. State-level affinity

state_affinity <- matrix(
  NA_real_,
  nrow = K,
  ncol = K,
  dimnames = list(as.character(state_order), as.character(state_order))
)

for (i in seq_along(state_order)) {
  s <- state_order[i]
  idx_s <- which(sample_df$hard_state == s)
  
  for (j in seq_along(state_order)) {
    t <- state_order[j]
    idx_t <- which(sample_df$hard_state == t)
    
    z <- affinity[idx_s, idx_t, drop = FALSE]
    
    if (s == t) {
      if (length(idx_s) >= 2L) {
        z <- z[row(z) != col(z)]
      } else {
        z <- NA_real_
      }
    }
    
    state_affinity[i, j] <- mean(z, na.rm = TRUE)
  }
}

state_affinity_long <- as.data.frame(
  as.table(state_affinity),
  responseName = "affinity"
) %>%
  transmute(
    origin_state = as.integer(as.character(Var1)),
    destination_state = as.integer(as.character(Var2)),
    affinity = affinity,
    x = match(destination_state, state_order),
    y = K + 1L - match(origin_state, state_order)
  )

# 12. State affinity panel

p_state <- ggplot(
  state_affinity_long,
  aes(x = x, y = y, fill = affinity)
) +
  geom_tile() +
  geom_vline(
    xintercept = macro_boundaries,
    colour = "white",
    linewidth = 0.1
  ) +
  geom_hline(
    yintercept = K + 1 - macro_boundaries,
    colour = "white",
    linewidth = 0.1
  ) +
  geom_rect(
    data = macro_block,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      colour = macrostate
    ),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = macro_cols) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    trans = "sqrt",
    begin = 0.08,
    end = 0.98,
    breaks = c(0, 0.1, 0.25, 0.5, 1),
    name = "Mean\naffinity",
    guide = guide_colourbar(
      barheight = unit(15, "mm"),
      barwidth = unit(2.1, "mm"),
      title.position = "top",
      title.hjust = 0.5,
      ticks = TRUE,
      frame.colour = NA
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(K),
    labels = state_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(K),
    labels = rev(state_labels),
    expand = c(0, 0)
  ) +
  annotate(
    "text",
    x = macro_block$xmid,
    y = K + 1.15,
    label = macro_block$macro_label,
    colour = macro_block$colour,
    family = "Arial",
    fontface = "plain",
    size = 1.85
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(0.5, K + 0.5),
    ylim = c(0.5, K + 0.5),
    clip = "off"
  ) +
  labs(
    x = "Microstate",
    y = "Microstate"
  ) +
  basic +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    
    # Enough headroom for the three macrostate labels
    # plot.margin = margin(t = 15, r = 0, b = 4, l = 0),
    legend.box.spacing = unit(1, "pt"),
    legend.spacing.x = unit(0, "pt"),
    
    legend.title = element_text(
      size = 5.5,
      lineheight = 0.9,
      margin = margin(b = 4)
    ),
    legend.text = element_text(
      size = 5,
      margin = margin(l = 1.5, unit = "pt")
    ),
    legend.spacing.y = unit(2, "mm")
  )

# 13. Arrow

make_arrow <- function(label) {
  ggplot() +
    annotate(
      "segment",
      x = 0.03,
      xend = 0.97,
      y = 0.50,
      yend = 0.50,
      linewidth = 0.28,
      lineend = "round",
      arrow = grid::arrow(
        angle = 18,
        length = unit(0.9, "mm"),
        type = "open"
      )
    ) +
    annotate(
      "text",
      x = 0.50,
      y = 0.56,
      label = label,
      family = "Arial",
      size = 2.2
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(l = 5, r = 0, b = 5, t = 5),
    )
}

arrow_kmeans <- make_arrow("k-means")

# 14. Final layout

final_plot <- (
  p_affinity +
    arrow_kmeans +
    p_state
) +
  plot_layout(
    widths = c(1.0, 0.15, 1.0),
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(
        fill = NA,
        colour = NA
      )
    )
  ) &
  theme(
    legend.position = "right",
    legend.box.spacing = unit(5.5, "pt")
  )

ggsave(file.path(output_dir, "states.svg"), device= svglite, final_plot, width = 3.9, height = 2,
       bg = "transparent")

# 15. Read the full frozen Sinkhorn plan; display its retained-state submatrix

names(ot_raw)[1] <- "origin_state"

ot_long0 <- ot_raw %>%
  mutate(
    origin_state = str_remove(
      as.character(origin_state),
      "^S"
    ),
    origin_state = as.integer(
      origin_state
    )
  ) %>%
  pivot_longer(
    cols = -origin_state,
    names_to = "destination_state",
    values_to = "mass"
  ) %>%
  mutate(
    destination_state = str_remove(
      destination_state,
      "^S"
    ),
    destination_state = as.integer(
      destination_state
    ),
    mass = as.numeric(
      mass
    )
  )

stopifnot(
  all(is.finite(ot_long0$mass)), all(ot_long0$mass >= 0),
  setequal(ot_long0$origin_state, roles$state),
  setequal(ot_long0$destination_state, roles$state),
  !anyDuplicated(paste(ot_long0$origin_state, ot_long0$destination_state))
)
# Display-only subsetting: retain original coupling values, no renormalization.
ot_long <- ot_long0 %>%
  dplyr::filter(origin_state %in% keep_states, destination_state %in% keep_states) %>%
  mutate(
    x = match(
      destination_state,
      state_order
    ),
    
    y = K + 1L -
      match(
        origin_state,
        state_order
      )
  )

# 16. Sinkhorn OT panel: retained-state display of the unchanged model
p_ot <- ggplot(
  ot_long,
  aes(x = x, y = y, fill = mass)
) + geom_tile() + geom_vline(
  xintercept = macro_boundaries,
  colour = "white",
  linewidth = 0.2
) + geom_hline(
  yintercept = K + 1 - macro_boundaries,
  colour = "white",
  linewidth = 0.2
) + geom_rect(
  data = macro_block,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    colour = macrostate
  ),
  inherit.aes = FALSE,
  fill = NA,
  linewidth = 0.55,
  show.legend = FALSE
) + scale_colour_manual(
  values = macro_cols
) + scale_fill_viridis_c(
  option = "inferno",
  name = "Sinkhorn\nOT mass",
  guide = guide_colourbar(
    barheight = unit(15, "mm"),
    barwidth = unit(2.1, "mm"),
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = NA
  )
) + scale_x_continuous(
  breaks = seq_len(K),
  labels = state_labels,
  expand = c(0, 0)
) + scale_y_continuous(
  breaks = seq_len(K),
  labels = rev(state_labels),
  expand = c(0, 0)
) + annotate(
  "text",
  x = macro_block$xmid,
  y = K + 1.15,
  label = macro_block$macro_label,
  colour = macro_block$colour,
  family = "Arial",
  fontface = "plain",
  size = 1.85
) + labs(
  x = "Destination microstate",
  y = "Origin microstate",
  tag = "d"
) + coord_fixed(
  xlim = c(0.5, K + 0.5),
  ylim = c(0.5, K + 0.5),
  clip = "off"
) + basic +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    # plot.margin = margin(t = 15, r = 4, b = 4, l = 4),
    legend.title = element_text(
      size = 5.5,
      lineheight = 0.9,
      margin = margin(b = 4)
    ),
    plot.tag.position  = c(0.0, 1.0),
    legend.text = element_text(
      size = 5,
      margin = margin(l = 1, unit = "pt")
    ),
    legend.spacing.y = unit(2, "mm"),
    legend.box.spacing = unit(1, "mm")
  )

ggsave(file.path(output_dir, "ot-flow.svg"), device = svglite, p_ot, width = 2.25, height = 2.05, dpi = 800)

# Section: Microstate-level profiles (original 07_state-profiles.r)
# Terminal-state histology composition
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)

# ── 1. Settings ───────────────────────────────────────────────────────────────

output_dir <- "result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

histology_levels <- c(
  "LUAD", "LUSC", "NSCLC", "SCLC", "Other", "Unknown"
)

histology_cols <- c(
  LUAD    = "#eb6f5d",
  LUSC    = "#71c9dd",
  NSCLC   = "#6376a0",
  SCLC    = "#33b39f",
  Other   = "#9da7c3",
  Unknown = "#a7dace"
)

# ── 2. Load data ──────────────────────────────────────────────────────────────

post <- read_csv(
  file.path(build_dir, "cache/model/SOURCE_POSTERIOR.csv"),
  show_col_types = FALSE
)

roles <- read_csv(
  file.path(build_dir, "analysis/tables/FINAL_STATE_ROLE_CLASSIFICATION.csv"),
  show_col_types = FALSE
)

stopifnot(
  all(c("group", "type", "hard_state") %in% names(post)),
  all(c("state", "macrostate_role") %in% names(roles))
)

# ── 3. Terminal states ────────────────────────────────────────────────────────

terminal_map <- roles %>%
  mutate(state = as.integer(state)) %>%
  filter(macrostate_role == "Terminal") %>%
  arrange(state) %>%
  mutate(
    terminal = paste0("T", row_number())
  ) %>%
  select(state, terminal)

# ── 4. Histology classification ───────────────────────────────────────────────

pre <- post %>%
  filter(tolower(as.character(group)) == "prediagnostic") %>%
  mutate(
    hard_state = as.integer(hard_state),
    hist_raw = toupper(trimws(as.character(type))),
    histology = case_when(
      is.na(hist_raw) |
        hist_raw %in% c(
          "", "NA", "NAN", "N/A", "NULL", "UNKNOWN"
        ) ~ "Unknown",
      hist_raw == "OTHER" ~ "Other",
      TRUE ~ hist_raw
    )
  ) %>%
  inner_join(
    terminal_map,
    by = c("hard_state" = "state")
  )

stopifnot(
  all(pre$histology %in% histology_levels)
)

# ── 5. Count histology by terminal ────────────────────────────────────────────

terminal_counts <- pre %>%
  count(terminal, histology, name = "n")

# Include zero-count categories
terminal_counts <- expand_grid(
  terminal = terminal_map$terminal,
  histology = histology_levels
) %>%
  left_join(
    terminal_counts,
    by = c("terminal", "histology")
  ) %>%
  mutate(
    n = replace_na(n, 0L),
    terminal = factor(
      terminal,
      levels = terminal_map$terminal
    ),
    histology = factor(
      histology,
      levels = histology_levels
    )
  )

# Check total counts
terminal_counts %>%
  group_by(terminal) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  print()

# ── 6. Stacked bar plot ───────────────────────────────────────────────────────

p_terminal <- ggplot(
  terminal_counts,
  aes(x = terminal, y = n, fill = histology)
) +
  geom_col(
    width = 0.7,
    position = position_stack(reverse = TRUE),
    colour = NA
  ) +
  scale_fill_manual(
    values = histology_cols,
    breaks = histology_levels,
    drop = FALSE,
    name = NULL,
    guide = guide_legend(
      nrow = 2,
      byrow = TRUE,
      keywidth = unit(2.5, "mm"),
      keyheight = unit(2.5, "mm")
    )
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = "Terminal state",
    y = "Cases (n)",
    tag = "c"
  ) +
  basic +
  theme(
    axis.text.x = element_text(
      size = 6,
      angle = 0
    ),
    panel.grid = element_blank(),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.title = element_blank(),
    legend.text = element_text(size = 6, margin = margin(l = 1, unit = "pt")),
    legend.key.spacing.x = unit(1, "mm"),
    legend.margin = margin(0, 0, 2, 0),
    legend.spacing.y = unit(2, "mm"),
    legend.box.spacing = unit(2, "mm")
  )

# Save
ggsave(
  file.path(output_dir, "term-profile.svg"),
  device = svglite,
  p_terminal,
  width = 1.6,
  height = 1.5,
  units = "in",
  dpi = 2000,
  bg = "transparent"
)

output_dir <- "result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_font <- "Arial"
plot_dpi <- 800
panel_tag <- "d"
sparse_n <- 5L
bar_height <- 0.7
smooth_bw <- 2.5
n_smooth_grid <- 800L
role_plot_order <- c("Source", "Intermediate", "Terminal")
role_labels <- c(Source = "Baseline", Intermediate = "Intermediate", Terminal = "Terminal")
smoking_levels <- c("None", "Quit", "Current")
smoking_cols <- c(None = "#c3c3c3", Quit = "#e37b7b", Current = "#e03131")

# Display order is the same as in the original multi-panel figure.
roles <- read_csv(file.path(build_dir, "analysis/tables/FINAL_STATE_ROLE_CLASSIFICATION.csv"), show_col_types = FALSE)
roles <- roles %>% transmute(state = as.integer(state), role = as.character(macrostate_role), n_source = as.integer(n))
stopifnot(!anyNA(roles$state), !anyDuplicated(roles$state),
          all(roles$role %in% role_plot_order), all(roles$n_source >= 0))
state_map <- roles %>%
  arrange(match(role, role_plot_order), state) %>%
  filter(n_source >= sparse_n | role == "Terminal") %>%
  group_by(role) %>%
  mutate(display_index = row_number()) %>%
  ungroup() %>%
  mutate(
    prefix = recode(role, Source = "B", Intermediate = "I", Terminal = "T"),
    state_label = paste0(prefix, display_index),
    plot_y = n() + 1L - row_number()
  )
keep_states <- state_map$state

audit_n <- state_map %>% count(role, name = "n_states")
print(audit_n)

# Load sources
post <- read_csv(file.path(build_dir, "cache/model/SOURCE_POSTERIOR.csv"), show_col_types = FALSE)
metadata <- read_csv(file.path(build_dir, "data/nested.csv"), show_col_types = FALSE)
J <- read.csv(file.path(build_dir, "cache/model/SINKHORN_NET_CURRENT.csv"),
              row.names = 1, check.names = FALSE) %>% as.matrix()
storage.mode(J) <- "double"
posterior_cols <- sprintf("P_state_%02d", 0:15)
stopifnot(
  all(c("sentrix_id", "group", "time", posterior_cols) %in% names(post)),
  all(c("sentrix_id", "smk_status") %in% names(metadata)),
  !anyDuplicated(post$sentrix_id), !anyDuplicated(metadata$sentrix_id),
  nrow(J) == ncol(J), identical(rownames(J), colnames(J)),
  all(is.finite(J)), all(J >= 0)
)
state_id_J <- as.integer(sub("^S", "", rownames(J)))
stopifnot(!anyNA(state_id_J), all(keep_states %in% state_id_J))

P <- as.matrix(post[, posterior_cols])
storage.mode(P) <- "double"
stopifnot(all(is.finite(P)), all(P >= 0), max(abs(rowSums(P) - 1)) < 1e-8)

# Prediagnostic samples; posterior weights retained for every displayed state
pre_base <- post %>%
  filter(tolower(as.character(group)) == "prediagnostic") %>%
  transmute(sentrix_id = as.character(sentrix_id), time = as.numeric(time),
            across(all_of(posterior_cols))) %>%
  left_join(metadata %>% transmute(sentrix_id = as.character(sentrix_id),
                                   smk_status = as.character(smk_status)), by = "sentrix_id") %>%
  mutate(
    smk_raw = tolower(trimws(smk_status)),
    smoking = case_when(
      smk_raw == "none" ~ "None", smk_raw == "quit" ~ "Quit",
      smk_raw == "current" ~ "Current", TRUE ~ "Unknown"
    )
  )
stopifnot(nrow(pre_base) > 0, all(is.finite(pre_base$time)), all(pre_base$time >= 0))
if (any(pre_base$smoking == "Unknown")) {
  smoking_levels <- c(smoking_levels, "Unknown")
  smoking_cols <- c(smoking_cols, Unknown = "#ececec")
}

pre_soft <- pre_base %>%
  pivot_longer(all_of(posterior_cols), names_to = "state_col", values_to = "weight") %>%
  mutate(state = as.integer(str_extract(state_col, "\\d+$"))) %>%
  filter(state %in% keep_states, is.finite(weight), weight >= 0) %>%
  left_join(state_map %>% select(state, state_label, role, plot_y), by = "state")
all_soft <- post %>%
  transmute(group = tolower(as.character(group)), across(all_of(posterior_cols))) %>%
  pivot_longer(all_of(posterior_cols), names_to = "state_col", values_to = "weight") %>%
  mutate(state = as.integer(str_extract(state_col, "\\d+$"))) %>%
  filter(state %in% keep_states, group %in% c("control", "prediagnostic"),
         is.finite(weight), weight >= 0) %>%
  left_join(state_map %>% select(state, state_label, role, plot_y), by = "state")

# (1) Posterior-weighted density; 0 years is LEFT, longer times RIGHT.
time_transform <- scales::pseudo_log_trans(sigma = 2, base = 10)
max_time <- max(pre_base$time)
stopifnot(max_time > 0)
time_plot_grid <- seq(0, time_transform$transform(max_time), length.out = n_smooth_grid)
time_grid_smooth <- time_transform$inverse(time_plot_grid)
smooth_time_mass <- pre_soft %>%
  group_by(state, state_label, role, plot_y) %>%
  group_modify(~ {
    z <- vapply(time_grid_smooth, function(g) {
      sum(.x$weight * (dnorm((g - .x$time) / smooth_bw) +
                         dnorm((g + .x$time) / smooth_bw)))
    }, numeric(1))
    zmax <- max(z, na.rm = TRUE)
    tibble(time_plot = time_plot_grid,
           intensity = if (is.finite(zmax) && zmax > 0) z / zmax else 0)
  }) %>%
  ungroup()

# (2) Posterior-weighted control proportion
control_prop <- all_soft %>%
  group_by(state, state_label, role, plot_y) %>%
  summarise(control_mass = sum(weight[group == "control"]),
            total_mass = sum(weight), .groups = "drop") %>%
  mutate(control_prop = control_mass / total_mass)
stopifnot(all(is.finite(control_prop$control_prop)),
          all(between(control_prop$control_prop, 0, 1)))

# (3) Original Sinkhorn directed current (In / Out)
state_flow <- tibble(state = state_id_J,
                     inflow = colSums(J), outflow = rowSums(J)) %>%
  inner_join(state_map %>% select(state, state_label, role, plot_y), by = "state")
flow_plot <- bind_rows(
  state_flow %>% transmute(plot_y, direction = "In", xmin = -inflow, xmax = 0),
  state_flow %>% transmute(plot_y, direction = "Out", xmin = 0, xmax = outflow)
)
flow_lim <- max(0.10, 1.05 * max(state_flow$inflow, state_flow$outflow))
flow_tick_label <- format(signif(flow_lim, 2), trim = TRUE)

# (4) Posterior-weighted smoking fractions
smoking_composition <- pre_soft %>%
  group_by(state, smoking) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  right_join(expand_grid(state = keep_states, smoking = smoking_levels),
             by = c("state", "smoking")) %>%
  left_join(state_map %>% select(state, plot_y), by = "state") %>%
  mutate(weight = replace_na(weight, 0),
         stack_order = match(smoking, rev(smoking_levels))) %>%
  arrange(state, stack_order) %>%
  group_by(state) %>%
  mutate(fraction = weight / sum(weight),
         xmax = cumsum(fraction), xmin = xmax - fraction) %>%
  ungroup() %>%
  mutate(smoking = factor(smoking, levels = smoking_levels))
stopifnot(all(is.finite(smoking_composition$fraction)))

# Aligned state-row geometry
y_limits <- c(0.5, nrow(state_map) + 0.5)
y_breaks <- state_map$plot_y
role_brackets <- state_map %>%
  group_by(role) %>%
  summarise(ymin = min(plot_y) - 0.23, ymax = max(plot_y) + 0.23,
            ymid = (ymin + ymax) / 2, .groups = "drop") %>%
  mutate(role_label = unname(role_labels[role]))
time_breaks <- c(0, 1, 2, 5, 10, 20)
time_breaks <- time_breaks[time_breaks <= max_time]

p_t <- ggplot(smooth_time_mass, aes(x = time_plot, y = plot_y, fill = intensity)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(colours = viridisLite::viridis(256), limits = c(0, 1), guide = "none") +
  scale_x_continuous(limits = c(0, time_transform$transform(max_time)),
                     breaks = time_transform$transform(time_breaks),
                     labels = time_breaks, expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks,
                     labels = state_map$state_label, expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Microstate", title = "Time-to-diagnosis\n(years)", tag = panel_tag) +
  basic +
  theme(text = element_text(family = plot_font),
        plot.margin = margin(l = 5, r = 1, b = 2, t = 4),
        axis.title.y = element_text(size = 8), axis.text = element_text(size = 6),
        panel.grid = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, margin = margin(t = 3, b = 1)))

p_role <- ggplot(role_brackets) +
  geom_segment(aes(x = 0.035, xend = 0.035, y = ymin, yend = ymax),
               linewidth = 0.18, lineend = "butt") +
  geom_text(aes(x = 0.09, y = ymid, label = role_label),
            family = plot_font, size = 2.15, hjust = 0) +
  scale_x_continuous(limits = c(0, 0.72), expand = c(0, 0)) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks,
                     expand = expansion(mult = c(0, 0))) +
  labs(title = " ") + theme_void(base_size = 9) +
  theme(text = element_text(family = plot_font),
        plot.title = element_text(size = 7, margin = margin(b = 5)),
        plot.margin = margin(l = 0, r = 5, b = 2, t = 4))

p_control <- ggplot(control_prop) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = control_prop,
      ymin = plot_y - bar_height / 2,
      ymax = plot_y + bar_height / 2
    ),
    fill = "#3d5e91",
    colour = NA
  ) +
  geom_vline(
    xintercept = 1/3,
    linetype = "22",
    linewidth = 0.2
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(
      breaks = 1/3,
      labels = "30%",
      name = NULL
    )
  ) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL, title = "Control\nProp.") + basic +
  theme(text = element_text(family = plot_font),
        axis.text.x.top = element_text(size = 5, margin = margin(b = 1)),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        plot.margin = margin(l = 0, r = 4, b = 2, t = 4),
        axis.text.x = element_text(size = 6), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 5)))

p_flow <- ggplot(flow_plot) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = plot_y - bar_height / 2, ymax = plot_y + bar_height / 2,
                fill = direction), colour = NA) +
  geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.22) +
  scale_fill_manual(values = c(In = "#168c80", Out = "#d1495b"), guide = "none") +
  scale_x_continuous(
    breaks = c(-0.11, 0, 0.11),
    labels = c("0.11", "0", "0.11"),
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(
      breaks = c(-0.05, 0.05),
      labels = c("In", "Out"),
      name = NULL
    )
  ) +
  coord_cartesian(xlim = c(-flow_lim, flow_lim), clip = "on") +
  scale_y_continuous(limits = y_limits, breaks = y_breaks,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL, title = "Net\ncurrent") + basic +
  theme(text = element_text(family = plot_font),
        plot.margin = margin(l = 5, r = 4, b = 2, t = 4),
        axis.text.x = element_text(size = 6),
        axis.text.x.top = element_text(size = 6, margin = margin(b = 1)),
        axis.ticks.x.top = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 1)))

p_smk <- ggplot(smoking_composition, aes(fill = smoking)) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = plot_y - bar_height / 2, ymax = plot_y + bar_height / 2), colour = NA) +
  scale_fill_manual(values = smoking_cols, breaks = smoking_levels, drop = FALSE,
                    name = "Smoking status", guide = guide_legend(name = NULL,
                                                           guide = guide_legend(nrow = 1, byrow = TRUE))) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL, title = "Smoking\nstatus") + basic +
  theme(text = element_text(family = plot_font),
        plot.margin = margin(l = 5, r = 0, b = 2, t = 4),
        axis.text.x = element_text(size = 6), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.text.y = element_blank(),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(2, "mm"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5, margin = margin(b = 5)))

# Four analytical panels + state role strip (no histology panel)
p4_b_soft <- (p_t | p_role | p_control | p_flow | p_smk) +
  plot_layout(widths = c(.45, .35, .17, .17, .17), guides = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(fill = NA, colour = NA))) &
  theme(
    # Panel tag
    plot.tag.location = "plot",
    plot.tag.position = c(0, 0.95),
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      margin = margin(l = -5)
    ),
    
    # Plot title
    plot.title = element_text(
      size = 7,
      hjust = 0.5,
      margin = margin(t = 3, b = 1)
    ),
    
    # Legend placement
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "right",
    legend.box.just = "right",
    
    # Legend spacing
    legend.box.spacing = unit(1.5, "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.x = unit(1, "mm"),
    
    # Legend title
    legend.title = element_text(
      family = plot_font,
      size = 6,
      margin = margin(r = 6)
    ),
    
    # Legend items
    legend.text = element_text(
      family = plot_font,
      size = 5,
      margin = margin(l = 1)
    ),
    
    legend.key.width = unit(2, "mm"),
    legend.key.height = unit(2, "mm")
  )

print(p4_b_soft)
ggsave(file.path(output_dir, "state-profiles.svg"), plot = p4_b_soft,
       device = svglite::svglite, width = 2.9, height = 1.9,
       units = "in", dpi = plot_dpi, bg = "transparent")

# Section: Terminal mass by time to diagnosis (original 08_tm-to-time.r)
library(ggbeeswarm)
library(tidyverse)
library(grid)

# ── 1. Settings ───────────────────────────────────────────────────────────────

plot_font <- "Arial"
output_dir <- "result"
plot_dpi <- 800

violin_fill <- "#D9E4E8"
point_col <- "#4F7080"
ctrl_col <- "#3D5E91"

point_size <- 0.55
point_alpha <- 0.25
median_point_size <- 1.5

# ── 2. Load current rebuild inputs ────────────────────────────────────────────

posterior_file <- file.path(build_dir, "cache", "model", "SOURCE_POSTERIOR.csv")
role_file <- file.path(build_dir, "analysis", "tables", "FINAL_STATE_ROLE_CLASSIFICATION.csv")

if (!file.exists(posterior_file) || !file.exists(role_file)) {
  stop("Check build_dir: SOURCE_POSTERIOR.csv or FINAL_STATE_ROLE_CLASSIFICATION.csv is missing.")
}

posterior <- read_csv(posterior_file, show_col_types = FALSE)

roles <- read_csv(role_file, show_col_types = FALSE) %>%
  mutate(
    state = as.integer(state),
    macrostate_role = as.character(macrostate_role)
  )

stopifnot(
  all(c("sentrix_id", "group", "time", "TM") %in% names(posterior)),
  all(c("state", "macrostate_role") %in% names(roles)),
  !anyDuplicated(posterior$sentrix_id),
  !anyDuplicated(roles$state),
  all(roles$macrostate_role %in% c("Source", "Intermediate", "Terminal"))
)

# ── 3. Identify terminal states and verify TM ─────────────────────────────────

terminal_states <- sort(roles$state[roles$macrostate_role == "Terminal"])

if (length(terminal_states) == 0L) {
  stop("No terminal states found in FINAL_STATE_ROLE_CLASSIFICATION.csv.")
}

terminal_cols <- paste0("P_state_", sprintf("%02d", terminal_states))

if (!all(terminal_cols %in% names(posterior))) {
  stop(
    "Terminal posterior column(s) missing: ",
    paste(setdiff(terminal_cols, names(posterior)), collapse = ", ")
  )
}

posterior <- posterior %>%
  mutate(TM_recalc = rowSums(across(all_of(terminal_cols))))

max_tm_error <- max(abs(posterior$TM - posterior$TM_recalc), na.rm = TRUE)

if (!is.finite(max_tm_error) || max_tm_error > 1e-8) {
  stop(
    "Stored TM does not match the sum of terminal posterior columns. Max error = ",
    signif(max_tm_error, 4)
  )
}

# ── 4. Prediagnostic source samples ───────────────────────────────────────────

D <- posterior %>%
  filter(
    tolower(as.character(group)) == "prediagnostic",
    is.finite(time),
    is.finite(TM)
  ) %>%
  transmute(
    sentrix_id = as.character(sentrix_id),
    time = as.numeric(time),
    time_plot = -as.numeric(time),
    TM = as.numeric(TM)
  )

stopifnot(
  nrow(D) > 0L,
  !anyDuplicated(D$sentrix_id),
  all(D$time >= 0),
  all(D$TM >= -1e-10 & D$TM <= 1 + 1e-10)
)

# ── 5. Continuous Spearman association ───────────────────────────────────────

rho_test <- suppressWarnings(
  cor.test(D$time_plot, D$TM, method = "spearman", exact = FALSE)
)

rho <- unname(rho_test$estimate)
p_value <- rho_test$p.value

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) "<0.001" else paste0("= ", formatC(p, format = "f", digits = 3))
}

stat_label <- sprintf(
  "Spearman %s = %.2f, p %s",
  "\u03c1", rho, format_p(p_value)
)

# ── 6. Time bins ──────────────────────────────────────────────────────────────

D <- D %>%
  mutate(
    time_bin = cut(
      time,
      # breaks = c(0, 1, 2, 5, 15, Inf),
      breaks = c(0, 1, 2, 5, Inf),
      right = FALSE,
      include.lowest = TRUE,
      # labels = c("<1", "1-2", "2-5", "5-15", "15+"),
      labels = c("<1", "1-2", "2-5", "5+")
    ),
    time_bin = factor(
      time_bin,
      # levels = c("15+", "5-15", "2-5", "1-2", "<1")
      levels = c("<1", "1-2", "2-5", "5+")
    )
  )

stopifnot(!anyNA(D$time_bin))

# ── 7. Summary within each interval ───────────────────────────────────────────

bin_summary <- D %>%
  group_by(time_bin) %>%
  summarise(
    n = n(),
    time_median = median(time, na.rm = TRUE),
    TM_mean = mean(TM, na.rm = TRUE),
    TM_median = median(TM, na.rm = TRUE),
    TM_q25 = quantile(TM, 0.25, na.rm = TRUE),
    TM_q75 = quantile(TM, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(bin_summary)

# ── 8. Control distribution ───────────────────────────────────────────────────

D_ctrl <- posterior %>%
  filter(
    tolower(as.character(group)) == "control",
    is.finite(TM)
  ) %>%
  transmute(
    TM = as.numeric(TM),
    side_group = factor("ctrl", levels = "ctrl")
  )

prob_tail_trans <- scales::new_transform(
  name = "arcsine_sqrt",
  transform = function(x) asin(sqrt(x)),
  inverse = function(x) sin(x)^2,
  domain = c(0, 1)
)

# ── 9. Plot ───────────────────────────────────────────────────────────────────

p_tm_time <- ggplot(
  D,
  aes(x = time_bin, y = TM)
) +
  
  geom_violin(
    width = 0.82,
    scale = "width",
    trim = TRUE,
    adjust = 0.75,
    fill = violin_fill,
    colour = NA
  ) +
  
  ggbeeswarm::geom_beeswarm(
    colour = point_col,
    size = point_size,
    alpha = point_alpha,
    shape = 16,
    cex = 0.72,
    priority = "density",
    corral = "gutter",
    corral.width = 0.34
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    colour = point_col,
    size = median_point_size,
    shape = 16
  ) +
  
  # Control distribution
  ggside::geom_ysideviolin(
    data = D_ctrl,
    aes(x = side_group, y = TM),
    inherit.aes = FALSE,
    orientation = "x",
    width = 0.8,
    color = NA,
    fill = ctrl_col
  ) +
  
  ggside::scale_ysidex_discrete(
    name = NULL,
    labels = c(ctrl = "ctrl")
  ) +
  
  scale_y_continuous(
    trans = prob_tail_trans,
    limits = c(0, 1),
    # breaks = c(
    #   0, 0.05, 0.10, 0.25,
    #   0.50, 0.75, 0.90, 0.95, 1
    # ),
    # labels = c(
    #   "0", ".05", ".10", ".25",
    #   ".50", ".75", ".90", ".95", "1"
    # ),
    breaks = c(
      0, 0.05, 0.25, 0.50, 0.75, 0.95, 1
    ),
    labels = c(
      "0", ".05", ".25", ".50", ".75", ".95", "1"
    ),
    expand = expansion(mult = c(0.00, 0.05))
  ) +
  
  labs(
    x = "Time-to-diagnosis (years)",
    y = "Terminal prob.",
    tag = "b"
  ) +
  
  basic +
  
  theme(
    text = element_text(family = plot_font),
    axis.text = element_text(size = 6),
    panel.grid = element_blank(),
    
    plot.tag.position = c(0, 0.95),
    
    ggside.panel.scale = 0.22,
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    
    legend.title = element_text(
      family = plot_font,
      size = 6
    ),
    legend.text = element_text(
      family = plot_font,
      size = 5
    ),
    legend.key.height = unit(8, "mm")
  )

p_tm_time

# Alias for Figure 4 assembly
p4_f <- p_tm_time

# ── 10. Save ──────────────────────────────────────────────────────────────────

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(output_dir, "source-time.svg"),
  device = svglite,
  p_tm_time,
  width = 1.6,
  height = 1.5,
  units = "in",
  dpi = plot_dpi,
  bg = "transparent"
)

# Section: Program scores across states (original 10_prog-state.r)
# All-state program pattern heatmap
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)
library(grid)
library(ggtext)

plot_font <- "Arial"
sparse_n <- 5L
panel_tag <- "f"

# 1. Program order / families --------------------------------------------------

display_order <- c("P1", "P2", "P6", "P3", "P9", "P7", "P5", "P8", "P4")

families <- list(
  ReprPC = c("P3", "P9"),
  Bivalent = c("P5", "P8")
)

family_colors <- c(
  ReprPC = "#292473",
  Bivalent = "#C2185B"
)

# 2. Visual settings -----------------------------------------------------------

# Stronger blue-white-red palette, closer to the supplied reference figure.
  low_col <- "#2c3173"
mid_col <- "#FAFAFA"
high_col <- "#de1c10"

moderate_cut <- 0.25
strong_cut <- 0.50
macro_gap <- 0.35

# 3. Load ----------------------------------------------------------------------

roles <- read_csv(
  file.path(build_dir, "analysis/tables/FINAL_STATE_ROLE_CLASSIFICATION.csv"),
  show_col_types = FALSE
) %>%
  mutate(state = as.integer(state), role = as.character(macrostate_role))

centroids <- read_csv(
  file.path(build_dir, "cache/model/STATE_CENTROIDS.csv"),
  show_col_types = FALSE
) %>%
  mutate(state = as.integer(state))

# 4. Display-state definition --------------------------------------------------

state_map <- roles %>%
  filter(n >= sparse_n | role == "Terminal") %>%
  arrange(factor(role, levels = c("Source", "Intermediate", "Terminal")), state) %>%
  group_by(role) %>%
  mutate(display_index = row_number()) %>%
  ungroup() %>%
  mutate(
    state_label = case_when(
      role == "Source" ~ paste0("B", display_index),
      role == "Intermediate" ~ paste0("I", display_index),
      role == "Terminal" ~ paste0("T", display_index)
    ),
    display_order = row_number(),
    macro_gap_index = case_when(
      role == "Source" ~ 2,
      role == "Intermediate" ~ 1,
      TRUE ~ 0
    ),
    plot_y = n() + 1L - display_order + macro_gap * macro_gap_index
  )

stopifnot(!anyDuplicated(state_map$state_label))

# 5. State x program table -----------------------------------------------------

plot_df <- centroids %>%
  inner_join(
    state_map %>% select(state, state_label, role, plot_y),
    by = "state"
  ) %>%
  select(state, state_label, role, plot_y, starts_with("program_")) %>%
  pivot_longer(
    cols = starts_with("program_"),
    names_to = "program_raw",
    values_to = "score"
  ) %>%
  mutate(
    program_index = as.integer(str_extract(program_raw, "[0-9]+")),
    program = paste0("P", program_index),
    plot_x = match(program, display_order),
    abs_score = abs(score),
    marker = case_when(
      abs_score >= strong_cut ~ "Strong",
      abs_score >= moderate_cut ~ "Moderate",
      TRUE ~ "None"
    )
  ) %>%
  filter(!is.na(plot_x))

stopifnot(nrow(plot_df) == nrow(state_map) * length(display_order))

# 6. Colour scale --------------------------------------------------------------

score_lim <- max(abs(plot_df$score), na.rm = TRUE)
if (!is.finite(score_lim) || score_lim <= 0) score_lim <- 1

# 7. Program-family annotations ------------------------------------------------

family_annot <- tibble(
  family = c("ReprPC", "Bivalent"),
  xmin = c(
    min(match(families$ReprPC, display_order)) - 0.38,
    min(match(families$Bivalent, display_order)) - 0.38
  ),
  xmax = c(
    max(match(families$ReprPC, display_order)) + 0.38,
    max(match(families$Bivalent, display_order)) + 0.38
  ),
  x = c(
    mean(match(families$ReprPC, display_order)),
    mean(match(families$Bivalent, display_order))
  ),
  y_line = max(state_map$plot_y) + 0.62,
  y_text = max(state_map$plot_y) + 1.12
)

program_tick_labels <- setNames(display_order, display_order)

program_tick_labels[families$ReprPC] <- paste0(
  "<span style='color:", unname(family_colors["ReprPC"]), "'>",
  families$ReprPC, "</span>"
)

program_tick_labels[families$Bivalent] <- paste0(
  "<span style='color:", unname(family_colors["Bivalent"]), "'>",
  families$Bivalent, "</span>"
)

# 8. Plot ----------------------------------------------------------------------

p_state_program <- ggplot(plot_df, aes(plot_x, plot_y)) +
  geom_tile(aes(fill = score), width = 1, height = 1, colour = NA) +
  geom_point(
    data = plot_df %>% filter(marker == "Moderate"),
    aes(shape = marker),
    size = 1.25, stroke = 0.5, colour = "black"
  ) +
  geom_point(
    data = plot_df %>% filter(marker == "Strong"),
    aes(shape = marker),
    size = 1.25, colour = "black"
  ) +
  geom_segment(
    data = family_annot,
    aes(x = xmin, xend = xmax, y = y_line, yend = y_line, colour = family),
    inherit.aes = FALSE,
    linewidth = 0.25,
    lineend = "round"
  ) +
  geom_text(
    data = family_annot,
    aes(x = x, y = y_text, label = family, colour = family),
    inherit.aes = FALSE,
    family = plot_font,
    size = 2.1
  ) +
  scale_fill_gradient2(
    low = low_col,
    mid = mid_col,
    high = high_col,
    midpoint = 0,
    limits = c(-score_lim, score_lim),
    oob = scales::squish,
    name = "Program score",
    guide = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(18, "mm"),
      barheight = unit(1.8, "mm"),
      ticks = FALSE,
      title.position = "top",
      title.hjust = 0.5,
      theme = theme(
        legend.title = element_text(
          size = 5.8,
          hjust = 0.5,
          margin = margin(b = 3)
        ),
        legend.text = element_text(size = 5.3)
      )
    )
  ) +
  scale_shape_manual(
    values = c(Moderate = 1, Strong = 16),
    guide = "none"
  ) +
  scale_colour_manual(values = family_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq_along(display_order),
    labels = unname(program_tick_labels[display_order]),
    limits = c(0.5, length(display_order) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = state_map$plot_y,
    labels = state_map$state_label,
    limits = c(0.5, max(state_map$plot_y) + 1.42),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, tag = panel_tag) +
  basic +
  theme(
    text = element_text(family = plot_font),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = ggtext::element_markdown(
      size = 6.2,
      colour = "black",
      margin = margin(t = 3)
    ),
    axis.text.y = element_text(
      size = 6.2,
      colour = "black",
      margin = margin(r = 2.8)
    ),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(
      family = plot_font,
      size = 5.8,
      hjust = 0.5,
      margin = margin(t = 1)
    ),
    legend.text = element_text(
      size = 5.3,
      margin = margin(t = 2)
    ),
    legend.box.spacing = unit(0, "mm"),
    plot.tag.position = c(0, 0.97),
    plot.margin = margin(l = 10, r = 5, b = 5, t = 5)
  )

print(p_state_program)

# 9. Save ----------------------------------------------------------------------

dir.create("result", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "result/state-program.svg",
  p_state_program,
  device = svglite,
  width = 1.6,
  height = 2.9,
  units = "in",
  bg = "transparent"
)


# Section: Program-axis density (original 13_prog-axis-density.r)
# Program-axis density | single ggplot, 3 true facets
# Figure: 6.66 cm x 4.30 cm | 
# Source: same score / posterior / role / route calculation as supplied R script.
# No patchwork and NO hand-drawn inter-facet divider lines.

# If needed, enable the project environment before running this file:

library(tidyverse)
library(grid)

if (!exists("build_dir", inherits = TRUE)) {
  stop("build_dir is not set. Source your project environment first.")
}

plot_font <- "Arial"
sparse_n <- 5L
zero_halfwidth <- 0.05

group_cols <- c(
  "T1/T2" = "#75686D",
  "T3/T4" = "#86647C",
  "T4" = "#527A7A",
  "T1/T2/T3" = "#75686D"
)

# 1. Load ----------------------------------------------------------------------

scores <- read_csv(
  file.path(build_dir, "cache", "model", "SOURCE_PROGRAM_SCORES.csv"),
  show_col_types = FALSE
)

post <- read_csv(
  file.path(build_dir, "cache", "model", "SOURCE_POSTERIOR.csv"),
  show_col_types = FALSE
) %>%
  mutate(group = tolower(as.character(group)))

roles <- read_csv(
  file.path(build_dir, "analysis", "tables", "FINAL_STATE_ROLE_CLASSIFICATION.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    state = as.integer(state),
    role = as.character(macrostate_role)
  )

route_scores <- read_csv(
  file.path(build_dir, "cache", "model", "SOURCE_ROUTE_SCORES.csv"),
  show_col_types = FALSE
)

# 2. Terminal map --------------------------------------------------------------

state_map <- roles %>%
  filter(n >= sparse_n | role == "Terminal") %>%
  arrange(factor(role, levels = c("Source", "Intermediate", "Terminal")), state) %>%
  group_by(role) %>%
  mutate(display_index = row_number()) %>%
  ungroup() %>%
  mutate(
    state_label = case_when(
      role == "Source" ~ paste0("B", display_index),
      role == "Intermediate" ~ paste0("I", display_index),
      role == "Terminal" ~ paste0("T", display_index)
    )
  )

terminal_map <- state_map %>%
  filter(role == "Terminal") %>%
  transmute(
    state,
    terminal = state_label,
    route_col = paste0("route_to_T", state)
  )

if (nrow(terminal_map) != 4L) {
  stop("Expected four terminal states; found ", nrow(terminal_map))
}

route_col_lookup <- setNames(
  terminal_map$route_col,
  terminal_map$terminal
)

stopifnot(
  all(c("T1", "T2", "T3", "T4") %in% names(route_col_lookup)),
  all(unname(route_col_lookup) %in% names(route_scores))
)

# 3. Program-axis scores + source route ----------------------------------------

sample_axes <- scores %>%
  transmute(
    sentrix_id = as.character(sentrix_id),
    P58 = (program_05 + program_08) / 2,
    P39 = (program_03 + program_09) / 2,
    P246 = (program_02 + program_04 + program_06) / 3
  ) %>%
  inner_join(
    post %>% transmute(sentrix_id = as.character(sentrix_id), group),
    by = "sentrix_id"
  ) %>%
  inner_join(
    route_scores %>%
      transmute(
        sentrix_id = as.character(sentrix_id),
        T1 = .data[[route_col_lookup[["T1"]]]],
        T2 = .data[[route_col_lookup[["T2"]]]],
        T3 = .data[[route_col_lookup[["T3"]]]],
        T4 = .data[[route_col_lookup[["T4"]]]]
      ),
    by = "sentrix_id"
  ) %>%
  filter(
    group == "prediagnostic",
    if_all(c(P58, P39, P246, T1, T2, T3, T4), is.finite)
  ) %>%
  mutate(
    T12 = T1 + T2,
    route_sum = T12 + T3 + T4
  ) %>%
  filter(route_sum > 0) %>%
  mutate(
    T12 = T12 / route_sum,
    T3 = T3 / route_sum,
    T4 = T4 / route_sum,
    route_branch = case_when(
      T12 >= T3 & T12 >= T4 ~ "T1/T2",
      T3 >= T12 & T3 >= T4 ~ "T3",
      TRUE ~ "T4"
    )
  )

print(sample_axes %>% count(route_branch))


# 4. Combine comparison data (unchanged branch definitions) --------------------

axis_levels <- c("mean(P5,P8)", "mean(P3,P9)", "mean(P2,P4,P6)")

d1 <- sample_axes %>%
  mutate(family = if_else(route_branch == "T1/T2", "T1/T2", "T3/T4"))
d2 <- sample_axes %>%
  mutate(family = if_else(route_branch == "T4", "T4", "T1/T2/T3"))
# P246 uses the same T1/T2 versus T3/T4 branch grouping as P58.
d3 <- d1

# 5. Pure data preparation per facet ------------------------------------------
# Density scaling, highlight cutoff, modes and arrow direction follow the
# make_panel() logic in the original file. The function returns data frames,
# NOT ggplot objects; one ggplot below renders all three facets together.

make_density_df <- function(x, group, y0, panel, n = 512L, adjust = 1) {
  x <- x[is.finite(x)]
  if (length(unique(x)) < 2L) {
    stop("At least two unique scores are required: ", panel, " / ", group)
  }
  den <- density(x, n = n, adjust = adjust, na.rm = TRUE)
  tibble(panel = panel, group = group, x = den$x,
         density = den$y, y0 = y0)
}

make_facet_data <- function(data, score_col, group_var, top_group, bottom_group,
                            panel, highlight_mode = c("focus", "both"),
                            focus_group = NULL,
                            arrow_mode = c("focus", "separation")) {
  highlight_mode <- match.arg(highlight_mode)
  arrow_mode <- match.arg(arrow_mode)
  
  d <- data %>%
    transmute(group = as.character(.data[[group_var]]),
              score = as.numeric(.data[[score_col]])) %>%
    filter(group %in% c(top_group, bottom_group), is.finite(score))
  
  if (!all(c(top_group, bottom_group) %in% unique(d$group))) {
    stop("Missing comparison group in facet: ", panel)
  }
  if ((highlight_mode == "focus" || arrow_mode == "focus") &&
      (is.null(focus_group) || !focus_group %in% c(top_group, bottom_group))) {
    stop("A valid focus_group is required in facet: ", panel)
  }
  
  panel_xr <- range(c(d$score, -zero_halfwidth, zero_halfwidth))
  panel_span <- diff(panel_xr)
  if (!is.finite(panel_span) || panel_span <= 0) panel_span <- 1
  panel_xlim <- c(panel_xr[1] - 0.025 * panel_span,
                  panel_xr[2] + 0.015 * panel_span)
  
  den_df <- bind_rows(
    make_density_df(d$score[d$group == top_group],
                    top_group, 0.95, panel),
    make_density_df(d$score[d$group == bottom_group],
                    bottom_group, 0.25, panel)
  )
  den_df <- den_df %>%
    mutate(y = y0 + density / max(density) * 0.36)
  
  summary_df <- d %>%
    group_by(group) %>%
    summarise(q1 = quantile(score, 0.25),
              med = median(score),
              q3 = quantile(score, 0.75), .groups = "drop") %>%
    mutate(
      panel = panel,
      y0 = if_else(group == top_group, 0.95, 0.25)
    )
  
  med_df <- summary_df %>% arrange(med, group)
  low_group <- med_df$group[1]
  high_group <- med_df$group[2]
  cutoff <- mean(med_df$med)
  
  if (highlight_mode == "both") {
    den_hi <- den_df %>%
      filter((group == low_group & x <= cutoff) |
               (group == high_group & x >= cutoff))
  } else {
    focus_med <- summary_df$med[match(focus_group, summary_df$group)]
    other_group <- setdiff(c(top_group, bottom_group), focus_group)
    other_med <- summary_df$med[match(other_group, summary_df$group)]
    den_hi <- den_df %>%
      filter(group == focus_group) %>%
      filter(if (focus_med >= other_med) x >= cutoff else x <= cutoff)
  }
  
  # Original arrow origin is the KDE mode, not a decorative right-edge arrow.
  modes <- den_df %>%
    group_by(group) %>%
    slice_max(order_by = density, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(group, mode_x = x, arrow_y = y0 + 0.4)
  arrow_len <- 0.16 * panel_span
  
  if (arrow_mode == "focus") {
    focus_mode <- modes$mode_x[match(focus_group, modes$group)]
    focus_y <- modes$arrow_y[match(focus_group, modes$group)]
    focus_med <- summary_df$med[match(focus_group, summary_df$group)]
    other_med <- summary_df$med[
      match(setdiff(c(top_group, bottom_group), focus_group), summary_df$group)
    ]
    direction_sign <- ifelse(focus_med >= other_med, 1, -1)
    arrow_df <- tibble(panel = panel, group = focus_group,
                       x = focus_mode,
                       xend = focus_mode + direction_sign * arrow_len,
                       y = focus_y)
  } else {
    # Both arrows: lower-median group points left, higher-median points right.
    arrow_df <- tibble(
      panel = panel, group = c(low_group, high_group),
      x = c(modes$mode_x[match(low_group, modes$group)],
            modes$mode_x[match(high_group, modes$group)]),
      xend = c(modes$mode_x[match(low_group, modes$group)] - arrow_len,
               modes$mode_x[match(high_group, modes$group)] + arrow_len),
      y = c(modes$arrow_y[match(low_group, modes$group)],
            modes$arrow_y[match(high_group, modes$group)])
    )
  }
  
  # Free-x faceting requires explicit bounds per facet. Include arrow tips.
  panel_xlim <- range(c(panel_xlim, arrow_df$x, arrow_df$xend))
  # Use the original plot's left-side group order, one label per ridge.
  label_df <- summary_df %>%
    mutate(
      x = med,
      y = y0 - 0.075
    )
  
  zero_df <- summary_df %>%
    transmute(panel, ymin = y0 - 0.035, ymax = y0 + 0.40)
  bounds_df <- tibble(panel = panel, x = panel_xlim,
                      y = rep(0.10, 2L))
  
  # Keep the same observed-score-based x window as the original code.
  # KDE tails outside it must not silently widen free-x facet axes.
  den_display <- den_df %>%
    filter(x >= panel_xlim[1], x <= panel_xlim[2])
  hi_display <- den_hi %>%
    filter(x >= panel_xlim[1], x <= panel_xlim[2])
  
  list(density = den_display, highlight = hi_display,
       summary = summary_df, arrows = arrow_df,
       labels = label_df, zero = zero_df, bounds = bounds_df)
}

prepared <- list(
  make_facet_data(d1, "P58", "family", "T3/T4", "T1/T2",
                  axis_levels[1], "focus", "T1/T2", "focus"),
  make_facet_data(d2, "P39", "family", "T4", "T1/T2/T3",
                  axis_levels[2], "focus", "T4", "focus"),
  make_facet_data(d3, "P246", "family", "T3/T4", "T1/T2",
                  axis_levels[3], "both", NULL, "separation")
)

combine_component <- function(name) {
  bind_rows(lapply(prepared, `[[`, name)) %>%
    mutate(panel = factor(panel, levels = axis_levels))
}

den_df <- combine_component("density")
den_hi <- combine_component("highlight")
summary_df <- combine_component("summary")
arrow_df <- combine_component("arrows")
label_df <- combine_component("labels")
zero_df <- combine_component("zero")
bounds_df <- combine_component("bounds")

# The printed checks also verify arrows inherit the computed KDE-mode origins.
print(arrow_df %>% select(panel, group, x, xend, y))

# 6. One ggplot, facet_wrap creates all panel divisions ------------------------
# No manual lines between panels. Only within-panel zero markers and axes.

p_program_density <- ggplot() +
  geom_blank(data = bounds_df, aes(x = x, y = y)) +
  geom_rect(
    data = zero_df,
    aes(xmin = -zero_halfwidth, xmax = zero_halfwidth,
        ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = "#F0EFF0", colour = NA
  ) +
  geom_segment(
    data = zero_df,
    aes(x = 0, xend = 0, y = ymin, yend = ymax),
    inherit.aes = FALSE, colour = "#949094",
    linetype = "dashed", linewidth = 0.20
  ) +
  geom_ribbon(
    data = den_df,
    aes(x = x, ymin = y0, ymax = y, fill = group, group = group),
    inherit.aes = FALSE, alpha = 0.30, colour = NA
  ) +
  geom_ribbon(
    data = den_hi,
    aes(x = x, ymin = y0, ymax = y, fill = group, group = group),
    inherit.aes = FALSE, alpha = 0.70, colour = NA
  ) +
  geom_line(
    data = den_df,
    aes(x = x, y = y, colour = group, group = group),
    inherit.aes = FALSE, linewidth = 0.30, lineend = "round"
  ) +
  geom_segment(
    data = summary_df,
    aes(x = q1, xend = q3, y = y0 - 0.02, yend = y0 - 0.02,
        colour = group),
    inherit.aes = FALSE, linewidth = 0.50, lineend = "round"
  ) +
  geom_point(
    data = summary_df,
    aes(x = med, y = y0 - 0.02, colour = group),
    inherit.aes = FALSE, size = 1.35
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE, colour = "#4D494B", linewidth = 0.35,
    arrow = grid::arrow(angle = 18, length = grid::unit(0.75, "mm"),
                        type = "closed", ends = "last")
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = group, colour = group),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1,
    family = plot_font,
    size = 2.0,
    show.legend = FALSE
  ) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_x_continuous(
    breaks = function(limits) {
      b <- pretty(limits, n = 3)
      b[b >= limits[1] & b <= limits[2]]
    },
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.07, 1.47),
    breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Program-axis score", y = NULL, tag = "f") +
  coord_cartesian(clip = "on") + 
  basic +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 6.0,
      margin = margin(t = 0, b = 0, unit = "mm")
    ),
    panel.spacing.x = grid::unit(1.8, "mm"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    plot.tag.position = c(0, 0.98),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    
  ); p_program_density

# plot.tag.position = c(0.005, 0.998),
# plot.margin = margin(t = 1.6, r = 1.5, b = 1.2, l = 3.5, unit = "mm"),

# 7. Output to exact physical dimensions --------------------------------------

ggsave(
  filename = "result/program-density.svg",
  plot = p_program_density, device = svglite::svglite,
  width = 6.8, height = 4.4, units = "cm", bg = "transparent"
)
# Optional raster verification at 600 dpi:
# ggsave("result/program-density-faceted-6.66x4.30cm.png",
#        p_program_density, width = 6.66, height = 4.30,
#        units = "cm", dpi = 600, bg = "white")

