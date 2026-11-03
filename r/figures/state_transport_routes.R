# Plot outputs for: 03 transport
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# Section: State transport flows (original 04_facet-flow.r)
# Figure 4d | Sparse-state display correction | 
source(file.path(repo, "r", "figures", "figure_environment.R"))
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(patchwork)
library(grid)

# Sinkhorn net-current route visualization
# - Same visual grammar as the original 04_ot.r
# - Node y-position = Sinkhorn net-current flow balance: (outflow-inflow)/(outflow+inflow)
# - Node x-position = one shared Sugiyama layout from the full retained net-current graph
# - Self flow is omitted
# - Two outputs: destination-time facets and eventual-terminal facets

# 1. Settings -----------------------------------------------------------------

sparse_n <- 5L
output_dir <- "result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
edge_strength <- 1

# - weak retained flows fade strongly by alpha
# - route-specific key states are bold + ◆ in terminal-route panels only

# Display-only pruning of weak connections.
# An edge is drawn only when it carries at least this fraction of the
# source state's total outgoing Sinkhorn net-current. The strongest
# outgoing edge from each source is always retained.
# This affects visualization only; the Sinkhorn current and shared layout are unchanged.
edge_share_min <- 0.05

edge_width_min <- 0.5
edge_width_max <- 2.0
edge_alpha_min <- 0.10
edge_alpha_max <- 0.82

facet_edge_scale <- 0.72
facet_node_size_min <- 0.60
facet_node_size_max <- 1.30
facet_label_size <- 2
facet_title_size <- 7
facet_label_gap <- 0.28

inactive_col <- "black"
flow_col <- "#D2D2D2"

horizon_levels <- c("Distant: ≥5y", "Mid: 3–5y", "Imminent: <3y")
role_levels <- c("Source", "Intermediate", "Terminal")

role_cols <- c(
  Source       = "#6F8792",
  Intermediate = "#74799A",
  Terminal     = "#8A687D"
)

role_labels <- c(
  Source = "Baseline",
  Intermediate = "Intermediate",
  Terminal = "Terminal"
)

# Route-specific visual anchors.
# These are DISPLAY labels, not raw state IDs.
route_key_states <- list(
  T1 = c("B3"),
  T2 = c("B3", "B5"),
  T3 = c("B2", "B4", "I3"),
  T4 = c("I2")
)
route_key_symbol <- "\u25C6"   # ◆; neutral route-marker, not a significance symbol

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

# 2. Load ---------------------------------------------------------------------

posterior <- read_csv(
  file.path(build_dir, "cache/model/SOURCE_POSTERIOR.csv"),
  show_col_types = FALSE
)

roles <- read_csv(
  file.path(build_dir, "analysis/tables/FINAL_STATE_ROLE_CLASSIFICATION.csv"),
  show_col_types = FALSE
)

J <- read.csv(
  file.path(build_dir, "cache/model/SINKHORN_NET_CURRENT_95.csv"),
  row.names = 1,
  check.names = FALSE
)

route_aff <- read_csv(
  file.path(build_dir, "analysis/tables/SINKHORN_TERMINAL_ROUTE_AFFINITY.csv"),
  show_col_types = FALSE
)

posterior$hard_state <- as.integer(posterior$hard_state)
roles$state <- as.integer(roles$state)
route_aff$state <- as.integer(route_aff$state)
route_aff$terminal_state <- as.integer(route_aff$terminal_state)

# 3. Filter globally, before assigning display labels --------------------------

check_source_state_counts(posterior, roles)
state_display <- make_state_display_map(roles, sparse_n)
state_display_map <- state_display$shown
keep_states <- state_display$keep_states
sparse_states <- state_display$audit$state[!state_display$audit$included]

source_states <- state_display_map$state[state_display_map$role == "Source"]
intermediate_states <- state_display_map$state[state_display_map$role == "Intermediate"]
terminal_states <- state_display_map$state[state_display_map$role == "Terminal"]
cat("Excluded raw sparse state IDs:", paste(sparse_states, collapse = ", "), "\n")

# 4. One contiguous display map for all time/terminal facets -------------------

display_map <- state_display_map %>% transmute(state, label = display_label)
label_lookup <- setNames(display_map$label, as.character(display_map$state))
terminal_label_lookup <- label_lookup[as.character(terminal_states)]

# 5. State-level residual time (facet assignment only) ------------------------

state_time <- posterior %>%
  filter(tolower(group) != "control", is.finite(time)) %>%
  group_by(state = hard_state) %>%
  summarise(
    pred_n = n(),
    median_time = median(time),
    .groups = "drop"
  ) %>%
  mutate(
    horizon = case_when(
      median_time < 3 ~ "Imminent: <3y",
      median_time < 5 ~ "Mid: 3–5y",
      TRUE ~ "Distant: ≥5y"
    ),
    horizon = factor(horizon, levels = horizon_levels)
  )

# 6. Nodes --------------------------------------------------------------------

nodes <- roles %>%
  filter(state %in% keep_states) %>%
  transmute(
    name = as.character(state),
    state,
    role = macrostate_role,
    n
  ) %>%
  left_join(state_time, by = "state") %>%
  left_join(display_map, by = "state")

stopifnot(all(is.finite(nodes$median_time)))

# 7. Net-current edge table ----------------------------------------------------

rownames(J) <- as.integer(str_remove(rownames(J), "^S"))
colnames(J) <- as.integer(str_remove(colnames(J), "^S"))

stopifnot(
  !anyDuplicated(rownames(J)), !anyDuplicated(colnames(J)),
  setequal(as.integer(rownames(J)), roles$state),
  setequal(as.integer(colnames(J)), roles$state),
  all(is.finite(as.matrix(J))), all(as.matrix(J) >= 0)
)
J <- J[
  as.character(keep_states),
  as.character(keep_states),
  drop = FALSE
]

edge_df <- as.data.frame(
  as.table(as.matrix(J)),
  stringsAsFactors = FALSE
) %>%
  transmute(
    source = as.integer(as.character(Var1)),
    target = as.integer(as.character(Var2)),
    mass = as.numeric(Freq)
  ) %>%
  filter(source != target, is.finite(mass), mass > 0)

role_lookup <- setNames(
  as.character(roles$macrostate_role),
  as.character(roles$state)
)

time_lookup <- setNames(nodes$median_time, as.character(nodes$state))
horizon_lookup <- setNames(as.character(nodes$horizon), as.character(nodes$state))

edge_df <- edge_df %>%
  mutate(
    source_role = unname(role_lookup[as.character(source)]),
    target_role = unname(role_lookup[as.character(target)]),
    target_time = unname(time_lookup[as.character(target)]),
    target_horizon = unname(horizon_lookup[as.character(target)]),
    target_horizon = factor(target_horizon, levels = horizon_levels)
  ) %>%
  group_by(source) %>%
  mutate(
    source_out_mass = sum(mass),
    source_share = if_else(source_out_mass > 0, mass / source_out_mass, 0),
    edge_display_keep = source_share >= edge_share_min | mass == max(mass)
  ) %>%
  ungroup()

stopifnot(nrow(edge_df) > 0)
stopifnot(all(is.finite(edge_df$target_time)))
stopifnot(all(!is.na(edge_df$target_horizon)))
stopifnot(any(edge_df$edge_display_keep))

cat("Retained Sinkhorn net-current edges:", nrow(edge_df), "\n")
cat(
  "Displayed edges after weak-edge pruning:",
  sum(edge_df$edge_display_keep), "/", nrow(edge_df),
  "(source-share cutoff =", edge_share_min, ")\n"
)

# 8. One shared layout for every facet ----------------------------------------

edges_graph <- edge_df %>%
  transmute(
    from = as.character(source),
    to = as.character(target),
    mass,
    source_role,
    target_role,
    target_time,
    target_horizon,
    edge_display_keep
  )

nodes_graph <- nodes %>%
  select(name, state, role, n, median_time, horizon, label)

g <- graph_from_data_frame(
  d = edges_graph,
  directed = TRUE,
  vertices = nodes_graph
)

sug <- layout_with_sugiyama(
  g,
  hgap = 1.3,
  vgap = 1,
  maxiter = 500
)

node_order <- tibble(
  name = V(g)$name,
  x = scales::rescale(sug$layout[, 1], to = c(-4.5, 4.5))
)

# Net-current flow balance: positive = exporter/source-like, negative = receiver/sink-like.
# This is used only as a visualization coordinate; terminal definitions remain LP-based.
J_mat <- as.matrix(J)
storage.mode(J_mat) <- "double"

outflow <- rowSums(J_mat)
inflow  <- colSums(J_mat)
flow_denom <- outflow + inflow
flow_balance <- ifelse(
  flow_denom > 0,
  (outflow - inflow) / flow_denom,
  0
)
names(flow_balance) <- rownames(J_mat)

nodes_plot <- nodes_graph %>%
  left_join(node_order, by = "name") %>%
  mutate(
    x = case_when(
      label == "T1" ~ x - 0.55,
      label == "T2" ~ x + 1.2,
      label == "T4" ~ x + 0.55,
      label == "B3" ~ x - 1.8,
      
      TRUE ~ x
    ),
    y = unname(flow_balance[name]),
    y = case_when(
      label == "B2" ~ y - 0.1,
      label == "I3" ~ y + 0.1,
      TRUE ~ y
    )
  )

g_tbl <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  left_join(
    nodes_plot %>% select(name, x, y),
    by = "name"
  ) %>%
  activate(edges) %>%
  mutate(
    target_horizon = factor(target_horizon, levels = horizon_levels),
    edge_width = scales::rescale(
      sqrt(mass),
      to = c(edge_width_min, edge_width_max)
    ),
    edge_alpha = scales::rescale(
      sqrt(mass),
      to = c(edge_alpha_min, edge_alpha_max)
    )
  )

# Label-side override for small-panel cleanup.
label_side_override <- c(
  "0" = "left",    # I1
  "8" = "right",
  "1" = "right",
  
  "3" = "right",
  
  "10" = "right",
  
  "6" = "right"
)

# Shared node styling helper ---------------------------------------------------

style_nodes <- function(g_sub, active_nodes, highlighted_labels = character()) {
  g_sub %>%
    activate(nodes) %>%
    mutate(
      route_highlight = label %in% highlighted_labels,
      label_plot = if_else(
        route_highlight,
        paste0(label, route_key_symbol),
        label
      ),
      label_face = if_else(route_highlight, "bold", "plain"),
      facet_node_size = scales::rescale(
        sqrt(n),
        to = c(facet_node_size_min, facet_node_size_max)
      ),
      facet_active = name %in% active_nodes,
      facet_display_role = case_when(
        route_highlight ~ "Highlight",
        facet_active ~ as.character(role),
        TRUE ~ "Inactive"
      ),
      facet_display_role = factor(
        facet_display_role,
        levels = c(role_levels, "Highlight", "Inactive")
      ),
      label_side = if_else(x <= 0, "left", "right"),
      label_side = if_else(
        name %in% names(label_side_override),
        unname(label_side_override[name]),
        label_side
      ),
      label_x = case_when(
        label_side == "left"  ~ x - facet_label_gap,
        label_side == "right" ~ x + facet_label_gap
      ),
      label_hjust = if_else(label_side == "left", 1, 0)
    )
}

base_graph_theme <- function() {
  list(
    scale_edge_width_identity(),
    scale_edge_alpha_identity(),
    scale_size_identity(guide = "none"),
    scale_colour_manual(
      values = c(
        role_cols,
        Highlight = "#de1c10",
        Inactive = inactive_col
      ),
      breaks = role_levels,
      labels = role_labels,
      name = "Macrostate",
      drop = FALSE,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "left",
        override.aes = list(
          shape = 16,
          size = 1.25,
          linetype = 0,
          linewidth = 0
        ),
        keywidth = unit(1.2, "mm"),
        keyheight = unit(2.5, "mm")
      )
    ),
    scale_x_continuous(
      expand = expansion(mult = c(0.11, 0.09))
    ),
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))),
    labs(x = NULL, y = NULL),
    basic,
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(
        size = facet_title_size,
        hjust = 0.5,
        margin = margin(b = 1)
      ),
      plot.margin = margin(
        t = 0, r = -3, b = 2, l = 0,
        unit = "pt"
      ),
      legend.title = element_text(
        size = 6,
        margin = margin(r = 3, unit = "pt")
      ),
      
      legend.text = element_text(
        size = 5.5,
        margin = margin(r = -1, unit = "pt")
      ),
      
      legend.key = element_blank(),
      legend.key.width = unit(3, "mm"),
      legend.key.height = unit(2, "mm"),
      legend.spacing.x = unit(0, "pt")
    )
  )
}

# 9. Destination-time panels ---------------------------------------------------

make_horizon_plot <- function(horizon, title = horizon) {
  facet_edge_df <- edge_df %>%
    filter(
      as.character(target_horizon) == horizon,
      edge_display_keep
    )
  
  facet_active_nodes <- unique(c(
    as.character(facet_edge_df$source),
    as.character(facet_edge_df$target)
  ))
  
  g_sub <- g_tbl %>%
    activate(edges) %>%
    filter(
      as.character(target_horizon) == horizon,
      edge_display_keep
    ) %>%
    mutate(facet_edge_width = edge_width * facet_edge_scale) %>%
    style_nodes(facet_active_nodes)
  
  p <- ggraph(g_sub, layout = "manual", x = x, y = y) +
    geom_edge_diagonal(
      aes(
        edge_width = facet_edge_width,
        edge_alpha = edge_alpha
      ),
      edge_colour = flow_col,
      strength = edge_strength,
      lineend = "round",
      start_cap = circle(0.1, "mm"),
      end_cap = circle(0.1, "mm"),
      show.legend = FALSE
    ) +
    geom_node_point(
      aes(size = facet_node_size, colour = facet_display_role),
      shape = 16,
      show.legend = TRUE,
      key_glyph = "point"
    ) +
    geom_node_text(
      aes(
        x = label_x,
        y = y,
        label = label_plot,
        hjust = label_hjust,
        colour = facet_display_role,
        fontface = label_face
      ),
      vjust = 0.5,
      size = facet_label_size,
      show.legend = FALSE
    ) +
    labs(title = title)
  
  for (z in base_graph_theme()) p <- p + z
  p
}

p_far <- make_horizon_plot(horizon_levels[1], "≥5y")
p_mid <- make_horizon_plot(horizon_levels[2], "3–5y")
p_near <- make_horizon_plot(horizon_levels[3], "<3y")

p_time <- (p_far + p_mid + p_near) +
  plot_layout(guides = "collect", widths = c(1, 1, 1)) &
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.margin = margin(l = -1, r = 0, b = 0, t = 7, unit = "pt"),
    legend.box.spacing = unit(4, "pt"),
    plot.margin = margin(l = 5, r = 5, b = 2, t = 2)
  )

# 10. Eventual-terminal route panels ------------------------------------------

# Each real net-current edge is apportioned to eventual terminal branches by
# the target state's net95 absorption affinity. No geometry-only edges are added.
branch_edges <- edge_df %>%
  filter(edge_display_keep) %>%
  select(source, target, mass, source_role, target_role) %>%
  inner_join(
    route_aff %>%
      select(state, terminal_state, net95_absorption_affinity),
    by = c("target" = "state")
  ) %>%
  mutate(
    branch_mass = mass * net95_absorption_affinity,
    terminal_label = unname(
      terminal_label_lookup[as.character(terminal_state)]
    )
  ) %>%
  filter(is.finite(branch_mass), branch_mass > 0)

if (nrow(branch_edges) > 0) {
  branch_edges <- branch_edges %>%
    mutate(
      edge_width = scales::rescale(
        sqrt(branch_mass),
        to = c(edge_width_min, edge_width_max)
      ),
      edge_alpha = scales::rescale(
        sqrt(branch_mass),
        to = c(edge_alpha_min, edge_alpha_max)
      )
    )
}

make_terminal_plot <- function(tlab) {
  be <- branch_edges %>% filter(terminal_label == tlab)
  
  # Build a facet graph with the same node set, then overwrite coordinates with
  # the one shared layout computed from the full retained net-current graph.
  eg <- be %>%
    transmute(
      from = as.character(source),
      to = as.character(target),
      mass = branch_mass,
      edge_width,
      edge_alpha
    )
  
  if (nrow(eg) == 0) {
    # Keep a node-only panel if a branch happens to have no retained edge.
    eg <- tibble(
      from = character(), to = character(), mass = numeric(),
      edge_width = numeric(), edge_alpha = numeric()
    )
  }
  
  gg <- graph_from_data_frame(
    d = eg,
    directed = TRUE,
    vertices = nodes_graph
  )
  
  active_nodes <- unique(c(eg$from, eg$to))
  
  highlighted_labels <- route_key_states[[tlab]]
  if (is.null(highlighted_labels)) highlighted_labels <- character()
  
  gt <- as_tbl_graph(gg) %>%
    activate(nodes) %>%
    left_join(nodes_plot %>% select(name, x, y), by = "name") %>%
    style_nodes(
      active_nodes,
      highlighted_labels = highlighted_labels
    )
  
  p <- ggraph(gt, layout = "manual", x = x, y = y)
  
  if (nrow(eg) > 0) {
    p <- p + geom_edge_diagonal(
      aes(
        edge_width = edge_width * facet_edge_scale,
        edge_alpha = edge_alpha
      ),
      edge_colour = flow_col,
      strength = edge_strength,
      lineend = "round",
      start_cap = circle(0.1, "mm"),
      end_cap = circle(0.1, "mm"),
      show.legend = FALSE
    )
  }
  
  p <- p +
    geom_node_point(
      aes(size = facet_node_size, colour = facet_display_role),
      shape = 16,
      show.legend = TRUE,
      key_glyph = "point"
    ) +
    geom_node_text(
      aes(
        x = label_x,
        y = y,
        label = label_plot,
        hjust = label_hjust,
        colour = facet_display_role,
        fontface = label_face
      ),
      vjust = 0.5,
      size = facet_label_size,
      show.legend = FALSE
    ) +
    labs(title = paste0(tlab, " route"))
  
  for (z in base_graph_theme()) p <- p + z
  p
}

terminal_display_levels <- unname(terminal_label_lookup[as.character(terminal_states)])
p_terminal_list <- lapply(terminal_display_levels, function(s) {
  make_terminal_plot(s) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
})

p_terminal <- wrap_plots(
  p_terminal_list,
  nrow = 1,
  guides = "collect"
) + plot_annotation(theme = theme(plot.background = element_rect(fill = NA, colour = NA))) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    
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
      size = 6,
      margin = margin(l = 1)
    ),
    
    legend.key.width = unit(2, "mm"),
    legend.key.height = unit(2, "mm"),
    
    plot.margin = margin(l = 3, r = 3, b = 2, t = 2)
  ) &
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

p_terminal <- wrap_elements(p_terminal) +
  labs(tag = "e") +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14
    ),
    plot.tag.position = c(0.03, 0.99)
  )

# 11. Save --------------------------------------------------------------------

print(p_time)
print(p_terminal)

ggsave(
  file.path(output_dir, "facet-flow.svg"),
  p_terminal,
  device = svglite,
  width = 4.2,
  height = 1.8,
  units = "in",
  bg = "transparent"
)


# Section: Terminal trajectory composition (original 09_ternary.r)
# Figure 4f | Source route ternary
# T1/T2 are collapsed only at the pre-terminal route level.
# Terminal identities remain separate elsewhere.


library(tidyverse)
library(ggtern)
library(grid)

plot_font <- "Arial"
sparse_n  <- 5L
panel_tag <- "e"

role_cols <- c(
  Source       = "#6F8792",
  Intermediate = "#74799A"
)

tern_major_breaks <- c(0.2, 0.4, 0.6, 0.8)
tern_minor_breaks <- seq(0.1, 0.9, by = 0.1)

# 1. Load route information

roles <- read_csv(
  file.path(build_dir, "analysis", "tables", "FINAL_STATE_ROLE_CLASSIFICATION.csv"),
  show_col_types = FALSE
) %>%
  mutate(state = as.integer(state))

route_aff <- read_csv(
  file.path(build_dir, "analysis", "tables", "SINKHORN_TERMINAL_ROUTE_AFFINITY.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    state = as.integer(state),
    terminal_state = as.integer(terminal_state)
  )

J <- read.csv(
  file.path(build_dir, "cache", "model", "SINKHORN_NET_CURRENT_95.csv"),
  row.names = 1,
  check.names = FALSE
) %>%
  as.matrix()

storage.mode(J) <- "double"

keep_states <- sort(
  roles$state[
    roles$n >= sparse_n |
      roles$macrostate_role == "Terminal"
  ]
)

source_states <- sort(
  roles$state[
    roles$macrostate_role == "Source" &
      roles$state %in% keep_states
  ]
)

intermediate_states <- sort(
  roles$state[
    roles$macrostate_role == "Intermediate" &
      roles$state %in% keep_states
  ]
)

terminal_states <- sort(
  roles$state[
    roles$macrostate_role == "Terminal" &
      roles$state %in% keep_states
  ]
)

if (length(terminal_states) != 4L) {
  stop("Expected four terminal states; found ", length(terminal_states))
}

state_map <- bind_rows(
  tibble(
    state = source_states,
    role = "Source",
    label = paste0("B", seq_along(source_states))
  ),
  tibble(
    state = intermediate_states,
    role = "Intermediate",
    label = paste0("I", seq_along(intermediate_states))
  ),
  tibble(
    state = terminal_states,
    role = "Terminal",
    label = paste0("T", seq_along(terminal_states))
  )
)

label_lookup <- setNames(state_map$label, state_map$state)
role_lookup  <- setNames(state_map$role, state_map$state)

# 2. Net-current edges

rownames(J) <- as.integer(str_remove(rownames(J), "^S"))
colnames(J) <- as.integer(str_remove(colnames(J), "^S"))

J <- J[
  as.character(keep_states),
  as.character(keep_states),
  drop = FALSE
]

diag(J) <- 0

edge_df <- as.data.frame(as.table(J), stringsAsFactors = FALSE) %>%
  transmute(
    source = as.integer(as.character(Var1)),
    target = as.integer(as.character(Var2)),
    mass = as.numeric(Freq)
  ) %>%
  filter(
    source != target,
    is.finite(mass),
    mass > 0
  )

# 3. Propagate edges toward terminal routes

branch <- edge_df %>%
  inner_join(
    route_aff %>%
      select(
        state,
        terminal_state,
        net95_absorption_affinity
      ),
    by = c("target" = "state")
  ) %>%
  mutate(
    branch_mass = mass * net95_absorption_affinity
  ) %>%
  filter(
    is.finite(branch_mass),
    branch_mass > 0
  )

bi_states <- c(
  source_states,
  intermediate_states
)

route_pref_long <- branch %>%
  filter(source %in% bi_states) %>%
  group_by(source, terminal_state) %>%
  summarise(
    branch_mass = sum(branch_mass),
    .groups = "drop"
  ) %>%
  complete(
    source = bi_states,
    terminal_state = terminal_states,
    fill = list(branch_mass = 0)
  ) %>%
  group_by(source) %>%
  mutate(
    total_branch_mass = sum(branch_mass),
    preference = if_else(
      total_branch_mass > 0,
      branch_mass / total_branch_mass,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    terminal_label = unname(
      label_lookup[as.character(terminal_state)]
    ),
    state_label = unname(
      label_lookup[as.character(source)]
    ),
    role = unname(
      role_lookup[as.character(source)]
    )
  )

# 4. Collapse T1 + T2 only at route level

route_pref <- route_pref_long %>%
  select(
    source,
    state_label,
    role,
    terminal_label,
    preference,
    total_branch_mass
  ) %>%
  pivot_wider(
    names_from = terminal_label,
    values_from = preference,
    values_fill = 0
  ) %>%
  distinct() %>%
  mutate(
    T12 = T1 + T2,
    route_sum = T12 + T3 + T4,
    
    role = factor(
      role,
      levels = c("Source", "Intermediate")
    ),
    
    state_label = factor(
      state_label,
      levels = c(
        paste0("B", seq_along(source_states)),
        paste0("I", seq_along(intermediate_states))
      )
    )
  ) %>%
  filter(
    total_branch_mass > 0,
    route_sum > 0
  ) %>%
  mutate(
    T12 = T12 / route_sum,
    T3  = T3  / route_sum,
    T4  = T4  / route_sum
  )

stopifnot(
  max(abs(route_pref$T12 + route_pref$T3 + route_pref$T4 - 1)) < 1e-8
)

# 5. Manual text coordinates
#
# IMPORTANT:
# These coordinates control TEXT ONLY.
#
# T12 = top vertex
# T3  = left vertex
# T4  = right vertex
#
# B2 and I3 share the same actual point at the T3 vertex, so their labels
# are deliberately separated inside the triangle.
#
label_position_override <- tribble(
  ~state_label, ~label_T12, ~label_T3, ~label_T4,
  
  "B5", 0.880, 0.055, 0.065,
  "B3", 0.710, 0.150, 0.140,
  
  "B7", 0.510, 0.300, 0.190,
  "B1", 0.392, 0.364, 0.244,
  "B6", 0.400, 0.080, 0.520,
  
  "I1", 0.196, 0.412, 0.392,
  "I2", 0.070, 0.100, 0.830,
  
  "B4", 0.092, 0.674, 0.234,
  
  "B2", 0.050, 0.890, 0.060,
  "I3", 0.120, 0.820, 0.060
)
# Make sure every plotted state has exactly one label coordinate.
plot_labels <- as.character(route_pref$state_label)

missing_labels <- setdiff(
  plot_labels,
  label_position_override$state_label
)

extra_labels <- setdiff(
  label_position_override$state_label,
  plot_labels
)

if (length(missing_labels)) {
  stop(
    "Missing manual label position(s): ",
    paste(missing_labels, collapse = ", ")
  )
}

if (length(extra_labels)) {
  warning(
    "Unused manual label position(s): ",
    paste(extra_labels, collapse = ", ")
  )
}

stopifnot(
  all(
    abs(
      label_position_override$label_T12 +
        label_position_override$label_T3 +
        label_position_override$label_T4 -
        1
    ) < 1e-3
  )
)

label_plot_df <- route_pref %>%
  mutate(
    state_label = as.character(state_label)
  ) %>%
  inner_join(
    label_position_override,
    by = "state_label"
  )

# 6. Plot

p_route_ternary <-
  ggtern(
    data = route_pref,
    aes(
      x = T3,
      y = T12,
      z = T4
    )
  ) +
  
  ggtern::geom_mask() +
  
  geom_point(
    aes(
      size = total_branch_mass,
      colour = role
    ),
    alpha = 0.9,
    shape = 16
  ) +
  
  geom_text(
    data = label_plot_df,
    aes(
      x = label_T3,
      y = label_T12,
      z = label_T4,
      label = state_label
    ),
    inherit.aes = FALSE,
    position = "identity",
    family = plot_font,
    colour = "black",
    size = 1.8,
    hjust = 0.5,
    vjust = 0.5,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  
  scale_T_continuous(
    name = "T1/T2",
    breaks = tern_major_breaks,
    minor_breaks = tern_minor_breaks,
    limits = c(0, 1)
  ) +
  
  scale_L_continuous(
    name = "T3",
    breaks = tern_major_breaks,
    minor_breaks = tern_minor_breaks,
    limits = c(0, 1)
  ) +
  
  scale_R_continuous(
    name = "T4",
    breaks = tern_major_breaks,
    minor_breaks = tern_minor_breaks,
    limits = c(0, 1)
  ) +
  
  scale_colour_manual(
    values = role_cols,
    breaks = c("Source", "Intermediate"),
    labels = c("Baseline", "Intermediate"),
    name = "State"
  ) +
  
  scale_size_continuous(
    range = c(1, 3),
    name = "Route flow",
    breaks = c(0.025, 0.050, 0.075, 0.100)
  ) +
  
  labs(
    tag = panel_tag
  ) +
  
  ggtern::theme_bw(
    base_size = 9,
    base_family = plot_font
  ) +
  
  ggtern::theme_showgrid() +
  ggtern::theme_nomask() +
  
  guides(
    size = guide_legend(
      order = 2,
      title.position = "top",
      direction = "vertical"
    ),
    colour = "none"
  ) +
  
  theme(
    tern.axis.clockwise = TRUE,
    tern.axis.arrow.show = FALSE,
    
    tern.panel.grid.major = element_line(
      colour = "#D0D0D0",
      linewidth = 0.1
    ),
    
    tern.panel.grid.minor = element_line(
      colour = "#ECECEC",
      linewidth = 0.05
    ),
    
    tern.axis.text = element_text(
      family = plot_font,
      size = 5.5
    ),
    
    tern.axis.title = element_text(
      family = plot_font,
      size = 9.5
    ),
    
    legend.position = "right",
    legend.justification = c(1, 0.95),
    legend.direction = "vertical",
    legend.box = "vertical",
    
    legend.margin = margin(
      b = 4,
      unit = "pt"
    ),
    
    legend.box.spacing = unit(
      0,
      "mm"
    ),
    
    legend.box.margin = margin(
      l = -25,
      r = 0,
      t = 0,
      b = 0,
      unit = "pt"
    ),
    
    legend.title = element_text(
      size = 6,
      margin = margin(
        b = 1,
        unit = "pt"
      )
    ),
    
    legend.text = element_text(
      size = 5.5
    ),
    
    legend.key.height = unit(
      2,
      "mm"
    ),
    
    legend.key.spacing.y = unit(
      0,
      "mm"
    ),
    
    legend.key.width = unit(
      2,
      "mm"
    ),
    
    legend.spacing.y = unit(
      0,
      "mm"
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.1
    ),
    
    tern.axis.line = element_line(
      colour = "black",
      linewidth = 0.1
    ),
    
    tern.axis.line.T = element_line(
      colour = "black",
      linewidth = 0.1
    ),
    
    tern.axis.line.L = element_line(
      colour = "black",
      linewidth = 0.1
    ),
    
    tern.axis.line.R = element_line(
      colour = "black",
      linewidth = 0.1
    ),
    
    plot.margin = margin(
      l = 0,
      r = 0,
      b = 0,
      t = 0
    ),
    
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      margin = margin(l = -5)
    ),
    
    plot.tag.position = c(
      0.1,
      0.95
    ),
    
    legend.background = element_blank(),
    
    tern.axis.title.T = element_text(
      size = 6.5
    ),
    
    tern.axis.title.L = element_text(
      size = 6.5
    ),
    
    tern.axis.title.R = element_text(
      size = 6.5
    )
  )

# 7. Output

p_source_route_ternary <- p_route_ternary

ggsave(
  "result/ternary.svg",
  p_route_ternary,
  device = svglite,
  width = 2,
  height = 1.9,
  units = "in",
  bg = "transparent"
)


# Section: Route-integrated trends (original 19_route-integrated-trend.r)
# Route-aware program-score trend: one regression line
# Selection: strongest B -> I -> T bottleneck; ties broken by edge product.
# If no Intermediate -> Terminal edge exists, use strongest direct B -> T.
# Selected relays are representative connected subpaths, not necessarily dominant terminal inflow.
# Single regression is fitted to unique microstate means, not repeated samples.
library(tidyverse)

output_dir <- "result"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_font <- "Arial"
min_flow <- 1e-6

roles <- read_csv(file.path(build_dir, "analysis/tables/FINAL_STATE_ROLE_CLASSIFICATION.csv"), show_col_types = FALSE)
post <- read_csv(file.path(build_dir, "cache/model/SOURCE_POSTERIOR.csv"), show_col_types = FALSE)
scores <- read_csv(file.path(build_dir, "cache/model/SOURCE_PROGRAM_SCORES.csv"), show_col_types = FALSE)
J <- read.csv(file.path(build_dir, "cache/model/SINKHORN_NET_CURRENT.csv"),
              row.names = 1, check.names = FALSE) |> as.matrix()
storage.mode(J) <- "double"
stopifnot(identical(rownames(J), colnames(J)), all(is.finite(J)), all(J >= 0))

# Use the same B/I/T naming rule as the Figure
state_map <- roles %>%
  transmute(state = as.integer(state), role = as.character(macrostate_role), n = as.integer(n)) %>%
  filter(n >= 5 | role == "Terminal") %>%
  arrange(factor(role, levels = c("Source", "Intermediate", "Terminal")), state) %>%
  group_by(role) %>%
  mutate(state_label = paste0(recode(role, Source = "B", Intermediate = "I", Terminal = "T"), row_number())) %>%
  ungroup()

# Directed Sinkhorn net current: row = outgoing, column = incoming
edges <- as_tibble(as.data.frame(as.table(J))) %>%
  transmute(from = as.integer(sub("^S", "", as.character(Var1))),
            to = as.integer(sub("^S", "", as.character(Var2))), flow = as.numeric(Freq)) %>%
  left_join(state_map %>% transmute(from = state, from_role = role), by = "from") %>%
  left_join(state_map %>% transmute(to = state, to_role = role), by = "to") %>%
  filter(!is.na(from_role), !is.na(to_role), flow > min_flow)

b_i <- edges %>%
  filter(from_role == "Source", to_role == "Intermediate") %>%
  transmute(b = from, i = to, bi = flow)
i_t <- edges %>%
  filter(from_role == "Intermediate", to_role == "Terminal") %>%
  transmute(i = from, t = to, it = flow)

relay_choice <- inner_join(b_i, i_t, by = "i") %>%
  mutate(bottleneck = pmin(bi, it), edge_product = bi * it) %>%
  group_by(t) %>%
  arrange(desc(bottleneck), desc(edge_product), b, i, .by_group = TRUE) %>%
  slice_head(n = 1) %>% ungroup()

direct_choice <- edges %>%
  filter(from_role == "Source", to_role == "Terminal") %>%
  group_by(to) %>%
  arrange(desc(flow), from, .by_group = TRUE) %>%
  slice_head(n = 1) %>% ungroup() %>%
  transmute(t = to, direct_b = from, direct_flow = flow)

route_choice <- state_map %>%
  filter(role == "Terminal") %>%
  transmute(t = state, route = state_label) %>%
  left_join(relay_choice, by = "t") %>%
  left_join(direct_choice, by = "t") %>%
  mutate(b = coalesce(b, direct_b),
         path_type = if_else(is.na(i), "Direct", "Via I"))
stopifnot(!anyNA(route_choice$b))

route_nodes <- bind_rows(
  route_choice %>% transmute(route, stage = "B", state = b, path_type),
  route_choice %>% filter(!is.na(i)) %>% transmute(route, stage = "I", state = i, path_type),
  route_choice %>% transmute(route, stage = "T", state = t, path_type)
) %>%
  left_join(state_map %>% select(state, state_label), by = "state")

# Summary of program scores among prediagnostic samples hard-assigned to each state
program_cols <- sprintf("program_%02d", 1:9)
stopifnot(all(program_cols %in% names(scores)))
pre_scores <- post %>%
  filter(tolower(as.character(group)) == "prediagnostic") %>%
  transmute(sentrix_id = as.character(sentrix_id), state = as.integer(hard_state)) %>%
  inner_join(scores %>% mutate(sentrix_id = as.character(sentrix_id)), by = "sentrix_id") %>%
  mutate(
    mean9 = rowMeans(across(all_of(program_cols))),
    P58 = (program_05 + program_08) / 2,
    P39 = (program_03 + program_09) / 2,
    P246 = (program_02 + program_04 + program_06) / 3
  )
state_scores <- pre_scores %>%
  group_by(state) %>%
  summarise(n_case = n(), across(c(mean9, P58, P39, P246), mean), .groups = "drop")

plot_df <- route_nodes %>%
  left_join(state_scores, by = "state") %>%
  mutate(
    stage = factor(stage, levels = c("B", "I", "T")),
    route = factor(route, levels = paste0("T", 1:4)),
    route_program = case_when(
      route %in% c("T1", "T2") ~ P58,
      route == "T3" ~ P246,
      route == "T4" ~ P39
    )
  ) %>%
  arrange(route, stage)
stopifnot(!anyNA(plot_df$n_case), all(is.finite(plot_df$mean9)), all(is.finite(plot_df$route_program)))

# How much current enters each terminal directly from B vs from I?
source_nodes <- state_map %>% filter(role == "Source") %>% pull(state)
inter_nodes <- state_map %>% filter(role == "Intermediate") %>% pull(state)
terminal_flow <- route_choice %>%
  transmute(
    route,
    direct_current = vapply(t, function(k) sum(J[paste0("S", source_nodes), paste0("S", k)]), numeric(1)),
    relay_current = vapply(t, function(k) sum(J[paste0("S", inter_nodes), paste0("S", k)]), numeric(1))
  ) %>%
  mutate(direct_pct = 100 * direct_current / (direct_current + relay_current),
         relay_pct = 100 - direct_pct)

print(route_choice %>% select(route, b, i, t, bi, it, direct_flow, path_type))
print(plot_df %>% select(route, stage, state_label, n_case, mean9, route_program), n = Inf)
print(terminal_flow)


# Plot: four route colors, one common regression line (microstate is the unit)
stopifnot(!anyDuplicated(plot_df$state))
plot_df <- plot_df %>% mutate(stage_ord = as.integer(stage) - 1L)
fit <- lm(mean9 ~ stage_ord, data = plot_df)
print(summary(fit))

fit <- lm(mean9 ~ stage_ord, data = plot_df)
trend_beta <- coef(fit)[["stage_ord"]]
trend_p <- coef(summary(fit))["stage_ord", "Pr(>|t|)"]


route_cols <- c(T1 = "#6550B3", T2 = "#009E9E", T3 = "#C99B00", T4 = "#737A84")
p <- ggplot(plot_df, aes(x = stage_ord, y = mean9)) +
  geom_hline(yintercept = 0, colour = "#C8C8C8", linetype = "22", linewidth = 0.15) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#353535", linewidth = 0.35) +
  geom_point(
    aes(colour = route),
    position = position_jitter(width = 0.07, height = 0, seed = 42),
    size = 1.05
  ) +
  scale_colour_manual(values = route_cols, breaks = paste0("T", 1:4), name = NULL,
                      guide = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
  scale_x_continuous(
    breaks = 0:2,
    labels = c("Baseline", "Intermediate", "Terminal"),
    limits = c(-0.2, 2.2),
    expand = c(0, 0)
  ) +
  annotate(
    "text",
    x = -0.16, y = Inf,
    label = sprintf("β = %+.3f\nP = %.4f", trend_beta, trend_p),
    hjust = -0.1, vjust = 1.4,
    size = 1.8, family = "Arial"
  ) +
  labs(x = "Macrostate", y = "Mean program score", tag = "g") +
  basic +
  theme(
    text = element_text(family = plot_font),
    panel.grid = element_blank(),
    legend.position = "top", legend.direction = "horizontal",
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    
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
      size = 6,
      margin = margin(l = 1)
    ),
    
    legend.key.width = unit(2, "mm"),
    legend.key.height = unit(2, "mm"),
    plot.tag.position  = c(-0.02, 0.95),
    plot.margin = margin(l = 13, r = 5, b = 5, t = 5),
    
  )

ggsave(file.path(output_dir, "route-program-single-trend.svg"), p,
       device = svglite::svglite, width = 5.2, height = 3.25,
       units = "cm", bg = "transparent")
print(p)
