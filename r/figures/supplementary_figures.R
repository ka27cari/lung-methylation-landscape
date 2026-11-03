# Plot outputs for: 07 supplementary
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO")
source(file.path(repo, "r", "figures", "figure_environment.R"))

# Section: Reproducibility and landscape sensitivity (original 20_s2-reproducible.r)
source(file.path(repo, "r", "figures", "figure_environment.R"))
# S2a-c: add before the three existing S2 panels (which become d-f).
# Discovery cohort only. K = 12:18, tau = 0.10:1.00, LP95.

read_build <- function(...) read.csv(file.path(build_dir, ...), check.names = FALSE)

# a. Why K = 16 and tau = 0.20?
k <- read_build("analysis", "tables", "K_SELECTION_STABILITY_AUDIT.csv")
t <- read_build("analysis", "tables", "POSTERIOR_TAU_CALIBRATION.csv")

t_start <- ggplot(t, aes(tau, loo_binary_nll)) +
  geom_ribbon(
    aes(ymin = loo_binary_nll - loo_binary_nll_se,
        ymax = loo_binary_nll + loo_binary_nll_se),
    alpha = .2,
    fill = "#56839A"
  ) +
  scale_x_continuous(expand = 0) +
  geom_vline(xintercept = .20, linetype = "22", color = "#666666", linewidth = 0.2) +
  geom_line(color = "#284B63", linewidth = 0.3) +
  labs(x = expression(tau), y = "LOO binary NLL", tag = "a") + basic +
  theme(
    plot.tag.position = c(-0.032, 0.95),
    plot.margin = margin(l = 13, r = 7, b = 5, t = 5)
  )

# Save
ggsave("result/tau-start.svg",
       device = svglite, t_start, width = 1.6, height = 1.3, dpi = 800)

# Source K-means starts: the supplied TSV belongs in 05_ncc/object/.
scan <- read_build("analysis", "tables", "K_START_ADMISSIBILITY_SCAN.csv")
scan <- scan[scan$passes %in% c(TRUE, "True", "TRUE"), ]
centroids <- read.delim("object/K_centroids.tsv", check.names = FALSE)
S <- read_build("cache", "model", "SOURCE_PROGRAM_SCORES.csv")
P <- read_build("cache", "model", "SOURCE_POSTERIOR.csv")

programs <- sprintf("program_%02d", 1:9)
X <- as.matrix(S[, programs])
p <- P[match(S$sentrix_id, P$sentrix_id), ]
pred <- p$group == "prediagnostic"
years <- -p$time[pred]
reference <- p$hard_state
tau_grid <- round(seq(.10, 1, by = .05), 2)
tau_display <- seq(.10, 1, by = .05)

choose2 <- function(x) x * (x - 1) / 2
baseline_same <- sum(choose2(table(reference)))
baseline_diff <- choose2(nrow(S)) - baseline_same

source_results <- list()
pair_results <- list()
offset <- 1

for (i in seq_len(nrow(scan))) {
  r <- scan[i, ]
  K <- r$K
  centres <- as.matrix(centroids[offset:(offset + K - 1), programs])
  offset <- offset + K
  
  distance <- vapply(seq_len(K), function(j)
    rowSums(sweep(X, 2, centres[j, ], "-")^2), numeric(nrow(X)))
  group <- max.col(-distance) - 1L
  
  overlap <- table(factor(reference, levels = 0:15),
                   factor(group, levels = 0:(K - 1)))
  same <- sum(choose2(overlap))
  sizes <- table(factor(group, levels = 0:(K - 1)))
  pair_results[[i]] <- data.frame(
    K = K, same = same / baseline_same,
    different = (sum(choose2(sizes)) - same) / baseline_diff)
  
  sc <- sapply(seq_len(K), function(j) median(distance[group == j - 1, j]))
  floor <- max(1e-4, .15 * median(sc[is.finite(sc) & sc > 0]))
  sc[!is.finite(sc) | sc < floor] <- floor
  terminals <- as.integer(strsplit(r$terminal_states, ";", fixed = TRUE)[[1]]) + 1L
  
  source_results[[i]] <- do.call(rbind, lapply(tau_grid, function(tau) {
    z <- -sweep(distance, 2, 2 * tau * sc, "/")
    z <- z - apply(z, 1, max)
    q <- exp(z) / rowSums(exp(z))
    tm <- rowSums(q[, terminals, drop = FALSE])
    
    test <- cor.test(
      tm[pred], years,
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      K = K,
      tau = tau,
      start = r$start,
      rho = unname(test$estimate),
      p = test$p.value
    )
  }))
}

source_detail <- bind_rows(source_results)

source_grid <- source_detail %>%
  group_by(K, tau) %>%
  summarise(
    rho = median(rho),
    p_worst = max(p),
    sign_consistent = all(rho > 0) || all(rho < 0),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_worst, method = "BH"),
    significant = p_adj < 0.05 & sign_consistent
  )

source_grid <- subset(source_grid, tau %in% tau_display)
source_grid$K <- factor(source_grid$K, levels = 18:12)

src_grid <- ggplot(source_grid, aes(factor(tau), K, fill = rho)) +
  geom_tile(colour = "white",
            linewidth = 0.25) +
  geom_tile(
    data = subset(source_grid, as.character(K) == "16" &
                    abs(tau - 0.20) < 1e-8),
    fill = NA,
    colour = "#D62728",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  labs(
    x = expression(tau),
    y = "K",
    fill = "Correlation", tag = "b") +
  basic + theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 0.5
    ),
    legend.box.spacing = unit(1, "mm"),
    legend.title.position = "right",
    legend.title = element_text(
      size = 6,
      angle = 270,
      margin = margin(l = 3)
    )
  ) + scale_fill_gradient(
    low = "#EAF2FB",
    high = "#2455A4",
    guide = guide_colorbar(
      title.position = "right",
      title.theme = element_text(
        angle = 270,
        vjust = 0.5
      )
    )
  )

# Save
ggsave("result/src-grid.svg",
       device = svglite, src_grid, width = 2.9, height = 1.3, dpi = 800)

# c. Original K16 same-group versus different-group pair co-clustering.
pair_detail <- bind_rows(pair_results)

pair_grid <- pair_detail %>%
  pivot_longer(
    c(same, different),
    names_to = "baseline",
    values_to = "probability"
  ) %>%
  group_by(K, baseline) %>%
  summarise(
    mean = mean(probability),
    lo = quantile(probability, .10),
    hi = quantile(probability, .90),
    .groups = "drop"
  )

pair_grid$baseline <- factor(
  pair_grid$baseline,
  levels = c("same", "different"),
  labels = c("Same", "Other")
)

pair_cols <- c(
  "Same" = "#284B63",
  "Other" = "black"
)

co_cluster <- ggplot(pair_grid, aes(K, mean, group = baseline,
                      colour = baseline, linetype = baseline)) +
  # geom_vline(xintercept = 16, linetype = "22", color = "#666666", linewidth = 0.2) +
  geom_point(size = 0.2) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    alpha = 0.2,
    colour = NA,
    fill = "#56839A",
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.3) +
   scale_colour_manual(values = pair_cols) +
  scale_fill_manual(values = pair_cols) +
  scale_x_continuous(
    breaks = 12:18, 
    expand = 0) +
  labs(
    x = "Alternative K",
    y = "Co-clustering",
    colour = NULL,
    linetype = NULL,
    tag = "c"
  ) +
  basic +
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.box.spacing = unit(1, "mm"),
    legend.key.width = unit(5, "mm")
  )

# Save
ggsave("result/co-clustering.svg",
       device = svglite, co_cluster, width = 2.5, height = 1.3, dpi = 800)

# Section: External/secondary validation (original 21_s3-validation.R)
source(file.path(repo, "r", "figures", "figure_environment.R"))
read_build <- function(...) read.csv(file.path(build_dir, ...), check.names = FALSE)
prog <- sprintf("program_%02d", 1:9)

S <- read_build("cache", "model", "SOURCE_PROGRAM_SCORES.csv")
L <- read_build("cache", "longitudinal", "LONGITUDINAL_POSTERIOR.csv")
P <- read_build("cache", "model", "SOURCE_POSTERIOR.csv")
P <- P[match(S$sentrix_id, P$sentrix_id), ]

X <- as.matrix(S[, prog])
Y <- as.matrix(L[, prog])
pred_source <- P$group == "prediagnostic"
years_source <- -P$time[pred_source]

ids <- unique(L$matched_index)
pred <- which(L$group == "prediagnostic")[match(ids, L$matched_index[L$group == "prediagnostic"])]
case <- which(L$group == "case")[match(ids, L$matched_index[L$group == "case"])]
far <- L$time[pred] > 5

# One *terminal-coherent consensus* per K. Constructed from admissible
# source-only starts by the same coassignment-basin rule as the primary K=16.
# See prepare_s3b_consensus.py; K=16 is checked against the original model.
consensus <- read.delim("object/K_consensus.tsv", check.names = FALSE)
audit <- read.delim("object/K_consensus_audit.tsv", check.names = FALSE)
if (!all(audit$status == "terminal_coherent_basin")) {
  stop("A K has no terminal-coherent consensus; inspect object/K_consensus_audit.tsv")
}
stopifnot(all(vapply(split(consensus, consensus$K), nrow, integer(1)) ==
                as.integer(names(split(consensus, consensus$K)))))

tau_grid <- round(seq(.10, 1, by = .05), 2)
tau_display <- round(seq(.10, .50, by = .05), 2)
rows <- list()

for (K in sort(unique(consensus$K))) {
  ck <- consensus[consensus$K == K, ]
  ck <- ck[order(ck$state), ]
  stopifnot(identical(as.integer(ck$state), 0:(K - 1L)))
  centres <- as.matrix(ck[, prog])
  distance <- function(x) sapply(seq_len(K), function(j)
    rowSums(sweep(x, 2, centres[j, ], "-")^2))
  ds <- distance(X)
  dl <- distance(Y)
  sc <- ck$scale
  terminals <- ck$state[ck$terminal == 1] + 1L
  
  tm <- function(d, tau) {
    z <- -sweep(d, 2, 2 * tau * sc, "/")
    z <- z - apply(z, 1, max)
    q <- exp(z) / rowSums(exp(z))
    rowSums(q[, terminals, drop = FALSE])
  }
  
  rows[[as.character(K)]] <- bind_rows(lapply(tau_grid, function(tau) {
    a <- tm(ds, tau)
    b <- tm(dl, tau)
    delta <- b[case] - b[pred]
    data.frame(K = K, tau = tau,
               linear_r = cor(a[pred_source], years_source),
               far = median(delta[far]),
               p_far = wilcox.test(delta[far], alternative = "greater", exact = FALSE)$p.value)
  }))
}

grid <- bind_rows(rows) %>%
  mutate(stars = case_when(far > 0 & p_far < .01 ~ "**",
                           far > 0 & p_far < .05 ~ "*",
                           far > 0 & p_far < .10 ~ ".", TRUE ~ ""))
plot_grid <- subset(grid, tau %in% tau_display)
plot_grid$K <- factor(plot_grid$K, levels = 18:12)

pb <- ggplot(plot_grid, aes(factor(tau), K, fill = linear_r)) +
  geom_tile() +
  geom_text(aes(label = stars), size = 2) +
  geom_tile(data = subset(plot_grid, K == "16" & tau == .20),
            aes(factor(tau), K), inherit.aes = FALSE,
            fill = NA, colour = "#D62728", linewidth = .8) +
  labs(x = "Tau", y = "K", tag= "b") +
  basic  + theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 270,
      hjust = 0,
      vjust = 0.5
    ),
    legend.box.spacing = unit(1, "mm"),
    legend.justification = "top"
  ) + scale_fill_gradient(
    low = "#EAF2FB",
    high = "#2455A4",
    name = "Pdx\nPearson r"
  )

# Save
ggsave(
  file.path("result/long-valid.svg"),
  device = svglite,
  pb,
  width = 2.4, height = 1.6, units = "in",bg = "transparent"
)

# c. CpG direction agreement by program and external effect magnitude.
M <- read_build("data", "m-ext.csv")
ph <- read_build("data", "ph-ext.csv")
orientation <- read_build("cache", "discovery", "SOURCE_CPG_ORIENTATION.csv")
top <- read_build("anchor", "nmf", "a0p0200", "stable", "STABLE_PROGRAM_TOP1000_CPG.csv")
selected <- read_build("cache", "model", "SELECTED_PROGRAMS.csv")

Z <- as.matrix(M[, ph$sample]); storage.mode(Z) <- "double"
case_ext <- ph$group == "LUAD"; ctrl_ext <- ph$group == "Control"
mu1 <- rowMeans(Z[, case_ext, drop = FALSE])
mu0 <- rowMeans(Z[, ctrl_ext, drop = FALSE])
v1 <- rowSums((Z[, case_ext, drop = FALSE] - mu1)^2) / (sum(case_ext) - 1)
v0 <- rowSums((Z[, ctrl_ext, drop = FALSE] - mu0)^2) / (sum(ctrl_ext) - 1)
sd <- sqrt(((sum(case_ext) - 1) * v1 + (sum(ctrl_ext) - 1) * v0) /
             (sum(case_ext) + sum(ctrl_ext) - 2))

effect <- data.frame(cpg = M$cpg, delta = mu1 - mu0,
                     magnitude = abs((mu1 - mu0) / sd)) %>%
  inner_join(orientation[, c("cpg", "direction")], by = "cpg") %>%
  mutate(agree = sign(delta) == direction) %>% filter(is.finite(magnitude))
members <- top %>%
  inner_join(selected[, c("stable_pool_rank", "program_index")], by = "stable_pool_rank") %>%
  inner_join(effect, by = "cpg")

concordance <- bind_rows(lapply(1:9, function(j) {
  d <- subset(members, program_index == j)
  bind_rows(lapply(c(100, 30, 10), function(percent) {
    cut <- quantile(d$magnitude, 1 - percent / 100)
    x <- subset(d, magnitude >= cut)
    data.frame(program = paste0("P", j), retained = percent,
               agreement = mean(x$agree))
  }))
}))
concordance$retained <- factor(concordance$retained, levels = c(100, 30, 10),
                               labels = c("All CpGs", "Top 30%", "Top 10%"))

program_cols <- setNames(
  grDevices::colorRampPalette(
    c("#392F72", "#369E9A", "#D2A15C"),
    space = "Lab"
  )(9),
  paste0("P", 1:9)
)

p <- ggplot(
  concordance,
  aes(retained, agreement, group = program, colour = program)
) +
  geom_line(linewidth = 0.2) +
  geom_point(size = 0.7) +
  geom_hline(yintercept = 0.5, linewidth = 0.2, linetype = "22", color = "#c3c3c3") +
  geom_text(
    data = subset(
      concordance,
      program %in% c("P5", "P8") & retained == "Top 10%"
    ),
    aes(label = program),
    hjust = -0.5,
    vjust = 0.5,
    size = 1.9,
    show.legend = FALSE
  ) +
  scale_x_discrete(
    expand = expansion(add = c(0.2, 0.4))
  ) +
  scale_colour_manual(
    values = program_cols,
    breaks = paste0("P", 1:9)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  labs(
    x = "External effect-size subset",
    y = "Directional agreement",
    colour = "Program",
    tag = "c"
  ) +
  basic +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "top",
    legend.box.spacing = unit(1, "mm"),
    legend.title = element_text(
      margin = margin(b = 4)
    ),
    axis.title.y = element_text(margin = margin(r = 0)),
    plot.tag.position = c(-.02, 0.95),
    plot.margin = margin(l = 12, r = 5, b = 5, t = 5),
    
  ) +
  guides(
    colour = guide_legend(ncol = 1)
  )

# Save
ggsave(
  file.path("result/direct-agree.svg"),
  device = svglite,
  p,
  width = 2.3, height = 1.6, units = "in",bg = "transparent"
)
