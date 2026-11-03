# 1. Set Environment ───────────────────────────────────────────────────────────
library(grid)
library(ragg)
library(mgcv)
library(haven)
library(timeROC)
library(survival)
library(survminer)
library(data.table)
library(forestploter)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Set priority
conflicts_prefer(dplyr::lag)

# 2. Generate Functions ────────────────────────────────────────────────────────
# Draw hazard ratio compared to naive non-smokers
draw_hr <- function(fit, new_data, n_sim = 1000, keep = TRUE) {
  require(MASS)
  
  new_data %<>% as.data.frame %>% mutate(
    ever     = as.integer(smk_status != "none"),
    past     = as.integer(smk_status %in% "quit"),
    smk_dur  = case_when(
      ever != 1L ~ 0,
      !is.na(smk_dur) ~ smk_dur,
      !is.na(smk_age) ~ pmax(0, age - smk_age),
      TRUE ~ 0
    ),
    ces_dur  = case_when(
      past != 1L ~ 0,
      !is.na(ces_dur) ~ ces_dur,
      !is.na(ces_age) ~ pmax(0, age - ces_age),
      TRUE ~ 0
    ),
    na_amnt  = as.integer(ever != 0L & is.na(smk_amnt)),
    smk_amnt = coalesce(smk_amnt, 0),
    na_cot   = as.integer(is.na(cotinine)),
    use_cot  = as.integer(!is.na(cotinine)),
    cotinine = coalesce(cotinine, 0),
    mono     = ever * na_amnt,
    bi       = ever * (1 - na_amnt)
  )
  
  # Define naive non-smoker baseline profile
  base <- new_data
  base$smk_dur  <- 0
  base$smk_amnt <- 0
  base$bi       <- 0
  base$mono     <- 0
  base$past     <- 0
  base$ces_dur  <- 0
  base$cotinine <- 0
  base$use_cot  <- 1
  base$na_cot   <- 0
  base$na_amnt  <- 0
  
  # Design matrices
  x_0 <- predict(fit, newdata = base,     type = "lpmatrix")
  x_1 <- predict(fit, newdata = new_data, type = "lpmatrix")
  
  # Simulate coefficients
  coefs <- mvrnorm(n_sim, mu = coef(fit), Sigma = vcov(fit))  # n_sim x p
  b <- t(coefs)
  
  # Absolute LP draws (pure risk score)
  eta <- x_1 %*% b
  
  # Relative HR vs. naive non-smokers
  hr <- exp((x_1 - x_0) %*% b)
  
  agg <- function(mat) {
    data <- data.table(
      mean  = rowMeans(mat, na.rm = TRUE),
      lower = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
    )
    
    if (keep) data$draw <- split(mat, row(mat)) else data$draw <- vector("list", nrow(data))
    
    n_row <- nrow(data)
    vars  <- c("mean", "lower", "upper", "draw")
    
    out <- lapply(
      seq_len(n_row),
      \(i) as.data.frame(data[i, ..vars])
    )
    
    return(out)
  }
  
  rs  <- agg(eta)  # pure risk score
  hrs <- agg(hr)   # HR vs. naive non-smoker
  
  new_data %>% select(sample_name, group) %>% mutate(risk = I(rs), hr = I(hrs))
}

# Compute cumulative concordance
cum_conc <- function(time, score) {
  ord <- order(time)
  
  t <- time[ord]
  s <- score[ord]
  m <- length(t)
  
  c_k <- rep(NA, m)
  
  for (k in 2:m) {
    ti <- t[1:k]
    si <- s[1:k]
    
    conc <- 0
    comp <- 0
    half <- 0
    
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        if (ti[i] == ti[j]) next
        comp <- comp + 1
        
        if (ti[i] < ti[j]) {
          if (si[i] > si[j]) conc <- conc + 1
          else if (si[i] == si[j]) half <- half + 1
        } else {
          if (si[j] > si[i]) conc <- conc + 1
          else if (si[i] == si[j]) half <- half + 1
        }
      }
    }
    
    c_k[k] <- if (comp > 0) (conc + 0.5 * half) / comp else NA
  }
  
  c_k[1] <- NA
  
  tibble(time = t, score = s, c_tau = c_k)
}

# Indent
indent <- function(x) paste0("\u00A0\u00A0", x)

# Fit cox
fit_cox <- function(formula, data) {
  fit <- tryCatch(
    coxph(formula, data = data, ties = "efron"), 
    error = function(e) NULL
  )
  
  tidy(fit, exponentiate = TRUE, conf.int = TRUE)
}

# Pick term
trm <- function(td, name) {
  out <- td %>% filter(term %in% name)
  
  out %>% slice(1) %>% select(estimate, conf.low, conf.high, p.value)
}

# Header row
header <- function(label) {
  tibble(variable = label, n = NA, estimate = NA, lower = NA, upper = NA, p.value = NA, ref = FALSE)
}

# Reference row
ref <- function(label, n) {
  tibble(variable = label, n = n, estimate = NA, lower = NA, upper = NA, p.value = NA, ref = TRUE)
}

# Effect row
eff <- function(label, n, td, name) {
  r <- trm(td, name)
  
  tibble(
    variable = label,
    n        = n,
    estimate = r$estimate[[1]],
    lower    = r$conf.low[[1]],
    upper    = r$conf.high[[1]],
    p.value  = r$p.value[[1]],
    ref      = FALSE
  )
}

# HR format
fmt_h <- function(est, low, up, ref = FALSE) {
  case_when(
    ref        ~ "Reference",
    is.na(est) ~ "",
    TRUE       ~ sprintf("%.2f (%.2f-%.2f)", est, low, up)
  )
}

# Significant format
fmt_p <- function(p, ref = FALSE) {
  case_when(
    ref      ~ "",
    is.na(p) ~ "",
    TRUE     ~ format.pval(p, digits = 3, eps = 1e-3)
  )
}

# 3. Cohort Study ──────────────────────────────────────────────────────────────
library(haven)
ncc <- read_sas("../00_manifest/ncc-cohort.sas7bdat")

dim(ncc)
ncc$Tcode1 %>% table


coh <- file.path(manifest_dir, "cohort.rds") %>% read_rds


ncc %<>% select(
  No, cancer_reg, vial, sex, age, lab = label, sample_name = label, dx = CT_Chest_Low_dose
) %>% unique
ncc %<>% mutate(sex = if_else(sex == 1, "M", "F"))

coh <- coh %<>% drop_na(smk_status)

coh <- left_join(coh, ncc)

# Generate variables
coh %<>% mutate(
  event    = !is.na(date_diagnosis),
  sex      = factor(sex, levels = c("M", "F")),
  entry    = age,
  exit     = age + time,
  ever     = smk_status != "none",
  past     = smk_status %in% "quit",
  current  = smk_status %in% "current",
  smk_dur  = case_when(
    isFALSE(ever) ~ 0, 
    !is.na(smk_dur) ~ smk_dur, 
    !is.na(smk_age) ~ pmax(0, age - smk_age), 
    TRUE ~ 0
  ),
  ces_dur  = case_when(
    isFALSE(past) ~ 0,
    !is.na(ces_dur) ~ ces_dur,
    !is.na(ces_age) ~ pmax(0, age - ces_age),
    TRUE ~ 0
  ),
  na_amnt = ever & is.na(smk_amnt),
  na_cot  = is.na(cotinine),
  na_ct   = is.na(ct),
  use_cot = !na_cot
)

coh %>% pull(event) %>% table
coh %>% filter(event != 1) %>% pull(vial) %>% table

# Add screening criteria eligibility
coh %<>% mutate(
  packyear = (smk_amnt / 20) * smk_dur,
  screen   = case_when(
    is.na(age) | is.na(smk_status) ~ NA,
    age < 54 | age > 74                   ~ FALSE,
    packyear < 30                         ~ FALSE,
    smk_status == "none"                  ~ FALSE,
    is.na(smk_amnt) | is.na(smk_dur) ~ NA,
    smk_status == "current"               ~ TRUE,
    smk_status == "quit" & is.na(ces_dur) ~ NA,
    smk_status == "quit" & ces_dur > 15   ~ FALSE,
    smk_status == "quit" & ces_dur <= 15  ~ TRUE,
    TRUE ~ NA
  )
)

coh %<>% mutate_if(is.logical, as.integer) %>% mutate(across(c(smk_amnt, cotinine), ~ coalesce(.x, 0)))

# Assign 2D surface indicators
coh$mono <- coh$ever * coh$na_amnt; coh$bi <- coh$ever * (1 - coh$na_amnt)

# Fit survival model
fit <- gam(
  time ~
    sex +
    s(age, bs = "cs", k = 10) +
    te(smk_dur, smk_amnt, bs = c("cs", "cs"), k = c(5, 5), by = bi) +
    s(smk_dur, bs = "cs", k = 5, by = mono) +
    s(ces_dur, bs = "cs", k = 5, by = past) +
    s(log1p(cotinine), bs = "cs", k = 5, by = use_cot) +
    na_cot + na_amnt,
  family  = cox.ph(),
  data    = coh,
  weights = event,
  method  = "REML",
  select  = TRUE
)

summary(fit)

# Risk assessment
risk <- predict(fit, newdata = coh, type = "link")

coh %<>% mutate(.pred = risk)

# Generate risk surface
setDT(coh)
brk <- c(1, 1.5, 2, 5)  # steps for risk

age <- median(coh$age, na.rm = TRUE)
cot <- median(coh$cotinine[coh$ever != 0], na.rm = TRUE)

# Set bands
coh[, age_band := fcase(age < 55, "30-54", age >= 55 & age < 60, "55-59", age >= 60, "≥60")]
coh$age_band %<>% fct_relevel("30-54", "55-59", "≥60")

coh[event != 0, time_bin := fcase(time > 10, ">10", time <= 10 & time > 5, "5-10", time <= 5 & time > 2, "2-5", default = "<2")]
coh$time_bin %<>% fct_relevel("<2", "2-5", "5-10", ">10")

# Grid for contours
gr <- with(
  coh, {
    dur  <- max(smk_dur[ever != 0L], na.rm = TRUE)
    amnt <- quantile(smk_amnt[ever != 0L], 0.995, na.rm = TRUE)
    
    expand.grid(
      smk_dur  = seq(0, dur,  by = 1),
      smk_amnt = seq(0, amnt, by = 1)
    )
  }
)

# Delta method for relative risk
beta <- coef(fit)
cov  <- vcov(fit)

n_g <- nrow(gr)

rs <- data.frame(
  time     = 1,  # placeholder
  sex      = factor(rep("M", n_g), levels = c("M", "F")),
  age      = rep(age, n_g),
  smk_dur  = gr$smk_dur,
  smk_amnt = gr$smk_amnt,
  past     = 0,
  ces_dur  = 0,
  cotinine = cot,
  use_cot  = 1,
  na_cot   = 0,
  na_amnt  = 0,
  bi       = 1,
  mono     = 0
)

base <- rs
base$smk_dur  <- 0
base$smk_amnt <- 0
base$bi       <- 0
base$mono     <- 0

x_1 <- predict(fit, newdata = rs,   type = "lpmatrix")
x_0 <- predict(fit, newdata = base, type = "lpmatrix")

x_d    <- x_1 - x_0
log_rr <- as.numeric(x_d %*% beta)

x_v <- x_d %*% cov
se  <- sqrt(pmax(rowSums(x_v * x_d), 0))

gr$mu    <- exp(log_rr)
gr$lower <- exp(log_rr - 1.645 * se)

# Copy grids as much as age bands
surface <- do.call(
  rbind, lapply(
    c("30-54", "55-59", "≥60"), function(b) {
      transform(gr, age_band = b)
    }
  )
)

surface$age_band %<>% fct_relevel("30-54", "55-59", "≥60")

# Extract current smokers
cur <- coh[
  sex == "M" & current == 1L & event != 0L & na_amnt == 0L, 
  .(age_band, smk_dur, smk_amnt, time, time_bin, morphology)
]

cur[, smk_dx := smk_dur + time]; cur[smk_dx < smk_dur, smk_dx := smk_dur]
cur[, amnt   := smk_amnt + runif(.N, 0, 0)]

# Remove incomplete rows
cur %<>% drop_na(morphology)

# Plot RR compared to non-smoker group
brd <- expand_grid(
  packyears = c(20, 30), 
  smk_dur   = seq(0.5, max(surface$smk_dur), by = 0.1)
) %>% mutate(
  label = paste(packyears, "PY"),
  type  = if_else(packyears != 30, "longdash", "dashed"),
  amnt  = 20 * packyears / smk_dur
) %>% filter(
  is.finite(amnt), amnt < max(surface$smk_amnt)
)

ctour <- ggplot(
  surface, aes(smk_dur, smk_amnt)
) + geom_textpath(
  data = brd,
  aes(x = smk_dur, y = amnt, label = label, group = label, linetype = type),
  color = "#a1a1a1",
  size  = 1.8,
  hjust = 0.15,
  vjust = -0.1,
  linewidth   = 0.35,
  show.legend = FALSE,
  inherit.aes = FALSE
) + geom_textcontour(
  aes(z = lower, linewidth = after_stat(if_else(level == 1, 0.25, 0.1)), label = after_stat(paste0("HR≥", level))),
  size = 1.8, breaks = brk, linecolour = "#000000", show.legend = FALSE
) + geom_point(
  data = cur,
  aes(x = smk_dur, y = amnt, color = morphology, size = time_bin),
  shape  = 19,
  stroke = 0
) + scale_colour_npg(
  alpha = 0.8
) + scale_size_manual(
  values = c("<2" = 2.8, "2-5" = 2.2, "5-10" = 1.5, ">10" = 1),
  drop = FALSE
) + scale_linetype_identity(
  guide = "none"
) + scale_linewidth_identity(
  guide = "none"
) + scale_y_continuous(
  expand = expansion(mult = c(0, 0), add = c(0, 0.75))
) + facet_wrap(
  ~ age_band, labeller = labeller(age_band = \(x) paste("Aged", x))
) + labs(
  x = "Smoking duration (years)", 
  y = "Cigarettes per day",
  size   = "Time-to-diagnosis (years)",
  colour = "Morphology"
)

ctour <- ctour + basic + guides(
  size  = guide_legend(order = 1),
  color = guide_legend(order = 2, override.aes = list(size = 2.1))
) + theme(
  axis.title.y = element_text(margin = margin(r = 4)),
  legend.position    = "right",
  legend.box         = "vertical",
  legend.justification = "top",
  legend.margin      = margin(t = -10, 0, 0)
) + labs(
  tag = "a"
)

# Save
ggsave("result/contour.svg", device = svglite, ctour, width = 5.5, height = 1.6, dpi = 800, bg = "transparent")

# Simplify morphology
sim <- coh %>% mutate(
  aged = age / 5,
  morp = case_when(
    event != 1 ~ NA,
    morphology %in% c("Adenocarcinoma (incl. variants)", "Squamous cell carcinoma", "Small cell carcinoma") ~ morphology,
    TRUE ~ "NOS, unspecified"
  ),
  dx
)

sim$smk_status %<>% fct_recode(Never = "none", Quit = "quit", Current = "current") %>% fct_relevel("Never", "Quit")

# Plot Kaplan-Meier
mod <- survfit(Surv(time, event) ~ screen, data = sim)
te  <- survdiff(Surv(time, event) ~ screen, data = sim)

p <- 1 - pchisq(te$chisq, df = 1)

obj <- ggsurvplot(
  fit           = mod,
  data          = sim,
  size          = 0.2,
  alpha         = 0.9,
  conf.int      = FALSE,
  censor        = TRUE,
  censor.size   = 0.85,
  xlab          = "Time-to-diagnosis (years)",
  ylab          = "Survival probability",
  break.time.by = 5,
  risk.table    = FALSE,
  legend.labs   = c("Ineligible", "Eligible"),
  palette       = c("#292473", "#e03131")
)

obj$plot <- obj$plot + labs(
  color = NULL, fill = NULL
) + scale_x_continuous(
  expand    = expansion(0), 
  transform = new_transform(
    name      = "expand",
    transform = function(x) x ^ 2,
    inverse   = function(x) x ^ (1 / 2)
  )
) + scale_y_continuous(
  limits = c(0.68, 1)
) + annotate(
  "text",
  x = 5, y = 0.68,
  label = "'Log-rank'~italic(p)<0.001",
  parse = TRUE,
  hjust = 0,
  vjust = 0,
  size  = 2.15
)

# Plot Kaplan Meier
km <- obj$plot + basic + theme(
  legend.position = "top",
  legend.margin = margin(b = -5, t = 5, 0, 0)
) + labs(tag = "b")

# Save
ggsave("result/kaplan.svg", device = svglite, km, width = 1.8, height = 1.7, dpi = 800, bg = "transparent")

# Explore about risk composed of age, sex and smoking status
exp <- fit_cox(Surv(time, event) ~ screen, sim)

# Age-specific cox
t <- levels(sim$time_bin)

band <- map_dfr(
  t, function(x) {
    data <- sim %>% mutate(ev = event == 1 & time_bin == x)
    fit  <- fit_cox(Surv(time, ev) ~ screen, data)
    
    eff(
      label = indent(x),
      n     = sum(data$ev == 1, na.rm = TRUE),
      td    = fit,
      name  = "screen"
    )
  }
)

# Morphology-specific Cox
type <- c("Adenocarcinoma (incl. variants)", "Small cell carcinoma", "Squamous cell carcinoma", "NOS, unspecified")

sub <- map_dfr(
  type, function(x) {
    data <- sim %>% mutate(ev = event == 1 & morp == x)
    fit  <- fit_cox(Surv(time, ev) ~ screen, data)
    
    eff(
      label = indent(x),
      n     = sum(data$ev == 1, na.rm = TRUE),
      td    = fit,
      name  = "screen"
    )
  }
)

# Count categories
cnts <- list(
  eligible = sum(sim$screen, na.rm = TRUE),
  through  = sum(!sim$screen, na.rm = TRUE),
  never    = sum(sim$smk_status == "Never"),
  quit     = sum(sim$smk_status == "Quit"),
  current  = sum(sim$smk_status == "Current")
)

# Build table
bd <- tribble(
  ~variable,                            ~n,            ~estimate,        ~lower,        ~upper,        ~p.value,        ~ref,
  "Screening eligibility",              NA,            NA,               NA,            NA,            NA,              FALSE,
  indent("Ineligible"),                 cnts$through,  NA,               NA,            NA,            NA,              TRUE,
  indent("Eligible"),                   cnts$eligible, 
  trm(exp, "screen")$estimate[[1]],
  trm(exp, "screen")$conf.low[[1]],
  trm(exp, "screen")$conf.high[[1]],
  trm(exp, "screen")$p.value[[1]],
  FALSE
) %>% bind_rows(
  header("Time-to-diagnosis"),
  band,
  header("Morphology-specific"),
  sub
) %>% mutate(
  n  = if_else(is.na(n), "", comma(n)),
  hr = fmt_h(estimate, lower, upper, ref),
  p  = fmt_p(p.value, ref)
)

# Adjust size
bd$space <- strrep(" ", 15)

# 기존 크기 비율 유지
bd$size <- scales::rescale(
  bd$estimate,
  to = c(1e-8, 0.5)
)

# 너무 작은 점에만 하한 적용
bd$size <- pmax(bd$size, 0.2)

# Forest plot theme
theme <- forest_theme(
  base_size = 6,
  base_family = "Arial",
  
  ci_pch = 15,
  ci_lwd = 0.7,
  
  refline_gp = gpar(
    lty = 2,
    lwd = 0.7
  ),
  
  core = list(
    padding = unit(c(4, 0.9), "mm"),
    bg_params = list(fill = "white", col = NA)
  ),
  
  colhead = list(
    padding = unit(c(4, 0.9), "mm"),
    fg_params = list(fontface = "plain"),
    bg_params = list(fill = "white", col = NA)
  ),
  
  xaxis_gp = gpar(
    lwd = 0.7,
    fontsize = 5.5
  )
)

# Make table
pd <- bd %>% select(variable, n, space, hr, p)
colnames(pd) <- c(
  "Variable",
  "n",
  "", 
  # "Hazard ratio (95% CI)",
  "HR (95% CI)",
  "p-value"
)

fr <- forest(
  pd,
  est       = bd$estimate,
  lower     = bd$lower,
  upper     = bd$upper,
  sizes     = bd$size,
  alpha     = 0.9,
  ci_column = 3,
  ref_line  = 1,
  x_trans   = "log",
  xlim      = c(0.5, 16),
  ticks_at  = c(0.5, 1, 2, 4, 8, 16),
  theme     = theme
)

fr <- edit_plot(
  fr,
  part  = "header",
  col   = c(2, 5),
  which = "text",
  gp    = gpar(fontface = "italic")
)

fr <- add_border(
  fr,
  part = "header",
  row = 1,
  where = "bottom",
  gp = gpar(lwd = 0.8, col = "#000000")
)

forestploter::get_wh(fr, unit = "in")

# SVG vector output
sz <- forestploter::get_wh(fr, unit = "in")

svglite::svglite(
  "result/forest.svg",
  width = sz[1] + 0.1,
  height = sz[2] + 0.1,
  bg = "white"
)

grid::grid.newpage()
grid::grid.draw(fr)
grDevices::dev.off()

# 4. Explore LDCT ──────────────────────────────────────────────────────────────
pt <- sim %>% filter(
  event != 0, !is.na(ct)
) %>% mutate(
  time_grp = if_else(
    time < 3,
    "<3 yrs",
    "≥3 yrs"
  ),
  ct_grp = if_else(
    ct != 0,
    "Nodule present",
    "Normal"
  ),
  time_grp = factor(
    time_grp,
    levels = c("≥3 yrs", "<3 yrs")
  ),
  ct_grp = factor(
    ct_grp,
    levels = c("Normal", "Nodule present")
  ),
  dx
)

pt <- pt[
  grepl("정상|nodule|결절", pt$dx, ignore.case = TRUE),
]

pt$time_grp %<>% fct_relevel("≥3 yrs"); pt$ct_grp %<>% fct_relevel("Negative")

cnts <- pt %>% count(screen, name = "n_tot")

# Organize positions
gap <- 2.4

ctr <- tibble(
  screen = c(0, 1),
  cx = c(0, 0), 
  cy = c(gap, -gap)
)

bands <- tibble(
  lower = c(0.90, 1.65),
  upper = c(1.50, 2.20),
  label = c("Timing", "LDCT")
)

# Inner layer
inn <- pt %>% count(
  screen, time_grp, name = "n"
) %>% left_join(
  cnts, by = "screen"
) %>% group_by(
  screen
) %>% arrange(
  time_grp, .by_group = TRUE
) %>% mutate(
  frac       = n / n_tot,
  end_prop   = cumsum(frac),
  start_prop = lag(end_prop, default = 0),
  start      = pi / 2 + 2 * pi * start_prop,
  end        = pi / 2 + 2 * pi * end_prop,
  r0         = bands$lower[1],
  r          = bands$upper[1]
) %>% ungroup %>% left_join(
  ctr, by = "screen"
) %>% na.omit

# Outer layer
out <- pt %>% count(
  screen, time_grp, ct_grp, name = "n_t"
) %>% left_join(
  inn, by = c("screen", "time_grp")
) %>% rename(
  t_start = start, t_end = end
) %>% group_by(
  screen, time_grp
) %>% arrange(
  ct_grp, .by_group = TRUE
) %>% mutate(
  frac       = n_t / sum(n_t),
  end_prop   = cumsum(frac),
  start_prop = lag(end_prop, default = 0),
  start      = t_start + (t_end - t_start) * start_prop,
  end        = t_start + (t_end - t_start) * end_prop,
  r0         = bands$lower[2],
  r          = bands$upper[2]
) %>% ungroup

# Merge adjacent arcs with identical fill
out_draw <- out %>%
  arrange(screen, start) %>%
  group_by(screen) %>%
  mutate(
    run = cumsum(
      ct_grp != lag(ct_grp, default = first(ct_grp))
    )
  ) %>%
  group_by(screen, run) %>%
  summarise(
    cx = first(cx),
    cy = first(cy),
    r0 = first(r0),
    r = first(r),
    start = min(start),
    end = max(end),
    ct_grp = first(ct_grp),
    .groups = "drop"
  )

# Center label
exp <- cnts %>% left_join(ctr, by = "screen")

exp$screen <- case_when(exp$screen == 0 ~ "Ineligible", exp$screen == 1 ~ "Eligible")
exp$label  <- paste(exp$screen, "participants", sprintf("(%d)", exp$n_tot), sep = "\n")

# Bottom labels
bot <- ctr %>% crossing(bands) %>% mutate(y = cy - (lower + upper) / 2)

# Draw status
donut <- ggplot() + geom_arc_bar(
  data = inn,
  aes(x0 = cx, y0 = cy, r0 = r0, r = r, start = start, end = end, fill = time_grp),
  color = NA
) + scale_fill_manual(
  values = c(
    "≥3 yrs" = "#d3d3d3",
    "<3 yrs" = "#5f9ea0"
  ),
  name = "Diagnosis timing",
  guide = guide_legend(order = 1)
) + new_scale_fill(
) + geom_arc_bar(
  data = out_draw,
  aes(
    x0 = cx, y0 = cy,
    r0 = r0, r = r,
    start = start, end = end,
    fill = ct_grp
  ),
  color = NA
) + scale_fill_manual(
  values = c(
    "Normal" = "#d3d3d3",
    "Nodule present" = "#6c8ead"
  ),
  name = bands$label[2],
  guide = guide_legend(order = 2)
) + new_scale_fill(
) + geom_text(
  data = exp,
  aes(x = cx, y = cy, label = label),
  size = 2,
  lineheight = 1.1
) + geom_label(
  data = bot,
  aes(x = cx, y = y, label = label),
  fill = NA,
  size = 2,
  linewidth = NA
) + coord_fixed(
  clip = "off"
) + basic + theme(
  axis.title       = element_blank(),
  axis.ticks       = element_blank(),
  axis.text        = element_blank(),
  axis.line        = element_blank(),
  legend.spacing.x = unit(-1, "mm"),
  legend.position  = "bottom",
  legend.direction = "vertical",
  legend.box       = "horizontal",
  legend.box.spacing = unit(-2, "mm")
) + labs(
  tag = "c"
)

# Save
ggsave("result/donut.svg", device = svglite, donut, width = 1.6, height = 3.3, dpi = 800, bg = "transparent")

# Generate data for C-index calculation
concor <- tibble(
  time  = coh$time,
  event = coh$event,
  score = risk,
  ever  = coh$ever,
  smk_status = coh$smk_status
) %>% filter(
  event != 0L, 
  !is.na(time), 
  !is.na(score), 
  !is.na(smk_status)
)

concor %<>% mutate(grp = if_else(ever != 0L, "ever", "never"))  # set comparison
concor$grp %<>% fct_relevel("never")

# Compute cumulative concordance
cc <- concor %>% group_by(grp) %>% group_modify(~ cum_conc(.x$time, .x$score))
cc %<>% filter(c_tau > 0)

cc$score %<>% rescale(to = c(0, 1))
bar_scale <- 0.4; cc$y <- cc$score * bar_scale

pnts <- c(2, 5, 10)
brks <- c(0, pnts, Inf)

cc %<>% drop_na(
  c_tau
) %>% mutate(
  bin = cut(time, breaks = brks, include.lowest = TRUE, right = FALSE)
)

marks <- cc %>% group_by(
  grp, bin
) %>% summarise(
  med = median(time), .groups = "drop"
) %>% left_join(
  cc, by = c("grp", "bin")
) %>% group_by(
  grp, bin
) %>% slice_min(
  order_by = abs(time - med), n = 1, with_ties = FALSE
) %>% ungroup %>% arrange(
  grp, time
) %>% mutate(
  label = sprintf("C[tau]==%.3f", c_tau)
)

# Adjust selected label positions
marks <- marks %>% mutate(
  dx = case_when(
    grp == "ever" & sprintf("%.3f", c_tau) == "0.651" ~ 0.8,
    grp == "ever" & sprintf("%.3f", c_tau) == "0.550" ~ -0.7,
    TRUE ~ 0
  ),
  dy = case_when(
    grp == "ever" & sprintf("%.3f", c_tau) == "0.481" ~ 0.06,
    TRUE ~ 0
  )
)

# Adjust arrows
shrink <- 0.25

seg <- marks %>% group_by(
  grp
) %>% arrange(
  time, .by_group = TRUE
) %>% mutate(
  d_t = lead(time),
  d_c = lead(c_tau)
) %>% filter(
  !is.na(d_t), !is.na(d_c), d_t > time
) %>% mutate(
  x    = time + shrink * (d_t - time),
  xend = d_t - shrink * (d_t - time),
  y    = c_tau + shrink * (d_c - c_tau),
  yend = d_c - shrink * (d_c - c_tau)
) %>% ungroup

# Plot concordance
brks <- seq(0, 1, by = 0.5); pal <- c(ever = "#8480c2", never = "#71c9dd")

idx <- ggplot(
  cc, aes(x = time)
) + geom_col(
  aes(y = y, fill = grp),
  width = 0.08,
  alpha = 0.35,
  linewidth = 0,
  position  = "identity",
  show.legend = FALSE
) + geom_point(
  aes(y = c_tau, color = grp), 
  size  = 1.3,
  alpha = 0.1,
  shape = 16,
  na.rm = TRUE
) + geom_segment(
  data = seg,
  aes(x = x, y = y, xend = xend, yend = yend, colour = grp),
  arrow = arrow(type = "closed", length = unit(1, "mm")),
  alpha = 0.9,
  linewidth = 0.25
) + geom_vline(
  xintercept  = pnts,
  linetype    = "dotted",
  linewidth   = 0.35,
  color       = "#c6c6c6",
  show.legend = FALSE
) + geom_point(
  data = marks, 
  aes(x = time, y = c_tau, color = grp),
  size  = 1.8,
  shape = 16
) + geom_text_repel(
  data = marks,
  aes(x = time, y = c_tau, label = label),
  size = 1.9,
  box.padding = 0.2,
  point.padding = 0.2,
  nudge_x = marks$dx,
  nudge_y = marks$dy,
  parse = TRUE,
  seed = 42
) + scale_x_continuous(
  name   = "Time-to-diagnosis (years)",
  expand = expansion(mult = c(0, 0), add = c(0.2, 0.1))
) + scale_y_continuous(
  name = expression(Concordance~(C[tau])),
  limits = c(0, 1),
  breaks = c(0, 0.5, 1),
  labels = c("0.0", "0.5", "1.0"),
  expand = expansion(mult = c(0, 0), add = c(0, 0.05)),
  sec.axis = sec_axis(
    ~ . / 0.4,
    breaks = brks,
    labels = brks,
    name = "Risk score (norm.)",
    guide = guide_axis(cap = "both")
  )
) + facet_wrap(
  ~ grp, labeller = labeller(grp = function(x) paste0("Smoking status: ", x))
) + scale_fill_manual(
  values = pal,
  name   = "Smoking exp."
) + scale_colour_manual(
  values = pal,
  guide  = NULL
)

idx <- idx + basic + theme(
  # plot.margin = margin(l = 15, r = 5, b = 5, t = 5),
  # plot.tag    = element_text(
  #   family = "Arial", face = "bold", size = 16, margin = margin(l = -15)
  # )
) + labs(
  tag = "d"
)

# Save
ggsave("result/concordance.svg", device = svglite, idx, width = 3.6, height = 1.6, dpi = 800, bg = "transparent")

# 5. Estimate Hazard Ratio ─────────────────────────────────────────────────────
hr <- lapply(
  meta, \(x) draw_hr(
    fit = fit, new_data = x
  ) %>% unnest(
    risk
  ) %>% select(
    sample_name, risk_score = mean
  )
)

# Save
write_rds(hr, "object/hazard-risk.rds")
