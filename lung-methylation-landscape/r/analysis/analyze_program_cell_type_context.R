suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GEOquery)
  library(Biobase)
  library(glmmTMB)
  library(pbapply)
  library(scales)
})

# =============================================================================
# PATH
# =============================================================================

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

manifest_dir <- file.path(base, "00_manifest")
pheno_dir <- file.path(base, "01_phenotype")
data_dir <- file.path(base, "02_preprocess")
project_dir <- file.path(base, "05_ncc")
build_dir <- file.path(project_dir, "landscape", "rebuild")
landscape_data_dir <- file.path(project_dir, "landscape", "data")
external_dir <- file.path(project_dir, "external")
object_dir <- file.path(project_dir, "object")
cache_dir <- file.path(object_dir, "cache")
source_cache_dir <- file.path(external_dir, "model_source_cache")

dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_cache_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(build_dir)) stop("Missing build directory: ", build_dir, call. = FALSE)
if (!dir.exists(landscape_data_dir)) stop("Missing landscape data directory: ", landscape_data_dir, call. = FALSE)

overwrite_reference <- FALSE
overwrite_projection <- FALSE
overwrite_aov <- FALSE
beta_epsilon <- 1e-6

cat("\nExternal data: ", external_dir, "\n", sep = "")
cat("Analysis cache: ", cache_dir, "\n\n", sep = "")

# =============================================================================
# CELL TYPES
# =============================================================================

cell_labels <- c(
  Neu = "Neutrophil", Mono = "Monocyte", Eos = "Eosinophil", Bas = "Basophil",
  Bcell = "B cell", Bnv = "B naive", Bmem = "B memory",
  CD4T = "CD4 T", CD4nv = "CD4 naive", CD4mem = "CD4 memory",
  CD8T = "CD8 T", CD8nv = "CD8 naive", CD8mem = "CD8 memory",
  Treg = "Regulatory T cell", NK = "NK cell"
)

lineage_map <- c(
  Neu = "Myeloid", Mono = "Myeloid", Eos = "Myeloid", Bas = "Myeloid",
  Bcell = "Lymphoid", Bnv = "Lymphoid", Bmem = "Lymphoid",
  CD4T = "Lymphoid", CD4nv = "Lymphoid", CD4mem = "Lymphoid",
  CD8T = "Lymphoid", CD8nv = "Lymphoid", CD8mem = "Lymphoid",
  Treg = "Lymphoid", NK = "Lymphoid"
)

fine_order <- c("Neu", "Mono", "Eos", "Bas", "Bnv", "Bmem", "CD4nv", "CD4mem", "CD8nv", "CD8mem", "Treg", "NK")
program_cols <- sprintf("program_%02d", 1:9)
program_labels <- paste0("P", 1:9)

harmonize_cell_type <- function(x) {
  key <- as.character(x) %>% str_trim() %>% str_to_lower() %>% str_replace_all("[^a-z0-9]", "")
  case_when(
    key %in% c("neu", "neut", "neutrophil", "neutrophils") ~ "Neu",
    key %in% c("mono", "monocyte", "monocytes") ~ "Mono",
    key %in% c("eos", "eosinophil", "eosinophils") ~ "Eos",
    key %in% c("bas", "baso", "basophil", "basophils") ~ "Bas",
    key %in% c("bcell", "bcells", "bcellbulk") ~ "Bcell",
    key %in% c("bnv", "bn", "bnaive", "bnaivecell", "bnaivecells") ~ "Bnv",
    key %in% c("bmem", "bmemory", "bmemorycell", "bmemorycells") ~ "Bmem",
    key %in% c("cd4t", "cd4", "cd4tcell", "cd4tcells") ~ "CD4T",
    key %in% c("cd4nv", "cd4naive", "cd4tnaive") ~ "CD4nv",
    key %in% c("cd4mem", "cd4memory", "cd4tmemory") ~ "CD4mem",
    key %in% c("cd8t", "cd8", "cd8tcell", "cd8tcells") ~ "CD8T",
    key %in% c("cd8nv", "cd8naive", "cd8tnaive") ~ "CD8nv",
    key %in% c("cd8mem", "cd8memory", "cd8tmemory") ~ "CD8mem",
    key %in% c("treg", "tregs", "tregulatory", "regulatoryt") ~ "Treg",
    key %in% c("nk", "nkcell", "nkcells", "naturalkiller") ~ "NK",
    str_detect(key, "mix") ~ "MIX",
    TRUE ~ NA_character_
  )
}

normalize_array_id <- function(x) {
  x <- basename(as.character(x))
  x <- str_replace(x, "\\.gz$", "")
  x <- str_replace(x, "\\.(idat|IDAT)$", "")
  str_replace(x, "_(Grn|Red)$", "")
}

pick_column <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  hit[[1]]
}

align_beta_to_meta <- function(beta, meta) {
  keep <- intersect(meta$sample_id, colnames(beta))
  meta <- meta %>% filter(sample_id %in% keep)
  beta <- beta[, meta$sample_id, drop = FALSE]
  stopifnot(identical(colnames(beta), meta$sample_id))
  list(beta = beta, meta = meta)
}

p_star <- function(p) case_when(
  is.na(p) ~ "ns",
  p < 0.001 ~ "***",
  p < 0.01 ~ "**",
  p < 0.05 ~ "*",
  TRUE ~ "ns"
)

# =============================================================================
# NPZ READER
# =============================================================================

read_npy <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", 6L)
  if (!identical(magic, c(as.raw(0x93), charToRaw("NUMPY")))) stop("Invalid NPY: ", path)
  ver <- as.integer(readBin(con, "raw", 2L))
  read_uint <- function(n) {
    z <- as.integer(readBin(con, "raw", n))
    sum(z * 256^(seq_len(n) - 1L))
  }
  hlen <- if (ver[1] == 1L) read_uint(2L) else read_uint(4L)
  header <- rawToChar(readBin(con, "raw", hlen))
  extract_one <- function(pattern) regmatches(header, regexec(pattern, header, perl = TRUE))[[1]][2]
  descr <- extract_one("['\"]descr['\"]\\s*:\\s*['\"]([^'\"]+)['\"]")
  fortran <- identical(extract_one("['\"]fortran_order['\"]\\s*:\\s*(True|False)"), "True")
  shape_txt <- extract_one("['\"]shape['\"]\\s*:\\s*\\(([^)]*)\\)")
  shape <- trimws(strsplit(shape_txt, ",", fixed = TRUE)[[1]])
  shape <- as.integer(shape[nzchar(shape)])
  n <- prod(shape)
  
  read_num <- function(size) {
    v <- readBin(con, numeric(), n = n, size = size, endian = "little")
    if (length(shape) == 1L) return(v)
    if (length(shape) == 2L) {
      if (fortran) return(matrix(v, nrow = shape[1], ncol = shape[2]))
      return(matrix(v, nrow = shape[1], ncol = shape[2], byrow = TRUE))
    }
    array(v, dim = shape)
  }
  
  if (descr %in% c("<f4", "=f4", "|f4")) return(read_num(4L))
  if (descr %in% c("<f8", "=f8", "|f8")) return(read_num(8L))
  
  if (grepl("^<U[0-9]+$", descr)) {
    width <- as.integer(sub("^<U", "", descr))
    payload <- readBin(con, "raw", n * width * 4L)
    b <- matrix(as.integer(payload), nrow = 4L)
    cp <- b[1, ] + 256 * b[2, ] + 65536 * b[3, ] + 16777216 * b[4, ]
    chars <- matrix(cp, nrow = width, ncol = n)
    out <- vapply(seq_len(n), function(i) {
      z <- chars[, i]
      z <- z[z != 0L]
      if (!length(z)) "" else intToUtf8(z)
    }, character(1))
    if (length(shape) == 1L) return(out)
    return(array(out, dim = shape))
  }
  stop("Unsupported NPY dtype: ", descr)
}

read_npz <- function(path, members) {
  available <- unzip(path, list = TRUE)$Name
  if (!all(members %in% available)) stop("Missing NPZ members: ", paste(setdiff(members, available), collapse = ", "))
  td <- tempfile("npz_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  unzip(path, files = members, exdir = td)
  setNames(lapply(members, function(x) read_npy(file.path(td, x))), members)
}

runtime_file <- file.path(build_dir, "cache", "model", "MODEL_RUNTIME_OUTPUT.npz")
runtime <- read_npz(runtime_file, c("selected_loadings.npy", "cpg.npy", "scaler_center.npy", "scaler_scale.npy", "source_ids.npy"))

H <- as.matrix(runtime[["selected_loadings.npy"]])
storage.mode(H) <- "double"
model_cpg <- as.character(runtime[["cpg.npy"]])
model_base_cpg <- str_remove(model_cpg, "_[^_]+$")
scaler_center <- as.numeric(runtime[["scaler_center.npy"]])
scaler_scale <- as.numeric(runtime[["scaler_scale.npy"]])
runtime_source_ids <- as.character(runtime[["source_ids.npy"]])

program_loading <- as_tibble(t(H))
names(program_loading) <- program_cols

program_loading <- program_loading %>%
  mutate(source_cpg = model_cpg, base_cpg = model_base_cpg, .before = 1) %>%
  pivot_longer(starts_with("program_"), names_to = "program", values_to = "loading") %>%
  mutate(program_index = as.integer(str_remove(program, "^program_")), program_label = paste0("P", program_index)) %>%
  group_by(program, program_index, program_label, base_cpg) %>%
  summarise(loading = sum(loading, na.rm = TRUE), .groups = "drop") %>%
  group_by(program, program_index, program_label) %>%
  mutate(loading_percentile = percent_rank(loading), core = loading_percentile >= 0.90, highlight = loading_percentile >= 0.98) %>%
  ungroup()

write_csv(program_loading, file.path(cache_dir, "PROGRAM_LOADING_BASE.csv.gz"))
target_cpg <- sort(unique(program_loading$base_cpg))

# =============================================================================
# SOURCE MODEL DATA
# =============================================================================

candidate_dirs <- c(
  file.path(landscape_data_dir, "model", "data"),
  file.path(landscape_data_dir, "model-new", "model", "data"),
  file.path(source_cache_dir, "model", "data")
)

valid <- candidate_dirs[file.exists(file.path(candidate_dirs, "m_tr.csv")) & file.exists(file.path(candidate_dirs, "nested.csv"))]

if (length(valid)) {
  source_dir <- valid[1]
} else {
  model_zip <- file.path(landscape_data_dir, "model-new.zip")
  if (!file.exists(model_zip)) stop("Missing: ", model_zip)
  unzip(model_zip, files = c("model/data/m_tr.csv", "model/data/nested.csv"), exdir = source_cache_dir, overwrite = FALSE)
  source_dir <- file.path(source_cache_dir, "model", "data")
}

m_tr_file <- file.path(source_dir, "m_tr.csv")
nested_file <- file.path(source_dir, "nested.csv")

read_methylation_csv <- function(path, sample_ids) {
  dt <- fread(path)
  cpg <- as.character(dt[[1]])
  cols <- names(dt)[-1]
  cols <- cols[cols %in% sample_ids]
  X <- as.matrix(dt[, ..cols])
  storage.mode(X) <- "double"
  list(X = t(X), ids = cols, cpg = cpg)
}

meta_source <- read_csv(nested_file, show_col_types = FALSE) %>% mutate(sentrix_id = as.character(sentrix_id))
src <- read_methylation_csv(m_tr_file, meta_source$sentrix_id)
meta_source <- meta_source %>% slice(match(src$ids, sentrix_id))

if (all(runtime_source_ids %in% src$ids)) {
  ord <- match(runtime_source_ids, src$ids)
  src$X <- src$X[ord, , drop = FALSE]
  src$ids <- src$ids[ord]
  meta_source <- meta_source[ord, , drop = FALSE]
}

ord_cpg <- match(model_cpg, src$cpg)
if (anyNA(ord_cpg)) stop("Runtime CpGs missing from m_tr.csv: ", sum(is.na(ord_cpg)))
src$X <- src$X[, ord_cpg, drop = FALSE]
src$cpg <- src$cpg[ord_cpg]

med <- apply(src$X, 2, median, na.rm = TRUE)
bad <- !is.finite(src$X)
if (any(bad)) {
  ij <- which(bad, arr.ind = TRUE)
  src$X[bad] <- med[ij[, 2]]
}

matched_pairs <- function(meta) {
  z <- meta %>% mutate(.row = row_number()) %>% group_by(matched_index) %>%
    summarise(ctrl = list(.row[group == "control"]), pred = list(.row[group == "prediagnostic"]), .groups = "drop") %>%
    filter(lengths(ctrl) == 1L, lengths(pred) == 1L)
  list(ctrl = unlist(z$ctrl), pred = unlist(z$pred))
}

source_orientation <- function(X, meta) {
  p <- matched_pairs(meta)
  d <- X[p$pred, , drop = FALSE] - X[p$ctrl, , drop = FALSE]
  m <- apply(d, 2, median, na.rm = TRUE)
  list(direction = ifelse(m >= 0, 1, -1), median_delta = m)
}

row_pct_pos <- function(A) {
  out <- matrix(0, nrow(A), ncol(A))
  for (i in seq_len(nrow(A))) {
    pos <- A[i, ] > 0
    if (any(pos)) out[i, pos] <- rank(A[i, pos], ties.method = "average") / sum(pos)
  }
  out
}

robust_fit <- function(X) {
  center <- apply(X, 2, median, na.rm = TRUE)
  q10 <- apply(X, 2, quantile, probs = 0.10, na.rm = TRUE, names = FALSE)
  q90 <- apply(X, 2, quantile, probs = 0.90, na.rm = TRUE, names = FALSE)
  scale <- q90 - q10
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

robust_apply <- function(X, s) sweep(sweep(X, 2, s$center, "-"), 2, s$scale, "/")

nearest_indices <- function(reference, query, k) {
  k <- min(k, nrow(reference))
  d2 <- outer(rowSums(query^2), rowSums(reference^2), "+") - 2 * tcrossprod(query, reference)
  d2[d2 < 0] <- 0
  z <- t(vapply(seq_len(nrow(d2)), function(i) order(d2[i, ])[seq_len(k)], integer(k)))
  if (is.null(dim(z))) z <- matrix(z, nrow = 1)
  z
}

choose_anchors <- function(X, meta, ori, frac = 0.02) {
  ctrl <- which(meta$group == "control")
  Xc <- X[ctrl, , drop = FALSE]
  m <- apply(Xc, 2, median, na.rm = TRUE)
  mad <- apply(abs(sweep(Xc, 2, m, "-")), 2, median, na.rm = TRUE) + 1e-6
  score <- 0.5 * rank(mad) / length(mad) + 0.5 * rank(abs(ori$median_delta)) / length(mad)
  order(score)[seq_len(max(20L, round(frac * ncol(X))))]
}

fit_cf <- function(X, meta, ori, frac = 0.02, k = 10L) {
  ctrl <- which(meta$group == "control")
  anc <- choose_anchors(X, meta, ori, frac)
  scaler <- robust_fit(X[ctrl, anc, drop = FALSE])
  list(ctrl = ctrl, anc = anc, scaler = scaler, C = robust_apply(X[ctrl, anc, drop = FALSE], scaler), controls = X[ctrl, , drop = FALSE], direction = ori$direction, k = min(k, length(ctrl)))
}

apply_cf <- function(X, cf) {
  Z <- robust_apply(X[, cf$anc, drop = FALSE], cf$scaler)
  nn <- nearest_indices(cf$C, Z, cf$k)
  ref <- t(vapply(seq_len(nrow(nn)), function(i) colMeans(cf$controls[nn[i, ], , drop = FALSE]), numeric(ncol(cf$controls))))
  row_pct_pos(pmax(sweep(X - ref, 2, cf$direction, "*"), 0))
}

source_cf <- function(X, meta, ori, frac = 0.02, k = 10L) {
  ctrl <- which(meta$group == "control")
  anc <- choose_anchors(X, meta, ori, frac)
  scaler <- robust_fit(X[ctrl, anc, drop = FALSE])
  C <- robust_apply(X[ctrl, anc, drop = FALSE], scaler)
  Z <- robust_apply(X[, anc, drop = FALSE], scaler)
  nn <- nearest_indices(C, Z, min(k + 1L, length(ctrl)))
  ctrl_pos <- setNames(seq_along(ctrl), as.character(ctrl))
  ref <- matrix(NA_real_, nrow(X), ncol(X))
  
  for (i in seq_len(nrow(X))) {
    js <- nn[i, ]
    if (as.character(i) %in% names(ctrl_pos)) js <- js[js != ctrl_pos[[as.character(i)]]]
    js <- js[seq_len(min(k, length(js)))]
    ref[i, ] <- colMeans(X[ctrl[js], , drop = FALSE])
  }
  
  row_pct_pos(pmax(sweep(X - ref, 2, ori$direction, "*"), 0))
}

fit_adapter <- function(X, Y) qr.solve(cbind(1, X), Y)
apply_adapter <- function(X, B) cbind(1, X) %*% B

ori_full <- source_orientation(src$X, meta_source)
A_source_full <- source_cf(src$X, meta_source, ori_full)
S_source_full <- A_source_full %*% t(H)

# =============================================================================
# DOWNLOAD / PREPARE GEO REFERENCES
# =============================================================================

download_large_file <- function(url, destfile, expected_bytes = NULL, max_attempts = 8L) {
  if (!requireNamespace("curl", quietly = TRUE)) stop("Package curl is required.")
  file_bytes <- function(x) if (file.exists(x)) as.numeric(file.info(x)$size) else 0
  
  for (attempt in seq_len(max_attempts)) {
    current <- file_bytes(destfile)
    if (!is.null(expected_bytes) && current == expected_bytes) return(invisible(destfile))
    
    message(sprintf("Download attempt %d/%d; existing %.1f MB", attempt, max_attempts, current / 1024^2))
    res <- curl::multi_download(url, destfile, resume = current > 0, progress = TRUE, multi_timeout = Inf, connecttimeout = 60L, followlocation = TRUE)
    
    size <- file_bytes(destfile)
    if (!is.null(expected_bytes) && size == expected_bytes) return(invisible(destfile))
    if (is.null(expected_bytes) && isTRUE(res$success[[1]])) return(invisible(destfile))
  }
  
  stop("Download incomplete: ", destfile)
}

read_gse167998_beta <- function(matrix_file, target_cpg, meta) {
  message("Inspecting GSE167998 processed matrix ...")
  
  hdr <- fread(matrix_file, nrows = 0, check.names = FALSE, showProgress = FALSE)
  cn <- names(hdr)
  probe_col <- cn[1]
  
  extract_array_id <- function(x) {
    x <- basename(x)
    id <- str_extract(x, "[0-9]{10,}_R[0-9]+C[0-9]+")
    ifelse(!is.na(id), id, normalize_array_id(x))
  }
  
  col_map <- tibble(
    source_col = cn[-1],
    array_id = extract_array_id(cn[-1]),
    gsm = str_extract(cn[-1], "GSM[0-9]+"),
    contains_beta = str_detect(str_to_lower(cn[-1]), "beta")
  )
  
  meta_key <- meta %>% select(sample_id, gsm) %>% distinct()
  
  col_map <- col_map %>%
    mutate(
      by_array = if_else(array_id %in% meta_key$sample_id, array_id, NA_character_),
      by_gsm = meta_key$sample_id[match(gsm, meta_key$gsm)],
      matched_sample = coalesce(by_array, by_gsm)
    ) %>%
    filter(!is.na(matched_sample)) %>%
    mutate(priority = case_when(contains_beta ~ 1L, source_col == matched_sample ~ 2L, TRUE ~ 3L)) %>%
    arrange(matched_sample, priority) %>%
    group_by(matched_sample) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(match(matched_sample, meta_key$sample_id))
  
  if (!nrow(col_map)) stop("Could not match any GSE167998 matrix columns to phenotype samples.")
  
  selected_cols <- col_map$source_col
  message("Matched ", length(selected_cols), " purified-cell samples.")
  
  dt <- fread(matrix_file, select = c(probe_col, selected_cols), check.names = FALSE, showProgress = TRUE)
  setnames(dt, probe_col, "probe_id")
  dt <- dt[probe_id %in% target_cpg]
  
  beta <- as.matrix(dt[, ..selected_cols])
  storage.mode(beta) <- "double"
  rownames(beta) <- dt$probe_id
  colnames(beta) <- col_map$matched_sample[match(selected_cols, col_map$source_col)]
  
  keep <- meta$sample_id[meta$sample_id %in% colnames(beta)]
  beta[, keep, drop = FALSE]
}

prepare_gse110554 <- function() {
  gdir <- file.path(external_dir, "GSE110554")
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  rds <- file.path(gdir, "PROGRAM_CPG_REFERENCE.rds")
  if (file.exists(rds) && !overwrite_reference) return(readRDS(rds))
  
  x <- getGEO("GSE110554", GSEMatrix = TRUE, getGPL = FALSE, destdir = gdir)
  if (is.list(x)) x <- x[[1]]
  ph <- pData(x)
  names(ph) <- make.names(names(ph), unique = TRUE)
  beta <- exprs(x)
  
  type_col <- pick_column(ph, c("cell.type.ch1", "cell_type.ch1"), "GSE110554 cell type")
  array_col <- pick_column(ph, c("array_well.ch1", "array.well.ch1"), "GSE110554 array ID")
  
  meta <- tibble(gsm = rownames(ph), raw_type = as.character(ph[[type_col]]), sentrix_id = as.character(ph[[array_col]])) %>%
    mutate(type = harmonize_cell_type(raw_type), lineage = unname(lineage_map[type]), sample_id = gsm, plate = str_replace(sentrix_id, "_.*$", "")) %>%
    filter(!is.na(type), type != "MIX", !is.na(lineage))
  
  beta <- beta[intersect(target_cpg, rownames(beta)), , drop = FALSE]
  z <- align_beta_to_meta(beta, meta)
  out <- list(accession = "GSE110554", beta = z$beta, pheno = z$meta)
  saveRDS(out, rds, compress = "xz")
  out
}

prepare_gse167998 <- function() {
  gdir <- file.path(external_dir, "GSE167998")
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  rds <- file.path(gdir, "PROGRAM_CPG_REFERENCE.rds")
  if (file.exists(rds) && !overwrite_reference) return(readRDS(rds))
  
  x <- getGEO("GSE167998", GSEMatrix = TRUE, getGPL = FALSE, destdir = gdir)
  if (is.list(x)) x <- x[[1]]
  ph <- pData(x)
  names(ph) <- make.names(names(ph), unique = TRUE)
  
  matrix_file <- file.path(gdir, "GSE167998_matrix_processed.txt.gz")
  if (!file.exists(matrix_file)) {
    download_large_file(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167998/suppl/GSE167998_matrix_processed.txt.gz",
      matrix_file,
      expected_bytes = 819745924
    )
  }
  
  type_col <- pick_column(ph, "description", "GSE167998 cell type")
  basename_col <- pick_column(ph, c("basename.ch1", "basename_ch1"), "GSE167998 basename")
  
  meta <- tibble(gsm = rownames(ph), raw_type = as.character(ph[[type_col]]), sentrix_id = basename(as.character(ph[[basename_col]]))) %>%
    mutate(type = harmonize_cell_type(raw_type), lineage = unname(lineage_map[type]), sample_id = sentrix_id, plate = str_replace(sentrix_id, "_.*$", "")) %>%
    filter(!is.na(type), type != "MIX", !is.na(lineage))
  
  beta <- read_gse167998_beta(matrix_file, target_cpg, meta)
  z <- align_beta_to_meta(beta, meta)
  out <- list(accession = "GSE167998", beta = z$beta, pheno = z$meta)
  saveRDS(out, rds, compress = "xz")
  out
}

ref <- list(GSE110554 = prepare_gse110554(), GSE167998 = prepare_gse167998())

# =============================================================================
# PROJECT P1-P9
# =============================================================================

project_reference <- function(ref_obj, gse) {
  beta <- ref_obj$beta
  meta <- ref_obj$pheno
  query_match <- match(model_base_cpg, rownames(beta))
  keep <- which(!is.na(query_match))
  
  Xsrc <- src$X[, keep, drop = FALSE]
  Xq <- t(beta[query_match[keep], , drop = FALSE])
  ori <- list(direction = ori_full$direction[keep], median_delta = ori_full$median_delta[keep])
  
  cf <- fit_cf(Xsrc, meta_source, ori)
  Asrc <- apply_cf(Xsrc, cf)
  Aq <- apply_cf(Xq, cf)
  
  Hm <- H[, keep, drop = FALSE]
  hnorm <- sqrt(rowSums(Hm^2))
  hnorm[hnorm <= 0 | !is.finite(hnorm)] <- 1
  Hm <- Hm / hnorm
  
  Ssrc <- Asrc %*% t(Hm)
  Sq <- Aq %*% t(Hm)
  adapter <- fit_adapter(Ssrc, S_source_full)
  Sraw <- apply_adapter(Sq, adapter)
  
  scale_safe <- scaler_scale
  scale_safe[!is.finite(scale_safe) | abs(scale_safe) < 1e-12] <- 1
  Sz <- sweep(sweep(Sraw, 2, scaler_center, "-"), 2, scale_safe, "/")
  
  colnames(Sraw) <- program_cols
  colnames(Sz) <- program_cols
  
  wide <- bind_cols(meta %>% select(sample_id, plate, type, lineage), as_tibble(Sz)) %>% mutate(gse = gse, .after = sample_id)
  raw_long <- bind_cols(tibble(sample_id = meta$sample_id), as_tibble(Sraw)) %>% pivot_longer(starts_with("program_"), names_to = "program", values_to = "score_raw")
  
  long <- wide %>%
    pivot_longer(starts_with("program_"), names_to = "program", values_to = "score_z") %>%
    mutate(program_index = as.integer(str_remove(program, "^program_")), program_label = paste0("P", program_index)) %>%
    left_join(raw_long, by = c("sample_id", "program"))
  
  list(wide = wide, long = long)
}

score_files <- file.path(cache_dir, c("PROGRAM_MODEL_SCORES_GSE110554.csv", "PROGRAM_MODEL_SCORES_GSE167998.csv"))

if (all(file.exists(score_files)) && !overwrite_projection) {
  score110 <- read_csv(score_files[1], show_col_types = FALSE)
  score167 <- read_csv(score_files[2], show_col_types = FALSE)
  score110_long <- read_csv(file.path(cache_dir, "PROGRAM_MODEL_SCORES_GSE110554_LONG.csv"), show_col_types = FALSE)
  score167_long <- read_csv(file.path(cache_dir, "PROGRAM_MODEL_SCORES_GSE167998_LONG.csv"), show_col_types = FALSE)
} else {
  p110 <- project_reference(ref$GSE110554, "GSE110554")
  p167 <- project_reference(ref$GSE167998, "GSE167998")
  score110 <- p110$wide
  score167 <- p167$wide
  score110_long <- p110$long
  score167_long <- p167$long
  
  write_csv(score110, file.path(cache_dir, "PROGRAM_MODEL_SCORES_GSE110554.csv"))
  write_csv(score167, file.path(cache_dir, "PROGRAM_MODEL_SCORES_GSE167998.csv"))
  write_csv(score110_long, file.path(cache_dir, "PROGRAM_MODEL_SCORES_GSE110554_LONG.csv"))
  write_csv(score167_long, file.path(cache_dir, "PROGRAM_MODEL_SCORES_GSE167998_LONG.csv"))
}

# =============================================================================
# PROGRAM-LEVEL STATISTICS
# =============================================================================

fit_type_test <- function(df, score_col) {
  d <- df %>% transmute(y = .data[[score_col]], type = factor(type), plate = factor(plate))
  f0 <- lm(y ~ plate, d)
  f1 <- lm(y ~ plate + type, d)
  a <- anova(f0, f1)
  tibble(p_value = a$`Pr(>F)`[2], partial_eta2 = max(0, 1 - sum(residuals(f1)^2) / sum(residuals(f0)^2)))
}

type_tests <- imap_dfr(list(GSE110554 = score110, GSE167998 = score167), function(df, gse) {
  map_dfr(seq_along(program_cols), function(i) fit_type_test(df, program_cols[i]) %>% mutate(gse = gse, program = program_cols[i], program_label = program_labels[i], .before = 1))
}) %>%
  group_by(gse) %>%
  mutate(p_bonf = p.adjust(p_value, "bonferroni"), sig_bonf = p_star(p_bonf)) %>%
  ungroup()

fit_binary <- function(df, score_col, group_col, levels2) {
  d <- df %>% transmute(y = .data[[score_col]], group = factor(.data[[group_col]], levels = levels2), plate = factor(plate)) %>% filter(!is.na(group))
  f0 <- lm(y ~ plate, d)
  f1 <- lm(y ~ plate + group, d)
  a <- anova(f0, f1)
  coef_name <- grep("^group", names(coef(f1)), value = TRUE)[1]
  tibble(effect = unname(coef(f1)[coef_name]), p_value = a$`Pr(>F)`[2])
}

fine_lineage <- score167 %>%
  mutate(lineage2 = case_when(type %in% c("Neu", "Mono", "Eos", "Bas") ~ "Myeloid", type %in% c("Bnv", "Bmem", "CD4nv", "CD4mem", "CD8nv", "CD8mem", "Treg", "NK") ~ "Lymphoid"))

lineage_tests <- map_dfr(seq_along(program_cols), function(i) {
  fit_binary(fine_lineage, program_cols[i], "lineage2", c("Myeloid", "Lymphoid")) %>%
    mutate(program = program_cols[i], program_label = program_labels[i], .before = 1)
}) %>%
  mutate(p_bonf = p.adjust(p_value, "bonferroni"), sig_bonf = p_star(p_bonf))

memory_base <- score167 %>%
  mutate(
    family = case_when(type %in% c("Bnv", "Bmem") ~ "B", type %in% c("CD4nv", "CD4mem") ~ "CD4", type %in% c("CD8nv", "CD8mem") ~ "CD8"),
    state = case_when(type %in% c("Bnv", "CD4nv", "CD8nv") ~ "naive", type %in% c("Bmem", "CD4mem", "CD8mem") ~ "memory")
  ) %>%
  filter(!is.na(family), !is.na(state)) %>%
  mutate(family = factor(family, levels = c("B", "CD4", "CD8")), state = factor(state, levels = c("naive", "memory")))

fit_memory <- function(df, score_col) {
  d <- df %>% transmute(y = .data[[score_col]], family, state, plate = factor(plate))
  f0 <- lm(y ~ plate + family, d)
  f1 <- lm(y ~ plate + family + state, d)
  a <- anova(f0, f1)
  tibble(effect = unname(coef(f1)["statememory"]), p_value = a$`Pr(>F)`[2])
}

memory_pooled <- map_dfr(seq_along(program_cols), function(i) {
  fit_memory(memory_base, program_cols[i]) %>% mutate(program = program_cols[i], program_label = program_labels[i], .before = 1)
}) %>%
  mutate(p_bonf = p.adjust(p_value, "bonferroni"), sig_bonf = p_star(p_bonf))

memory_pairs <- map_dfr(c("B", "CD4", "CD8"), function(fam) {
  d <- memory_base %>% filter(family == fam)
  map_dfr(seq_along(program_cols), function(i) {
    fit_binary(d, program_cols[i], "state", c("naive", "memory")) %>% mutate(family = fam, program = program_cols[i], program_label = program_labels[i], .before = 1)
  })
}) %>%
  mutate(p_bonf_global = p.adjust(p_value, "bonferroni")) %>%
  group_by(program_label) %>%
  mutate(p_holm_program = p.adjust(p_value, "holm"), sig_holm_program = p_star(p_holm_program)) %>%
  ungroup()

fit_one_vs_rest <- function(df, score_col, target) {
  d <- df %>% transmute(y = .data[[score_col]], group = factor(if_else(type == target, "target", "rest"), levels = c("rest", "target")), plate = factor(plate))
  f0 <- lm(y ~ plate, d)
  f1 <- lm(y ~ plate + group, d)
  a <- anova(f0, f1)
  tibble(effect = unname(coef(f1)["grouptarget"]), p_value = a$`Pr(>F)`[2])
}

one_vs_rest <- map_dfr(seq_along(program_cols), function(i) {
  map_dfr(fine_order, function(ct) {
    fit_one_vs_rest(score167, program_cols[i], ct) %>% mutate(program = program_cols[i], program_label = program_labels[i], type = ct, .before = 1)
  })
}) %>%
  mutate(p_bonf_global = p.adjust(p_value, "bonferroni")) %>%
  group_by(program_label) %>%
  mutate(p_holm_program = p.adjust(p_value, "holm"), sig_holm_program = p_star(p_holm_program)) %>%
  ungroup()

strongest <- one_vs_rest %>%
  group_by(program_label) %>%
  slice_max(abs(effect), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(program_label, driver_type = type, driver_effect = effect, driver_p = p_holm_program, driver_sig = sig_holm_program)

driver_summary <- type_tests %>%
  filter(gse == "GSE167998") %>%
  select(program_label, type_partial_eta2 = partial_eta2, type_p_bonf = p_bonf, type_sig = sig_bonf) %>%
  left_join(lineage_tests %>% select(program_label, lineage_effect = effect, lineage_p_bonf = p_bonf, lineage_sig = sig_bonf), by = "program_label") %>%
  left_join(memory_pooled %>% select(program_label, memory_effect = effect, memory_p_bonf = p_bonf, memory_sig = sig_bonf), by = "program_label") %>%
  left_join(strongest, by = "program_label")

write_csv(type_tests, file.path(cache_dir, "PROGRAM_CELLTYPE_TESTS.csv"))
write_csv(lineage_tests, file.path(cache_dir, "PROGRAM_LINEAGE_TESTS.csv"))
write_csv(memory_pooled, file.path(cache_dir, "PROGRAM_MEMORY_POOLED_TESTS.csv"))
write_csv(memory_pairs, file.path(cache_dir, "PROGRAM_MEMORY_PAIR_TESTS.csv"))
write_csv(one_vs_rest, file.path(cache_dir, "PROGRAM_ONE_VS_REST.csv"))
write_csv(driver_summary, file.path(cache_dir, "PROGRAM_DRIVER_SUMMARY.csv"))

# =============================================================================
# CpG BETA-AOV FOR CLOUD
# This is the slow step. It is required once for the cloud and is cached as RDS.
# =============================================================================

fit_beta_lrt_one <- function(i, beta, meta) {
  d <- meta %>% transmute(y = as.numeric(beta[i, ]), type = factor(type), plate = factor(plate)) %>% filter(is.finite(y), !is.na(type), !is.na(plate))
  d$y <- pmin(pmax(d$y, beta_epsilon), 1 - beta_epsilon)
  
  f0 <- try(glmmTMB(y ~ 1 + (1 | plate), dispformula = ~1, family = beta_family(link = "logit"), data = d), silent = TRUE)
  f1 <- try(glmmTMB(y ~ type + (1 | plate), dispformula = ~type, family = beta_family(link = "logit"), data = d), silent = TRUE)
  
  if (inherits(f0, "try-error") || inherits(f1, "try-error")) return(tibble(base_cpg = rownames(beta)[i], p_type = NA_real_))
  a <- try(anova(f0, f1), silent = TRUE)
  
  tibble(base_cpg = rownames(beta)[i], p_type = if (inherits(a, "try-error")) NA_real_ else as.numeric(a$`Pr(>Chisq)`[2]))
}

run_beta_aov <- function(ref_obj, gse) {
  f <- file.path(cache_dir, paste0("CPG_BETA_AOV_", gse, ".rds"))
  if (file.exists(f) && !overwrite_aov) return(readRDS(f))
  
  message("Running CpG beta-AOV: ", gse, " (", nrow(ref_obj$beta), " CpGs)")
  out <- pblapply(seq_len(nrow(ref_obj$beta)), function(i) fit_beta_lrt_one(i, ref_obj$beta, ref_obj$pheno)) %>%
    bind_rows() %>%
    mutate(gse = gse, .before = 1)
  
  saveRDS(out, f, compress = "xz")
  out
}

aov_all <- bind_rows(run_beta_aov(ref$GSE110554, "GSE110554"), run_beta_aov(ref$GSE167998, "GSE167998")) %>%
  group_by(gse) %>%
  mutate(
    n_test = sum(is.finite(p_type)),
    logp_type = -log10(pmax(p_type, .Machine$double.xmin)),
    bonf_threshold = -log10(0.05 / n_test),
    significant_bonf = is.finite(p_type) & p_type < 0.05 / n_test
  ) %>%
  ungroup()

thresholds <- aov_all %>% group_by(gse) %>% summarise(n_test = first(n_test), bonf_threshold = first(bonf_threshold), .groups = "drop")

typewise <- aov_all %>%
  select(gse, base_cpg, p_type, logp_type, significant_bonf) %>%
  pivot_wider(names_from = gse, values_from = c(p_type, logp_type, significant_bonf), names_sep = "__") %>%
  filter(is.finite(p_type__GSE110554), is.finite(p_type__GSE167998)) %>%
  mutate(
    cell_group = case_when(
      significant_bonf__GSE110554 & significant_bonf__GSE167998 ~ "Both",
      xor(significant_bonf__GSE110554, significant_bonf__GSE167998) ~ "Least",
      TRUE ~ "None"
    )
  )

cloud_df <- program_loading %>%
  inner_join(typewise, by = "base_cpg") %>%
  mutate(broad = logp_type__GSE110554, fine = logp_type__GSE167998)

write_csv(aov_all, file.path(cache_dir, "CPG_BETA_AOV_ALL.csv.gz"))
write_csv(thresholds, file.path(cache_dir, "CPG_BONFERRONI_THRESHOLDS.csv"))
write_csv(cloud_df, file.path(cache_dir, "PROGRAM_CPG_CLOUD.csv.gz"))

# =============================================================================
# CLOUD MODE
# =============================================================================

cloud_trans <- pseudo_log_trans(sigma = 1, base = 10)

find_2d_mode <- function(df, adjust = 0.45, grid_n = 180L) {
  x <- cloud_trans$transform(df$broad)
  y <- cloud_trans$transform(df$fine)
  hx <- bw.nrd0(x) * adjust
  hy <- bw.nrd0(y) * adjust
  if (!is.finite(hx) || hx <= 0) hx <- 0.15
  if (!is.finite(hy) || hy <= 0) hy <- 0.15
  gx <- seq(min(x) - 1.5 * hx, max(x) + 1.5 * hx, length.out = grid_n)
  gy <- seq(min(y) - 1.5 * hy, max(y) + 1.5 * hy, length.out = grid_n)
  z <- exp(-0.5 * (outer(gx, x, "-") / hx)^2) %*% t(exp(-0.5 * (outer(gy, y, "-") / hy)^2))
  ij <- arrayInd(which.max(z), dim(z))
  tibble(mode_broad = cloud_trans$inverse(gx[ij[1]]), mode_fine = cloud_trans$inverse(gy[ij[2]]))
}

program_modes <- cloud_df %>%
  filter(core, is.finite(broad), is.finite(fine)) %>%
  group_by(program_label) %>%
  group_modify(~ find_2d_mode(.x)) %>%
  ungroup() %>%
  mutate(program_index = as.integer(str_remove(program_label, "^P"))) %>%
  arrange(program_index)

write_csv(program_modes, file.path(cache_dir, "PROGRAM_CLOUD_MODES.csv"))

# =============================================================================
# PLATE-ADJUSTED BOXPLOT DATA
# =============================================================================

plate_adjust <- function(d) {
  if (n_distinct(d$plate) <= 1L) {
    d$score_plot <- d$score_z
    return(d)
  }
  fit <- lm(score_z ~ factor(plate), d)
  d$score_plot <- residuals(fit) + mean(d$score_z, na.rm = TRUE)
  d
}

box_data <- score167_long %>%
  filter(type %in% fine_order) %>%
  group_by(program, program_label) %>%
  group_modify(~ plate_adjust(.x)) %>%
  ungroup()

write_csv(box_data, file.path(cache_dir, "PROGRAM_BOXPLOT_DATA.csv"))

# =============================================================================
# RELATIVE P1-P9 PCA / BIPLOT DATA
# =============================================================================

cell_info <- tibble(
  type = fine_order,
  cell_name = unname(cell_labels[fine_order]),
  cell_group = c(rep("Myeloid", 4), "Naive", "Memory", "Naive", "Memory", "Naive", "Memory", "Treg", "NK")
)

score167_meta <- score167 %>% left_join(cell_info, by = "type")
Xraw <- score167_meta %>% select(all_of(program_cols)) %>% as.matrix()
raw_pca <- prcomp(Xraw, center = TRUE, scale. = FALSE)
raw_pc1 <- 100 * raw_pca$sdev[1]^2 / sum(raw_pca$sdev^2)

Xrel <- sweep(Xraw, 1, rowMeans(Xraw, na.rm = TRUE), "-")
pca <- prcomp(Xrel, center = TRUE, scale. = FALSE)

sample_pc <- as_tibble(pca$x[, 1:2, drop = FALSE])
names(sample_pc) <- c("PC1", "PC2")
sample_pc <- bind_cols(sample_pc, score167_meta %>% select(sample_id, type, cell_name, cell_group))

centroids <- sample_pc %>%
  group_by(type, cell_name, cell_group) %>%
  summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop")

loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE]) %>% rownames_to_column("program")
names(loadings)[2:3] <- c("PC1", "PC2")

loadings <- loadings %>%
  mutate(program_index = as.integer(str_remove(program, "^program_")), program_label = paste0("P", program_index))

arrow_multiplier <- 0.72 * min(
  diff(range(sample_pc$PC1)) / diff(range(loadings$PC1)),
  diff(range(sample_pc$PC2)) / diff(range(loadings$PC2))
)

loadings <- loadings %>% mutate(PC1_plot = PC1 * arrow_multiplier, PC2_plot = PC2 * arrow_multiplier)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

write_csv(sample_pc, file.path(cache_dir, "PROGRAM_BIPLOT_SAMPLES.csv"))
write_csv(centroids, file.path(cache_dir, "PROGRAM_BIPLOT_CENTROIDS.csv"))
write_csv(loadings, file.path(cache_dir, "PROGRAM_BIPLOT_LOADINGS.csv"))
write_csv(tibble(PC1 = variance[1], PC2 = variance[2], raw_PC1_common_factor = raw_pc1), file.path(cache_dir, "PROGRAM_BIPLOT_VARIANCE.csv"))

cat("\nDONE\nExternal files: ", external_dir, "\nAll analysis results: ", cache_dir, "\n", sep = "")