# Export Figure 1c aggregate counts from existing cohort and LDCT reports.
# No model fitting, image rendering, or writes to source data directories.

for (pkg in c("dplyr", "haven", "tidyselect")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package not installed: ", pkg)
}

# Usage: Rscript r/source_data/export_fig1c_counts.R <cohort.rds> <ncc-cohort.sas7bdat> <output.csv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Provide: cohort.rds ncc-cohort.sas7bdat output.csv")
cohort_file <- args[[1]]
sas_file <- args[[2]]
out_file <- args[[3]]
if (!file.exists(cohort_file) || !file.exists(sas_file)) stop("Input file not found")

coh <- readRDS(cohort_file)
required_coh <- c("smk_status", "date_diagnosis", "time", "smk_amnt",
                  "smk_dur", "smk_age", "ces_dur", "ces_age")
missing_coh <- setdiff(required_coh, names(coh))
if (length(missing_coh)) stop("Missing cohort columns: ", paste(missing_coh, collapse = ", "))

sas_names <- names(haven::read_sas(sas_file, n_max = 0))
report_candidates <- c("CT_Chest_Low_dose", "LD_CHEST")
report_col <- intersect(report_candidates, sas_names)
if (!length(report_col)) stop("No LDCT report column (CT_Chest_Low_dose or LD_CHEST) in SAS file")
report_col <- report_col[[1]]
needed_sas <- c("No", "cancer_reg", "vial", "sex", "age", "label", report_col)
missing_sas <- setdiff(needed_sas, sas_names)
if (length(missing_sas)) stop("Missing SAS columns: ", paste(missing_sas, collapse = ", "))

ncc <- haven::read_sas(sas_file, col_select = tidyselect::all_of(needed_sas))
ncc <- dplyr::distinct(dplyr::transmute(
  ncc, No = .data$No, cancer_reg = .data$cancer_reg,
  vial = .data$vial, sex = dplyr::if_else(.data$sex == 1, "M", "F"),
  age = .data$age, lab = .data$label, sample_name = .data$label,
  dx = as.character(.data[[report_col]])
))

# Use the same natural join as the original Figure 1 analysis.
coh <- dplyr::filter(coh, !is.na(.data$smk_status))
join_keys <- intersect(names(coh), names(ncc))
if (!length(join_keys)) stop("No matching cohort/SAS join columns")
message("Join keys: ", paste(join_keys, collapse = ", "))
coh <- dplyr::left_join(coh, ncc, by = join_keys)
if (!"dx" %in% names(coh)) stop("LDCT report was not joined")

# Classify only reports explicitly describing normal findings or nodules.
# Other abnormalities and ambiguous reports are not treated as nodules.
dx <- as.character(coh$dx)
normal <- !is.na(dx) & grepl("정상", dx, ignore.case = TRUE)
nodule <- !is.na(dx) & grepl("nodule|결절", dx, ignore.case = TRUE)
coh$ct <- dplyr::case_when(
  nodule & !normal ~ 1L,
  normal & !nodule ~ 0L,
  TRUE ~ NA_integer_
)

coh <- dplyr::mutate(
  coh,
  event = !is.na(.data$date_diagnosis),
  ever = .data$smk_status != "none",
  past = .data$smk_status %in% "quit",
  smk_dur = dplyr::case_when(
    !.data$ever ~ 0,
    !is.na(.data$smk_dur) ~ .data$smk_dur,
    !is.na(.data$smk_age) ~ pmax(0, .data$age - .data$smk_age),
    TRUE ~ 0
  ),
  ces_dur = dplyr::case_when(
    !.data$past ~ 0,
    !is.na(.data$ces_dur) ~ .data$ces_dur,
    !is.na(.data$ces_age) ~ pmax(0, .data$age - .data$ces_age),
    TRUE ~ 0
  ),
  packyear = (.data$smk_amnt / 20) * .data$smk_dur,
  screen = dplyr::case_when(
    is.na(.data$age) | is.na(.data$smk_status) ~ NA,
    .data$age < 54 | .data$age > 74 ~ FALSE,
    .data$packyear < 30 ~ FALSE,
    .data$smk_status == "none" ~ FALSE,
    is.na(.data$smk_amnt) | is.na(.data$smk_dur) ~ NA,
    .data$smk_status == "current" ~ TRUE,
    .data$smk_status == "quit" & is.na(.data$ces_dur) ~ NA,
    .data$smk_status == "quit" & .data$ces_dur > 15 ~ FALSE,
    .data$smk_status == "quit" & .data$ces_dur <= 15 ~ TRUE,
    TRUE ~ NA
  )
)

pt <- dplyr::filter(coh, .data$event, !is.na(.data$ct), !is.na(.data$time))
pt <- dplyr::mutate(
  pt,
  diagnosis_interval = dplyr::if_else(.data$time < 3, "<3 yrs", ">=3 yrs"),
  ldct_category = dplyr::if_else(.data$ct == 0L, "Normal", "Nodule present"),
  eligibility = dplyr::case_when(
    is.na(.data$screen) ~ "Unknown",
    .data$screen ~ "Eligible",
    TRUE ~ "Ineligible"
  )
)
counts <- dplyr::count(pt, .data$eligibility, .data$diagnosis_interval,
                       .data$ldct_category, name = "n")
counts <- as.data.frame(counts)
counts <- counts[order(counts$eligibility, counts$diagnosis_interval, counts$ldct_category), , drop = FALSE]


utils::write.csv(counts, out_file, row.names = FALSE, fileEncoding = "UTF-8")
message("Output: ", normalizePath(out_file, mustWork = FALSE))
message("Cases with unambiguous normal/nodule reports: ", nrow(pt))
message("Eligible: ", sum(pt$screen %in% TRUE), "; Ineligible: ", sum(pt$screen %in% FALSE),
        "; Unknown: ", sum(is.na(pt$screen)))
message("Excluded ambiguous reports (normal AND nodule terms): ",
        sum(coh$event & normal & nodule, na.rm = TRUE))
message("Review small cells before sharing the counts publicly. Reports are text-classified; verify against the constructed original ct variable.")
