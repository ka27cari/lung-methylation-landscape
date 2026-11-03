# Write the installed R package versions required by the analysis scripts.
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = getwd())
r_dir <- file.path(repo, "r")
files <- list.files(r_dir, pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No R scripts found under ", r_dir)
code <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
patterns <- c("(?:library|require|requireNamespace)\\s*\\(\\s*['\"]?([A-Za-z][A-Za-z0-9.]*)",
              "\\b([A-Za-z][A-Za-z0-9.]*)::")
packages <- character()
for (pattern in patterns) {
  hits <- regmatches(code, gregexpr(pattern, code, perl = TRUE))[[1]]
  if (length(hits)) packages <- c(packages, sub(pattern, "\\1", hits, perl = TRUE))
}
packages <- sort(unique(packages))
installed <- installed.packages()[, c("Package", "Version"), drop = FALSE]
rows <- data.frame(package = packages,
                   installed_version = installed[match(packages, installed[, "Package"]), "Version"],
                   stringsAsFactors = FALSE)
write.csv(rows, file.path(repo, "R_PACKAGE_INVENTORY.csv"), row.names = FALSE, na = "")
cat("Wrote", nrow(rows), "package entries; R", as.character(getRversion()), "\n")
