# 1. Set Environment ───────────────────────────────────────────────────────────
library(sesame)
library(GEOquery)
library(data.table)
library(missMethyl)

options(timeout = Inf)

# Load sources
repo <- Sys.getenv("LUNG_METHYLATION_REPO", unset = "")
if (!nzchar(repo)) stop("Set LUNG_METHYLATION_REPO to the code repository")
source(file.path(repo, "r", "analysis", "analysis_environment.R"))

# Load data
load("object/modeling-set.rda")

anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Set priority
conflicts_prefer(base::union, base::intersect, base::setdiff)

# 2. Data Preprocessing ────────────────────────────────────────────────────────
tsv_path     <- sesameAnno_download("EPICv2.hg38.manifest.tsv.gz")
manifest     <- sesameAnno_buildManifestGRanges(tsv_path)
address_file <- sesameAnno_buildAddressFile(tsv_path)

path    <- find_path("Data")  # find path for .idat files
n_cores <- detectCores() - 1  # number of cores to use

# Run SeSAMe pipeline
beta <- openSesame(
  x             = path,
  prep          = "QCDPB",
  manifest      = address_file,
  platform      = "EPICv2",
  BPPARAM       = SnowParam(n_cores),
  min_beads     = 3,
  collapseToPfx = FALSE
)

# Save beta matrix
write_rds(beta, file.path(data_dir, "beta-raw.rds"))

# 3. Preparation of External Data ──────────────────────────────────────────────
dir.create("external/GSE275371", recursive = TRUE, showWarnings = FALSE)

url  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275371/suppl/GSE275371_Matrix_processed.txt.gz"
dest <- "external/GSE275371/GSE275371_Matrix_processed.txt.gz"

download.file(
  url      = url,
  destfile = dest,
  mode     = "wb",
  method   = "libcurl",
  quiet    = FALSE
)

# Read
ext <- fread(dest, data.table = FALSE)

# Beta columns
nm_idx <- match("Name", names(ext))

spl <- names(ext)[2:(nm_idx - 1)]

beta_ext <- ext[, spl, drop = FALSE] %>% as.matrix
rownames(beta_ext) <- ext$probeID; storage.mode(beta_ext) <- "numeric"

# Phenotype
ph_ext <- tibble(
  sample = colnames(beta_ext)
) %>%  mutate(
  prefix = str_extract(sample, "^[A-Za-z]+"),
  group  = case_when(
    str_detect(sample, "^N[Cc]_")  ~ "Control",
    str_detect(sample, "^(LC|L)_") ~ "LUAD",
    TRUE ~ "Unknown"
  ),
  id = parse_number(sample)
)

# Liftover
anno %<>% as.data.frame %>% rownames_to_column("EPIC_v2") %>% mutate(
  base      = sub("_.*$", "", EPIC_v2),
  EPIC_v1   = na_if(EPICv1_Loci, ""),
  Methyl450 = na_if(Methyl450_Loci, "")
) %>% select(
  EPIC_v2, base, EPIC_v1, Methyl450, chr, pos, Type, ProbeSeqA, ProbeSeqB
)

tgt <- tibble(
  EPIC_v2 = rownames(m_tr),
  order   = seq_along(rownames(m_tr))
) %>% left_join(
  anno, by = "EPIC_v2"
)

ids <- rownames(beta_ext)

map <- tgt %>% mutate(
  ext_cpg = case_when(
    base %in% ids      ~ base,
    EPIC_v1 %in% ids   ~ EPIC_v1,
    Methyl450 %in% ids ~ Methyl450,
    TRUE ~ NA
  ),
  map_type = case_when(
    ext_cpg == base      ~ "base",
    ext_cpg == EPIC_v1   ~ "EPIC_v1",
    ext_cpg == Methyl450 ~ "Methyl450",
    TRUE ~ "unmapped"
  )
)

use <- map %>% drop_na(ext_cpg) %>% arrange(order)

beta_ext <- beta_ext[use$ext_cpg, , drop = FALSE]
rownames(beta_ext) <- use$EPIC_v2

# Transform into M-value
m_ext <- logit(beta_ext)

# Save
save(ph_ext, m_ext, file = "object/external-data.rda")
