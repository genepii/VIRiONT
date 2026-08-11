#!/usr/bin/env Rscript

options(OutDec = ".")

argv <- commandArgs(TRUE)

if (length(argv) < 5) {
  stop("Usage: Rscript 11_compute_metrics.R <raw_csv> <dehost_csv> <trimm_csv> <geno_csv> <output_csv>")
}

rawtablepath    <- as.character(argv[1])
dehosttablepath <- as.character(argv[2])
trimtablepath   <- as.character(argv[3])
genotablepath   <- as.character(argv[4])
outputtablepath <- as.character(argv[5])

# Fonction sécurisée pour lire un fichier CSV2 même s'il est vide
read_safe <- function(filepath) {
  if (!file.exists(filepath) || file.info(filepath)$size == 0) {
    return(data.frame())
  }
  df <- tryCatch({
    read.csv2(filepath, header = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
    return(data.frame())
  })
  return(df)
}

raw_count    <- read_safe(rawtablepath)
dehost_count <- read_safe(dehosttablepath)
trimm_count  <- read_safe(trimtablepath)
geno_count   <- read_safe(genotablepath)

# Fonction pour standardiser chaque dataframe à exactement 8 colonnes (évite les erreurs de rbind)
standardize_df <- function(df) {
  if (nrow(df) == 0) return(data.frame())
  ncol_current <- ncol(df)
  if (ncol_current < 8) {
    for (i in (ncol_current + 1):8) {
      df[, i] <- NA
    }
  } else if (ncol_current > 8) {
    df <- df[, 1:8]
  }
  colnames(df) <- c("readlength", "sample", "step", "assignedref", "pident", "depth", "variants", "assigned_reads")
  return(df)
}

raw_count    <- standardize_df(raw_count)
dehost_count <- standardize_df(dehost_count)
trimm_count  <- standardize_df(trimm_count)
geno_count   <- standardize_df(geno_count)

ALLDATA <- rbind(raw_count, dehost_count, trimm_count, geno_count)

if (nrow(ALLDATA) == 0) {
  cat("⚠️ Aucune donnée métrologique à traiter.\n")
  empty_df <- data.frame(
    sample = character(), step = character(), assignedref = character(),
    read_count = numeric(), minlengthread = numeric(), maxlengthread = numeric(),
    meanread_length = numeric(), medianread_length = numeric(),
    pident_blast = character(), mean_depth_coverage = character(),
    clair3_variants = character(), assigned_reads = character()
  )
  write.table(empty_df, outputtablepath, sep = ";", row.names = FALSE, quote = FALSE)
  quit(save = "no", status = 0)
}

# =====================================================================
# 🧹 NETTOYAGE STRICT DES LIGNES PARASITES ET VIDES
# =====================================================================
ALLDATA <- ALLDATA[!is.na(ALLDATA$sample) & trimws(ALLDATA$sample) != "" & grepl("^barcode_", ALLDATA$sample), ]

if (nrow(ALLDATA) == 0) {
  cat("⚠️ Aucune ligne échantillon valide trouvée.\n")
  empty_df <- data.frame(
    sample = character(), step = character(), assignedref = character(),
    read_count = numeric(), minlengthread = numeric(), maxlengthread = numeric(),
    meanread_length = numeric(), medianread_length = numeric(),
    pident_blast = character(), mean_depth_coverage = character(),
    clair3_variants = character(), assigned_reads = character()
  )
  write.table(empty_df, outputtablepath, sep = ";", row.names = FALSE, quote = FALSE)
  quit(save = "no", status = 0)
}

ALLDATA$count <- 1
ALLDATA$readlength <- as.numeric(ALLDATA$readlength)

# --- Agrégation des métriques ---
countread_data <- aggregate(ALLDATA$count, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), sum)
colnames(countread_data) <- c("sample", "step", "assignedref", "read_count")

minlengthread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0) else return(min(x))
})
colnames(minlengthread_data) <- c("sample", "step", "assignedref", "minlengthread")

maxlengthread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0) else return(max(x))
})
colnames(maxlengthread_data) <- c("sample", "step", "assignedref", "maxlengthread")

meanread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0.0) else return(mean(x))
})
colnames(meanread_data) <- c("sample", "step", "assignedref", "meanread_length")

medianread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0.0) else return(median(x))
})
colnames(medianread_data) <- c("sample", "step", "assignedref", "medianread_length")

pident_data   <- aggregate(ALLDATA$pident, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(pident_data) <- c("sample", "step", "assignedref", "pident_blast")

depth_data    <- aggregate(ALLDATA$depth, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(depth_data) <- c("sample", "step", "assignedref", "mean_depth_coverage")

variants_data <- aggregate(ALLDATA$variants, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(variants_data) <- c("sample", "step", "assignedref", "clair3_variants")

assigned_data <- aggregate(ALLDATA$assigned_reads, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(assigned_data) <- c("sample", "step", "assignedref", "assigned_reads")

# --- Fusion propre ---
METRIC_data <- merge(countread_data, minlengthread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, maxlengthread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, meanread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, medianread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, pident_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, depth_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, variants_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, assigned_data, by = c("sample", "step", "assignedref"))

# Formatage explicite de la moyenne
METRIC_data$meanread_length <- sprintf("%.1f", as.numeric(METRIC_data$meanread_length))

# Tri
METRIC_data <- METRIC_data[order(METRIC_data$sample, METRIC_data$step), ]

# Réordonnancement final des colonnes
METRIC_data <- METRIC_data[, c(
  "sample", "step", "assignedref", "read_count", 
  "minlengthread", "maxlengthread", "meanread_length", "medianread_length", 
  "pident_blast", "mean_depth_coverage", "clair3_variants", "assigned_reads"
)]

write.table(METRIC_data, outputtablepath, sep = ";", row.names = FALSE, quote = FALSE)
cat("✅ Tableau métrique réaligné avec succès :", outputtablepath, "\n")