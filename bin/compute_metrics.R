#!/usr/bin/env Rscript

# Options pour forcer le point comme séparateur décimal
options(OutDec = ".")

argv <- commandArgs(TRUE)

rawtablepath    <- as.character(argv[1])
dehosttablepath <- as.character(argv[2])
trimtablepath   <- as.character(argv[3])
genotablepath   <- as.character(argv[4])
outputtablepath <- as.character(argv[5])

raw_count    <- read.csv2(rawtablepath, header = FALSE, stringsAsFactors = FALSE)
dehost_count <- read.csv2(dehosttablepath, header = FALSE, stringsAsFactors = FALSE)
trimm_count  <- read.csv2(trimtablepath, header = FALSE, stringsAsFactors = FALSE)
geno_count   <- read.csv2(genotablepath, header = FALSE, stringsAsFactors = FALSE)

ALLDATA <- rbind(raw_count, dehost_count, trimm_count, geno_count)

colnames(ALLDATA) <- c("readlength", "sample", "step", "assignedref", "pident", "depth", "variants", "assigned_reads")

ALLDATA$count <- 1
ALLDATA$readlength <- as.numeric(ALLDATA$readlength)

# Agrégation des métriques
countread_data <- aggregate(ALLDATA$count, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), sum)
colnames(countread_data) <- c("sample", "step", "assignedref", "read_count")

minlengthread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), min, na.rm = TRUE)
colnames(minlengthread_data) <- c("sample", "step", "assignedref", "minlengthread")

maxlengthread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), max, na.rm = TRUE)
colnames(maxlengthread_data) <- c("sample", "step", "assignedref", "maxlengthread")

meanread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), mean, na.rm = TRUE)
colnames(meanread_data) <- c("sample", "step", "assignedref", "meanread_length")

medianread_data <- aggregate(ALLDATA$readlength, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), median, na.rm = TRUE)
colnames(medianread_data) <- c("sample", "step", "assignedref", "medianread_length")

pident_data <- aggregate(ALLDATA$pident, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(pident_data) <- c("sample", "step", "assignedref", "pident_blast")

depth_data <- aggregate(ALLDATA$depth, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(depth_data) <- c("sample", "step", "assignedref", "mean_depth_coverage")

variants_data <- aggregate(ALLDATA$variants, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(variants_data) <- c("sample", "step", "assignedref", "clair3_variants")

assigned_data <- aggregate(ALLDATA$assigned_reads, by = list(ALLDATA$sample, ALLDATA$step, ALLDATA$assignedref), function(x) x[1])
colnames(assigned_data) <- c("sample", "step", "assignedref", "assigned_reads")

# Fusion propre des données
METRIC_data <- merge(countread_data, minlengthread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, maxlengthread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, meanread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, medianread_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, pident_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, depth_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, variants_data, by = c("sample", "step", "assignedref"))
METRIC_data <- merge(METRIC_data, assigned_data, by = c("sample", "step", "assignedref"))

# Formatage explicite du nombre flottant
METRIC_data$meanread_length <- sprintf("%.1f", METRIC_data$meanread_length)

# Tri par échantillon et par étape
METRIC_data <- METRIC_data[order(METRIC_data$sample, METRIC_data$step), ]

# Écriture propre avec séparateur de colonnes ';' et point décimal
write.table(METRIC_data, outputtablepath, sep = ";", row.names = FALSE, quote = FALSE)
cat("✅ Tableau métrique réaligné avec succès :", outputtablepath, "\n")