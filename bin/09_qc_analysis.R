#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT QC ANALYSIS (BASE R UNISSANT PERFORMANCES ET COMPATIBILITÉ)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 09_qc_analysis.R <summary_tsv> <min_length> <max_length>")
}

summary_tsv_file <- args[1]
min_len_param    <- as.numeric(args[2])
max_len_param    <- as.numeric(args[3])

# Fonction pour extraire les longueurs réelles depuis un FASTQ (via AWK)
get_fastq_lengths <- function(fq_path) {
  if (!file.exists(fq_path) || file.size(fq_path) == 0) return(numeric(0))
  cmd <- paste0("zcat -f '", fq_path, "' | awk 'NR%4==2 {print length($0)}'")
  lens <- as.numeric(system(cmd, intern = TRUE))
  return(lens[!is.na(lens)])
}

# Lecture du résumé de génotypage
summary_dt <- if (file.exists(summary_tsv_file)) {
  read.delim(summary_tsv_file, stringsAsFactors = FALSE)
} else {
  data.frame()
}

# Récupération exclusive des échantillons ayant généré un consensus VALIDE
valid_samples <- c()
if (nrow(summary_dt) > 0 && "status" %in% colnames(summary_dt)) {
  valid_samples <- unique(summary_dt$sample[tolower(summary_dt$status) == "validated"])
}

qc_list <- list()

# Détection de tous les échantillons bruts fusionnés
raw_files <- list.files(".", pattern = ".*_merged\\.fastq(\\.gz)?", full.names = TRUE)

for (rf in raw_files) {
  sid <- gsub("_merged\\.fastq(\\.gz)?", "", basename(rf))
  sid <- gsub("^\\./", "", sid)

  # FILTRE STRICT : Seuls les échantillons validés au génotypage sont traités
  if (!(sid %in% valid_samples)) {
    next
  }

  # 1. STEP 01_RAW
  lens_raw <- get_fastq_lengths(rf)
  if (length(lens_raw) > 0) {
    qc_list[[length(qc_list) + 1]] <- data.frame(
      sample = sid, step = "01_RAW", assignedref = "NONE",
      read_count = length(lens_raw), minlengthread = min(lens_raw),
      maxlengthread = max(lens_raw), meanread_length = round(mean(lens_raw), 1),
      medianread_length = as.character(round(median(lens_raw), 1)),
      pident_blast = "NA", mean_depth_coverage = "NA", clair3_variants = "NA", assigned_reads = "NA",
      stringsAsFactors = FALSE
    )
  }

  # 2. STEP 02_DEHOSTING
  dh_file <- list.files(".", pattern = paste0("^", sid, ".*dehosted.*fastq(\\.gz)?"), full.names = TRUE)[1]
  if (!is.na(dh_file) && file.exists(dh_file)) {
    lens_dh <- get_fastq_lengths(dh_file)
    if (length(lens_dh) > 0) {
      qc_list[[length(qc_list) + 1]] <- data.frame(
        sample = sid, step = "02_DEHOSTING", assignedref = "NONE",
        read_count = length(lens_dh), minlengthread = min(lens_dh),
        maxlengthread = max(lens_dh), meanread_length = round(mean(lens_dh), 1),
        medianread_length = as.character(round(median(lens_dh), 1)),
        pident_blast = "NA", mean_depth_coverage = "NA", clair3_variants = "NA", assigned_reads = "NA",
        stringsAsFactors = FALSE
      )
    }
  }

  # 3. STEP 03_FILTERED_TRIMMED
  tr_file <- list.files(".", pattern = paste0("^", sid, ".*trimmed.*fastq(\\.gz)?"), full.names = TRUE)[1]
  lens_tr <- numeric(0)
  if (!is.na(tr_file) && file.exists(tr_file)) {
    lens_tr <- get_fastq_lengths(tr_file)
    if (length(lens_tr) > 0) {
      qc_list[[length(qc_list) + 1]] <- data.frame(
        sample = sid, step = "03_FILTERED_TRIMMED", assignedref = "NONE",
        read_count = length(lens_tr), minlengthread = min(lens_tr),
        maxlengthread = max(lens_tr), meanread_length = round(mean(lens_tr), 1),
        medianread_length = as.character(round(median(lens_tr), 1)),
        pident_blast = "NA", mean_depth_coverage = "NA", clair3_variants = "NA", assigned_reads = "NA",
        stringsAsFactors = FALSE
      )
    }
  }

  # 4. STEP 05_GENOTYPING
  sample_summary <- summary_dt[summary_dt$sample == sid & tolower(summary_dt$status) == "validated", ]
  
  # Couverture moyenne via Mosdepth
  mosdepth_summary <- list.files(".", pattern = paste0("^", sid, ".*\\.mosdepth\\.summary\\.txt$"), full.names = TRUE)[1]
  mean_cov <- "NA"
  if (!is.na(mosdepth_summary) && file.exists(mosdepth_summary)) {
    ms_data <- read.delim(mosdepth_summary, stringsAsFactors = FALSE)
    total_row <- ms_data[ms_data$chrom == "total", ]
    if (nrow(total_row) > 0) {
      mean_cov <- round(as.numeric(total_row$mean[1]), 2)
    }
  }

  # Nombre de variants réels Clair3
  vcf_file <- list.files(".", pattern = paste0("^", sid, ".*\\.vcf(\\.gz)?$"), full.names = TRUE)[1]
  n_var <- "NA"
  if (!is.na(vcf_file) && file.exists(vcf_file)) {
    vcf_cmd <- paste0("bcftools view -H '", vcf_file, "' | wc -l")
    n_var <- trimws(system(vcf_cmd, intern = TRUE))
  }

  # Longueurs mesurées réelles des reads conservés
  real_mean_len <- if (length(lens_tr) > 0) round(mean(lens_tr), 1) else "NA"
  real_med_len  <- if (length(lens_tr) > 0) as.character(round(median(lens_tr), 1)) else "NA"

  if (nrow(sample_summary) > 0) {
    for (i in 1:nrow(sample_summary)) {
      qc_list[[length(qc_list) + 1]] <- data.frame(
        sample = sid,
        step = "05_GENOTYPING",
        assignedref = as.character(sample_summary$genotype[i]),
        read_count = as.numeric(sample_summary$total_reads[i]),
        minlengthread = min_len_param,
        maxlengthread = max_len_param,
        meanread_length = real_mean_len,
        medianread_length = real_med_len,
        pident_blast = round(as.numeric(sample_summary$pident[i]), 3),
        mean_depth_coverage = mean_cov,
        clair3_variants = n_var,
        assigned_reads = as.numeric(sample_summary$total_reads[i]),
        stringsAsFactors = FALSE
      )
    }
  }
}

# Assemblage et sauvegarde
if (length(qc_list) > 0) {
  final_df <- do.call(rbind, qc_list)
  write.table(final_df, "QC_Metrics_Summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("Aucun échantillon valide à exporter dans le QC.\n")
}