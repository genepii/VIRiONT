#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT QC ANALYSIS (CALCUL DE PROFONDEUR BASÉ SUR LE CONSENSUS VALIDÉ)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 09_qc_analysis.R <summary_tsv> <min_length> <max_length>")
}

summary_tsv_file <- args[1]
min_len_param    <- as.numeric(args[2])
max_len_param    <- as.numeric(args[3])

# Extraction des longueurs réelles des reads depuis un FASTQ (via AWK)
get_fastq_lengths <- function(fq_path) {
  if (!file.exists(fq_path) || file.size(fq_path) == 0) return(numeric(0))
  cmd <- paste0("zcat -f '", fq_path, "' | awk 'NR%4==2 {print length($0)}'")
  lens <- as.numeric(system(cmd, intern = TRUE))
  return(lens[!is.na(lens)])
}

# 1. Lecture du résumé de génotypage issu de 08_genotyping
summary_dt <- if (file.exists(summary_tsv_file)) {
  read.delim(summary_tsv_file, stringsAsFactors = FALSE)
} else {
  data.frame()
}

# Récupération exclusive des échantillons validés
valid_samples <- c()
if (nrow(summary_dt) > 0 && "status" %in% colnames(summary_dt)) {
  valid_samples <- unique(summary_dt$sample[tolower(summary_dt$status) == "validated"])
}

qc_list <- list()

# Détection de tous les fichiers FASTQ bruts fusionnés
raw_files <- list.files(".", pattern = ".*_merged\\.fastq(\\.gz)?", full.names = TRUE)

for (rf in raw_files) {
  sid <- gsub("_merged\\.fastq(\\.gz)?", "", basename(rf))
  sid <- gsub("^\\./", "", sid)

  # FILTRE : Seuls les échantillons validés au génotypage sont inclus
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
  
  mean_read_len_num <- if (length(lens_tr) > 0) mean(lens_tr) else NA
  real_mean_len_str <- if (!is.na(mean_read_len_num)) round(mean_read_len_num, 1) else "NA"
  real_med_len_str  <- if (length(lens_tr) > 0) as.character(round(median(lens_tr), 1)) else "NA"

  if (nrow(sample_summary) > 0) {
    for (i in 1:nrow(sample_summary)) {
      geno          <- as.character(sample_summary$genotype[i])
      reads_geno    <- as.numeric(sample_summary$total_reads[i])
      consensus_len <- as.numeric(sample_summary$length[i])  # <-- TAILLE EXACTE DU GÉNOTE DU TABLEAU
      
      # Calcul de la profondeur réelle basée sur la taille du génome validé
      mean_cov_val <- "NA"
      if (!is.na(consensus_len) && consensus_len > 0 && !is.na(reads_geno) && !is.na(mean_read_len_num)) {
        mean_cov_val <- round((reads_geno * mean_read_len_num) / consensus_len, 2)
      }

      # Comptage des variants Clair3
      vcf_file <- list.files(".", pattern = paste0("^", sid, ".*", geno, ".*\\.vcf(\\.gz)?$"), full.names = TRUE)[1]
      if (is.na(vcf_file)) {
        vcf_file <- list.files(".", pattern = paste0("^", sid, ".*\\.vcf(\\.gz)?$"), full.names = TRUE)[1]
      }
      n_var <- "NA"
      if (!is.na(vcf_file) && file.exists(vcf_file)) {
        vcf_cmd <- paste0("bcftools view -H '", vcf_file, "' 2>/dev/null | wc -l")
        n_var <- trimws(system(vcf_cmd, intern = TRUE))
      }

      qc_list[[length(qc_list) + 1]] <- data.frame(
        sample = sid,
        step = "05_GENOTYPING",
        assignedref = geno,
        read_count = reads_geno,
        minlengthread = min_len_param,
        maxlengthread = max_len_param,
        meanread_length = real_mean_len_str,
        medianread_length = real_med_len_str,
        pident_blast = round(as.numeric(sample_summary$pident[i]), 3),
        mean_depth_coverage = mean_cov_val,
        clair3_variants = n_var,
        assigned_reads = reads_geno,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Export de la table récapitulative
output_file <- "RUN_METRICS_SUMMARY_TABLE.tsv"
qc_columns <- c(
  "sample", "step", "assignedref", "read_count", "minlengthread",
  "maxlengthread", "meanread_length", "medianread_length",
  "pident_blast", "mean_depth_coverage", "clair3_variants", "assigned_reads"
)

if (length(qc_list) > 0) {
  final_df <- do.call(rbind, qc_list)
} else {
  final_df <- data.frame(matrix(ncol = length(qc_columns), nrow = 0))
  colnames(final_df) <- qc_columns
}

write.table(final_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)