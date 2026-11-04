#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    input_csv = NULL,
    min_freq_minor = 8.0,
    sample_id = "sample",
    out_csv = ""
  )
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--input" && i < length(args)) {
      params$input_csv <- args[i + 1]; i <- i + 2
    } else if (arg == "--min-freq-minor" && i < length(args)) {
      params$min_freq_minor <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--sample-id" && i < length(args)) {
      params$sample_id <- args[i + 1]; i <- i + 2
    } else if (arg == "--out" && i < length(args)) {
      params$out_csv <- args[i + 1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

write_empty <- function(path, cols) {
  if (is.null(path) || path == "") return()
  empty_df <- data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(empty_df) <- cols
  write.table(empty_df, path, sep = ";", row.names = FALSE, quote = FALSE)
}

# Détection du type de fichier en lisant l'en-tête
is_indel_file <- grepl("indels", params$out_csv, ignore.case = TRUE) || grepl("indels", params$input_csv, ignore.case = TRUE)

COLS_OUT_VARIANTS <- c(
  "SAMPLE", "GENE", "GENOTYPE", "NUM_CODON", "POS_TYPE",
  "REF_AA", "ALT_AA", "Mutation_name", "freq", "base_status",
  "count", "total_count", "REF_POS"
)

COLS_OUT_INDELS <- c(
  "SAMPLE", "GENE", "GENOTYPE", "REF_POS", "CODON",
  "INDEL_TYPE", "INDEL_SIZE_BP", "CONSEQUENCE",
  "REF_ALLELE", "ALT_ALLELE", "INDEL_NAME",
  "FREQ", "STATUS", "ALT_READS", "TOTAL_DEPTH"
)

expected_cols <- if (is_indel_file) COLS_OUT_INDELS else COLS_OUT_VARIANTS

if (is.null(params$input_csv) || !file.exists(params$input_csv) || file.info(params$input_csv)$size == 0) {
  write_empty(params$out_csv, expected_cols)
  quit(save = "no", status = 0)
}

df <- tryCatch(
  read.csv2(params$input_csv, stringsAsFactors = FALSE),
  error = function(e) data.frame()
)

if (nrow(df) == 0) {
  write_empty(params$out_csv, expected_cols)
  quit(save = "no", status = 0)
}

# 1. Filtre alertes
df_mut <- df[grepl("^alerte", df$warning), ]
if (nrow(df_mut) == 0) {
  write_empty(params$out_csv, expected_cols)
  quit(save = "no", status = 0)
}

# 2. Règle clinique : majo OU >= seuil (8 %)
if (is_indel_file) {
  df_mut$FREQ <- as.numeric(df_mut$FREQ)
  df_selected <- df_mut[df_mut$STATUS == "majo" | df_mut$FREQ >= params$min_freq_minor, ]
} else {
  df_mut <- df_mut[df_mut$REF_AA != df_mut$ALT_AA, ]
  df_mut$freq <- as.numeric(df_mut$freq)
  df_selected <- df_mut[df_mut$base_status == "majo" | df_mut$freq >= params$min_freq_minor, ]
}

if (nrow(df_selected) == 0) {
  write_empty(params$out_csv, expected_cols)
  quit(save = "no", status = 0)
}

# 3. Formatage de sortie
df_selected$SAMPLE <- params$sample_id

for (col in expected_cols) {
  if (!col %in% colnames(df_selected)) {
    df_selected[[col]] <- NA
  }
}

df_final <- df_selected[, expected_cols, drop = FALSE]

# 4. Dédoublonnage et tri
if (is_indel_file) {
  df_final <- df_final[!duplicated(df_final[, c("SAMPLE", "GENE", "REF_POS", "INDEL_NAME", "STATUS")]), ]
  df_final$REF_POS_NUM <- as.numeric(df_final$REF_POS)
  df_final <- df_final[order(df_final$GENE, df_final$REF_POS_NUM), ]
  df_final$REF_POS_NUM <- NULL
} else {
  df_final <- df_final[!duplicated(df_final[, c("SAMPLE", "GENE", "REF_POS", "Mutation_name", "base_status")]), ]
  df_final$REF_POS_NUM <- as.numeric(df_final$REF_POS)
  df_final <- df_final[order(df_final$GENE, df_final$REF_POS_NUM), ]
  df_final$REF_POS_NUM <- NULL
}

write.table(df_final, params$out_csv, sep = ";", row.names = FALSE, quote = FALSE)