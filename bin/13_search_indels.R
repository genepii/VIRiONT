#!/usr/bin/env Rscript

# ==============================================================================
# 13_search_indels.R - Détection et Annotation Dédiée des INDELs & Frameshifts VHB
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    vcf = NULL,
    tables_dir = NULL,
    ref_fasta = NULL,
    genotype = "D",
    freq_min = 0.0,
    min_dp = 20,
    min_alt = 5,
    window_pos = 0,
    sample_id = "sample",
    out_indels = ""
  )
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--vcf" && i < length(args)) {
      params$vcf <- args[i + 1]; i <- i + 2
    } else if (arg == "--tables-dir" && i < length(args)) {
      params$tables_dir <- args[i + 1]; i <- i + 2
    } else if (arg == "--ref-fasta" && i < length(args)) {
      params$ref_fasta <- args[i + 1]; i <- i + 2
    } else if (arg == "--genotype" && i < length(args)) {
      params$genotype <- args[i + 1]; i <- i + 2
    } else if (arg == "--freq-min" && i < length(args)) {
      params$freq_min <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--min-dp" && i < length(args)) {
      params$min_dp <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--min-alt" && i < length(args)) {
      params$min_alt <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--window-pos" && i < length(args)) {
      params$window_pos <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--sample-id" && i < length(args)) {
      params$sample_id <- args[i + 1]; i <- i + 2
    } else if (arg == "--out-indels" && i < length(args)) {
      params$out_indels <- args[i + 1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

cli_params       <- parse_args(args)
vcf_input        <- cli_params$vcf
tablemut_input   <- cli_params$tables_dir
genotype         <- cli_params$genotype
filter_threshold <- cli_params$freq_min
min_dp           <- cli_params$min_dp
min_alt          <- cli_params$min_alt
window_pos       <- cli_params$window_pos
sample_id        <- cli_params$sample_id
file_out_f       <- cli_params$out_indels

if (!is.na(filter_threshold) && filter_threshold < 1.0 && filter_threshold > 0) {
  filter_threshold <- filter_threshold * 100.0
}

write_empty <- function(path) {
  if (is.null(path) || path == "") return()
  empty_df <- data.frame(
    REFERENCE=character(), REF_POS=numeric(), REF=character(),
    total_count=numeric(), base_status=character(), base=character(),
    count=numeric(), NUM_CODON=numeric(), POS_TYPE=character(),
    REF_AA=character(), ALT_AA=character(), Mutation_name=character(),
    freq=numeric(), warning=character(), GENE=character(),
    GENOTYPE=character(), stringsAsFactors=FALSE
  )
  write.table(empty_df, path, sep=";", row.names=FALSE, quote=FALSE)
}

# 1. Chargement de la table de mutations du génotype
gt_letter <- toupper(substr(gsub("[^A-Za-z]", "", genotype), 1, 1))
if (gt_letter == "") gt_letter <- "D"

table_name <- paste0("mutation_GT", gt_letter, ".csv")
mut_table_path <- file.path(tablemut_input, table_name)

if (gt_letter %in% c("G", "H") || !file.exists(mut_table_path)) {
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

table_mut <- read.csv2(mut_table_path, stringsAsFactors = FALSE)
for (col in colnames(table_mut)) {
  if (is.character(table_mut[[col]])) table_mut[[col]] <- trimws(table_mut[[col]])
}

# 2. Parsing du VCF : Extraction EXCLUSIVE des INDELs avec filtrage qualité
if (is.null(vcf_input) || !file.exists(vcf_input) || file.info(vcf_input)$size == 0) {
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

vcf_lines <- readLines(vcf_input)
vcf_clean <- vcf_lines[!grepl("^#", vcf_lines)]

if (length(vcf_clean) == 0) {
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

raw_vcf <- read.delim(text = vcf_clean, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(raw_vcf)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
if (ncol(raw_vcf) >= 9) colnames(raw_vcf)[9] <- "FORMAT"
if (ncol(raw_vcf) >= 10) colnames(raw_vcf)[10] <- "SAMPLE"

parsed_indels <- list()

for (i in seq_len(nrow(raw_vcf))) {
  ref_str <- toupper(raw_vcf$REF[i])
  raw_alt <- raw_vcf$ALT[i]
  if (is.na(raw_alt) || raw_alt == "." || raw_alt == "") next

  alt_list <- unlist(strsplit(raw_alt, ","))
  
  dp_val <- 0
  ad_vals <- c()
  if ("FORMAT" %in% colnames(raw_vcf) && "SAMPLE" %in% colnames(raw_vcf)) {
    fmt_keys <- unlist(strsplit(raw_vcf$FORMAT[i], ":"))
    fmt_vals <- unlist(strsplit(raw_vcf$SAMPLE[i], ":"))
    if ("DP" %in% fmt_keys) {
      idx <- which(fmt_keys == "DP")
      if (idx <= length(fmt_vals)) dp_val <- as.numeric(fmt_vals[idx])
    }
    if ("AD" %in% fmt_keys) {
      idx <- which(fmt_keys == "AD")
      if (idx <= length(fmt_vals)) ad_vals <- as.numeric(unlist(strsplit(fmt_vals[idx], ",")))
    }
  }
  actual_dp <- max(dp_val, sum(ad_vals))

  # VERROU 1 : Élimination des sites sous-couverts (bords d'amplicons)
  if (actual_dp < min_dp) next

  for (k in seq_along(alt_list)) {
    alt_str <- toupper(alt_list[k])
    
    # CONDITION FORMELLE : Ne retenir QUE les vrais indels (variation de taille)
    if (nchar(ref_str) != nchar(alt_str)) {
      cnt <- if (length(ad_vals) >= (k + 1)) ad_vals[k + 1] else 0
      
      # VERROU 2 : Filtre anti-bruit (support minimal en reads mutés)
      if (cnt >= min_alt) {
        freq_ind <- round((cnt / actual_dp) * 100.0, 2)
        stat <- if (freq_ind >= 50.0) "majo" else "mino_indel"

        parsed_indels[[length(parsed_indels) + 1]] <- data.frame(
          REFERENCE   = raw_vcf$CHROM[i],
          REF_POS     = as.numeric(raw_vcf$POS[i]),
          REF         = ref_str,
          ALT         = alt_str,
          count       = cnt,
          total_count = actual_dp,
          freq        = freq_ind,
          base_status = stat,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(parsed_indels) == 0) {
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

indels_df <- do.call(rbind, parsed_indels)

# 3. Annotation Clinique des Indels & Frameshifts
region_configs <- list(
  list(pattern = "RT", pfx = "rt", gene = "Domaine RT"),
  list(pattern = "^Core", pfx = "c", gene = "Core"),
  list(pattern = "Domaine S", pfx = "s", gene = "Domaine S"),
  list(pattern = "PreS1", pfx = "preS1_", gene = "Domaine PreS1"),
  list(pattern = "PreS2", pfx = "preS2_", gene = "Domaine PreS2"),
  list(pattern = "HBx", pfx = "x", gene = "Domaine HBx"),
  list(pattern = "PreCore", pfx = "pc_", gene = "PreCore"),
  list(pattern = "BCP", pfx = "bcp_", gene = "BCP")
)

results <- list()

for (r in seq_len(nrow(indels_df))) {
  ind      <- indels_df[r, ]
  pos      <- ind$REF_POS
  ref_len  <- nchar(ind$REF)
  alt_len  <- nchar(ind$ALT)
  len_diff <- alt_len - ref_len
  is_fs    <- (len_diff %% 3 != 0)
  span_pos <- seq(pos, pos + max(0, ref_len - 1))

  matched_any <- FALSE

  for (cfg in region_configs) {
    sub_mut <- table_mut[grepl(cfg$pattern, table_mut$Region, ignore.case = TRUE), ]
    if (nrow(sub_mut) == 0) next

    n1 <- as.numeric(sub_mut$NUC1_P3) + window_pos
    n2 <- as.numeric(sub_mut$NUC2_P3) + window_pos
    n3 <- as.numeric(sub_mut$NUC3_P3) + window_pos

    hits <- which((n1 %in% span_pos) | (n2 %in% span_pos) | (n3 %in% span_pos))

    if (length(hits) > 0) {
      matched_any <- TRUE
      hit_row <- sub_mut[hits[1], ]
      pos_eco <- hit_row$Position_EcoR1

      if (is_fs) {
        pos_type <- "FRAMESHIFT"
        ref_aa   <- "INDEL"
        alt_aa   <- "FS"
        type_tag <- ifelse(len_diff > 0, paste0("ins", len_diff, "bp"), paste0("del", abs(len_diff), "bp"))
        mut_name <- paste0(cfg$pfx, "fs", pos_eco, "_", type_tag)
        warn     <- "alerte! FRAMESHIFT"
      } else {
        num_aa   <- abs(len_diff) / 3
        pos_type <- ifelse(len_diff > 0, "IN_FRAME_INS", "IN_FRAME_DEL")
        ref_aa   <- "INDEL"
        alt_aa   <- ifelse(len_diff > 0, paste0("+", num_aa, "AA"), paste0("-", num_aa, "AA"))
        type_tag <- ifelse(len_diff > 0, paste0("ins", num_aa, "AA"), paste0("del", num_aa, "AA"))
        mut_name <- paste0(cfg$pfx, type_tag, "_pos", pos_eco)
        warn     <- "alerte! INDEL"
      }

      results[[length(results) + 1]] <- data.frame(
        REFERENCE     = ind$REFERENCE,
        REF_POS       = ind$REF_POS,
        REF           = ind$REF,
        total_count   = ind$total_count,
        base_status   = ind$base_status,
        base          = ind$ALT,
        count         = ind$count,
        NUM_CODON     = pos_eco,
        POS_TYPE      = pos_type,
        REF_AA        = ref_aa,
        ALT_AA        = alt_aa,
        Mutation_name = mut_name,
        freq          = ind$freq,
        warning       = warn,
        GENE          = cfg$gene,
        GENOTYPE      = paste0("GT", gt_letter),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!matched_any) {
    type_tag <- ifelse(len_diff > 0, paste0("ins", len_diff, "bp"), paste0("del", abs(len_diff), "bp"))
    results[[length(results) + 1]] <- data.frame(
      REFERENCE     = ind$REFERENCE,
      REF_POS       = ind$REF_POS,
      REF           = ind$REF,
      total_count   = ind$total_count,
      base_status   = ind$base_status,
      base          = ind$ALT,
      count         = ind$count,
      NUM_CODON     = NA,
      POS_TYPE      = ifelse(is_fs, "FRAMESHIFT", ifelse(len_diff > 0, "IN_FRAME_INS", "IN_FRAME_DEL")),
      REF_AA        = "INDEL",
      ALT_AA        = ifelse(is_fs, "FS", ifelse(len_diff > 0, paste0("+", abs(len_diff)/3, "AA"), paste0("-", abs(len_diff)/3, "AA"))),
      Mutation_name = paste0("indel_", pos, "_", type_tag),
      freq          = ind$freq,
      warning       = ifelse(is_fs, "alerte! FRAMESHIFT", "alerte! INDEL"),
      GENE          = "AUTRE_REGION",
      GENOTYPE      = paste0("GT", gt_letter),
      stringsAsFactors = FALSE
    )
  }
}

if (length(results) == 0) {
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

final_df <- do.call(rbind, results)
final_filtered <- subset(final_df, freq >= filter_threshold)

if (nrow(final_filtered) > 0) {
  write.table(final_filtered, file_out_f, sep = ";", row.names = FALSE, quote = FALSE)
  cat(sprintf("[13_search_indels] %d indels/frameshifts exportes vers %s (min_dp=%d, min_alt=%d)\n",
              nrow(final_filtered), file_out_f, min_dp, min_alt))
} else {
  write_empty(file_out_f)
}