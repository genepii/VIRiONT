#!/usr/bin/env Rscript

# =====================================================================
# VIRiONT_NF - Screening des Mutations VHB (Sortie Structurée & Unique)
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    vcf = NULL,
    tables_dir = NULL,
    genotype = "GTD",
    freq_min = 0.08,
    window_pos = 0,
    sample_id = "sample",
    out_variants = ""
  )
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--vcf" && i < length(args)) {
      params$vcf <- args[i + 1]; i <- i + 2
    } else if (arg == "--tables-dir" && i < length(args)) {
      params$tables_dir <- args[i + 1]; i <- i + 2
    } else if (arg == "--genotype" && i < length(args)) {
      params$genotype <- args[i + 1]; i <- i + 2
    } else if (arg == "--freq-min" && i < length(args)) {
      params$freq_min <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--window-pos" && i < length(args)) {
      params$window_pos <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--sample-id" && i < length(args)) {
      params$sample_id <- args[i + 1]; i <- i + 2
    } else if (arg == "--out-variants" && i < length(args)) {
      params$out_variants <- args[i + 1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

if (length(args) >= 2 && grepl("^--", args[1])) {
  cli_params <- parse_args(args)
  vcf_input            <- cli_params$vcf
  tablemut_input       <- cli_params$tables_dir
  genotype             <- cli_params$genotype
  filter_threshold     <- cli_params$freq_min
  window_pos           <- cli_params$window_pos
  sample_id            <- cli_params$sample_id
  file_output_VARIANTS <- cli_params$out_variants
} else {
  stop("Usage: Rscript 16_search_mutation.R --vcf <vcf> --tables-dir <dir> --genotype <GTD> --sample-id <id>")
}

# --- CRÉATION DES DOSSIERS ARBORESCENTS ---
dir.create("all_results", showWarnings = FALSE, recursive = TRUE)
dir.create("filtered", showWarnings = FALSE, recursive = TRUE)

clean_gt <- ifelse(startsWith(genotype, "GT"), genotype, paste0("GT", genotype))
prefix   <- paste0(sample_id, "_", clean_gt)

file_output_ALL_CSV      <- file.path("all_results", paste0(prefix, "_all_mutations.csv"))
file_output_FILTERED_CSV <- file.path("filtered", paste0(prefix, "_filtered_mutations.csv"))

if (is.null(file_output_VARIANTS) || file_output_VARIANTS == "") {
  file_output_VARIANTS <- file_output_FILTERED_CSV
}

# --- CODE GÉNÉTIQUE ---
GENETIC_CODE <- c(
  "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L", "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
  "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*", "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
  "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L", "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
  "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q", "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
  "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M", "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
  "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K", "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
  "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V", "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
  "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E", "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
)

getAA <- function(codon) {
  if (is.na(codon) || grepl("N", codon) || nchar(codon) != 3) return("ININT")
  aa <- GENETIC_CODE[codon]
  if (is.na(aa)) return("ININT") else return(unname(aa))
}

searchMUT_CODON <- function(vcf, mut_table, mutation) {
  mut_sub <- subset(mut_table, Region == mutation)
  empty_df <- data.frame(REFERENCE=character(), REF_POS=numeric(), REF=character(), total_count=numeric(), 
                          base_status=character(), base=character(), count=numeric(), NUM_CODON=numeric(), 
                          POS_TYPE=character(), REF_AA=character(), ALT_AA=character(), Mutation_name=character(), 
                          freq=numeric(), warning=character(), stringsAsFactors=FALSE)
  
  if (nrow(mut_sub) == 0 || nrow(vcf) == 0) return(empty_df)
  
  records <- list()
  for (i in seq_len(nrow(mut_sub))) {
    pos_eco <- mut_sub$Position_EcoR1[i]
    for (pt in c("NUC1_P3", "NUC2_P3", "NUC3_P3")) {
      val <- mut_sub[[pt]][i]
      if (!is.na(val) && val != "") {
        records[[length(records) + 1]] <- data.frame(
          Position_EcoR1 = pos_eco,
          POS_TYPE = pt,
          REF_POS = as.numeric(val) + window_pos,
          NUM_CODON = pos_eco,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(records) == 0) return(empty_df)
  mut_long <- do.call(rbind, records)
  
  table_vcf_mut <- merge(vcf, mut_long, by = "REF_POS")
  
  if (nrow(table_vcf_mut) == 0 && min(vcf$REF_POS) > 0) {
    offset <- min(mut_long$REF_POS) - min(vcf$REF_POS)
    vcf_shift <- vcf
    vcf_shift$REF_POS <- vcf_shift$REF_POS + offset
    table_vcf_mut <- merge(vcf_shift, mut_long, by = "REF_POS")
  }
  
  if (nrow(table_vcf_mut) == 0) return(empty_df)
  
  table_vcf_mut$CODON_REF <- "ATG"
  table_vcf_mut$CODON_ALT <- table_vcf_mut$CODON_REF
  
  for (i in seq_len(nrow(table_vcf_mut))) {
    b <- substr(table_vcf_mut$base[i], 1, 1)
    pt <- table_vcf_mut$POS_TYPE[i]
    cr <- "ATG"
    if (pt == "NUC1_P3") substr(cr, 1, 1) <- b
    if (pt == "NUC2_P3") substr(cr, 2, 2) <- b
    if (pt == "NUC3_P3") substr(cr, 3, 3) <- b
    table_vcf_mut$CODON_ALT[i] <- cr
  }
  
  table_vcf_mut$REF_AA <- sapply(table_vcf_mut$CODON_REF, getAA)
  table_vcf_mut$ALT_AA <- sapply(table_vcf_mut$CODON_ALT, getAA)
  table_vcf_mut$warning <- ifelse(table_vcf_mut$REF_AA == table_vcf_mut$ALT_AA, "OK", "alerte!")
  table_vcf_mut$Mutation_name <- paste0(table_vcf_mut$REF_AA, as.character(table_vcf_mut$NUM_CODON), table_vcf_mut$ALT_AA)
  
  cols_order <- c("REFERENCE", "REF_POS", "REF", "total_count", "base_status", "base", "count", "NUM_CODON", "POS_TYPE", "REF_AA", "ALT_AA", "Mutation_name", "freq", "warning")
  cols_avail <- intersect(cols_order, colnames(table_vcf_mut))
  return(table_vcf_mut[order(table_vcf_mut$REF_POS), cols_avail])
}

searchMUT_NT <- function(vcf, mut_table, mutation) {
  mut_sub <- subset(mut_table, Region == mutation)
  empty_df <- data.frame(REFERENCE=character(), REF_POS=numeric(), Position_EcoR1=numeric(), REF=character(), 
                          total_count=numeric(), base_status=character(), base=character(), count=numeric(), 
                          Mutation_name=character(), freq=numeric(), stringsAsFactors=FALSE)
  
  if (nrow(mut_sub) == 0 || nrow(vcf) == 0) return(empty_df)
  
  mut_sub$REF_POS <- as.numeric(mut_sub$NUC1_P3) + window_pos
  table_vcf_mut <- merge(vcf, mut_sub, by = "REF_POS")
  
  if (nrow(table_vcf_mut) == 0) {
    offset <- min(mut_sub$REF_POS, na.rm=T) - min(vcf$REF_POS, na.rm=T)
    vcf_shift <- vcf
    vcf_shift$REF_POS <- vcf_shift$REF_POS + offset
    table_vcf_mut <- merge(vcf_shift, mut_sub, by = "REF_POS")
  }
  
  if (nrow(table_vcf_mut) == 0) return(empty_df)
  
  table_vcf_mut$Mutation_name <- paste0(table_vcf_mut$REF, as.character(table_vcf_mut$Position_EcoR1), table_vcf_mut$base)
  cols_order <- c("REFERENCE", "REF_POS", "Position_EcoR1", "REF", "total_count", "base_status", "base", "count", "Mutation_name", "freq")
  cols_avail <- intersect(cols_order, colnames(table_vcf_mut))
  return(table_vcf_mut[order(table_vcf_mut$REF_POS), cols_avail])
}

write_empty <- function(path) {
  if (is.null(path) || is.na(path) || path == "") return()
  empty_df <- data.frame(REFERENCE=character(), GENE=character(), GENOTYPE=character(), REF_POS=numeric(), Position_EcoR1=numeric(), REF=character(), total_count=numeric(), base_status=character(), base=character(), count=numeric(), Mutation_name=character(), freq=numeric(), ALT=character(), NUM_CODON=numeric(), POS_TYPE=character(), REF_AA=character(), ALT_AA=character(), warning=character())
  write.csv2(empty_df, path, row.names=FALSE, quote=FALSE)
}

if (is.null(vcf_input) || !file.exists(vcf_input) || file.info(vcf_input)$size == 0) {
  write_empty(file_output_ALL_CSV)
  write_empty(file_output_FILTERED_CSV)
  quit(save = "no", status = 0)
}

vcf_lines <- readLines(vcf_input)
vcf_clean_lines <- vcf_lines[!startsWith(vcf_lines, "##")]

if (length(vcf_clean_lines) <= 1) {
  write_empty(file_output_ALL_CSV)
  write_empty(file_output_FILTERED_CSV)
  quit(save = "no", status = 0)
}

raw_vcf <- read.delim(text = vcf_clean_lines, sep = "\t", stringsAsFactors = FALSE)
colnames(raw_vcf)[1] <- "CHROM"

if (nrow(raw_vcf) == 0) {
  write_empty(file_output_ALL_CSV)
  write_empty(file_output_FILTERED_CSV)
  quit(save = "no", status = 0)
}

parsed_records <- list()
for (i in seq_len(nrow(raw_vcf))) {
  fmt_keys <- unlist(strsplit(raw_vcf$FORMAT[i], ":"))
  sample_col <- ifelse("SAMPLE" %in% colnames(raw_vcf), raw_vcf$SAMPLE[i], raw_vcf[[ncol(raw_vcf)]][i])
  fmt_vals <- unlist(strsplit(sample_col, ":"))
  
  dp_val <- 100; af_val <- 100.0
  if ("DP" %in% fmt_keys) dp_val <- as.numeric(fmt_vals[which(fmt_keys == "DP")])
  if ("AF" %in% fmt_keys) af_val <- round(as.numeric(fmt_vals[which(fmt_keys == "AF")]) * 100, 2)
  
  parsed_records[[i]] <- data.frame(
    REFERENCE   = raw_vcf$CHROM[i],
    REF_POS     = as.numeric(raw_vcf$POS[i]),
    REF         = raw_vcf$REF[i],
    ALT         = raw_vcf$ALT[i],
    base        = raw_vcf$ALT[i],
    base_status = "majo",
    count       = round(dp_val * (af_val / 100)),
    total_count = dp_val,
    freq        = af_val,
    stringsAsFactors = FALSE
  )
}
table_vcf <- do.call(rbind, parsed_records)

# --- BASE DE DONNÉES GÉNOTYPE ---
table_name <- paste0("mutation_", clean_gt, ".csv")
mut_table_path <- file.path(tablemut_input, table_name)

if (!file.exists(mut_table_path)) {
  mut_table_path <- file.path(tablemut_input, "mutation_GTD.csv")
}

if (!file.exists(mut_table_path)) {
  write_empty(file_output_ALL_CSV)
  write_empty(file_output_FILTERED_CSV)
  quit(save = "no", status = 0)
}

table_mut <- read.csv2(mut_table_path, stringsAsFactors = FALSE)
table_mut$Region <- as.character(table_mut$Region)

add_info <- function(df, gene_name, gt_val) {
  if (nrow(df) > 0) {
    df$GENE     <- gene_name
    df$GENOTYPE <- gt_val
    if (!"ALT" %in% colnames(df) && "base" %in% colnames(df)) df$ALT <- df$base
  }
  return(df)
}

# Extraction par région
raw_list <- list(
  add_info(searchMUT_NT(table_vcf, table_mut, "BCP"), "BCP", clean_gt),
  add_info(searchMUT_NT(table_vcf, table_mut, "PreCore"), "PreCore", clean_gt),
  add_info(searchMUT_CODON(table_vcf, table_mut, "Domaine RT"), "Domaine RT", clean_gt),
  add_info(searchMUT_CODON(table_vcf, table_mut, "Core"), "Core", clean_gt),
  add_info(searchMUT_CODON(table_vcf, table_mut, "Domaine S"), "Domaine S", clean_gt),
  add_info(searchMUT_CODON(table_vcf, table_mut, "Domaine PreS1"), "Domaine PreS1", clean_gt),
  add_info(searchMUT_CODON(table_vcf, table_mut, "Domaine PreS2"), "Domaine PreS2", clean_gt),
  add_info(searchMUT_CODON(table_vcf, table_mut, "Domaine HBx"), "Domaine HBx", clean_gt)
)

combine_and_order <- function(df_list) {
  valid_dfs <- Filter(function(x) nrow(x) > 0, df_list)
  if (length(valid_dfs) > 0) {
    all_cols <- unique(unlist(lapply(valid_dfs, colnames)))
    front_cols <- c("REFERENCE", "GENE", "GENOTYPE")
    all_cols <- c(front_cols, setdiff(all_cols, front_cols))
    
    aligned <- lapply(valid_dfs, function(df) {
      missing_cols <- setdiff(all_cols, colnames(df))
      for (col in missing_cols) df[[col]] <- NA
      return(df[, all_cols, drop = FALSE])
    })
    return(do.call(rbind, aligned))
  }
  return(NULL)
}

# --- TABLE 1: BRUT (ALL RESULTS) ---
df_all <- combine_and_order(raw_list)
if (!is.null(df_all)) {
  write.csv2(df_all, file_output_ALL_CSV, row.names = FALSE, quote = FALSE)
} else {
  write_empty(file_output_ALL_CSV)
}

# --- TABLE 2: FILTRÉ (MAJORITAIRE >= freq_min) ---
filter_tbl <- function(df) {
  if (nrow(df) > 0 && "freq" %in% colnames(df)) {
    return(subset(df, as.numeric(gsub(",", ".", as.character(freq))) >= filter_threshold))
  }
  return(df)
}

filtered_list <- lapply(raw_list, filter_tbl)
df_filtered   <- combine_and_order(filtered_list)

if (!is.null(df_filtered)) {
  write.csv2(df_filtered, file_output_FILTERED_CSV, row.names = FALSE, quote = FALSE)
} else {
  write_empty(file_output_FILTERED_CSV)
}

cat("Screening Rscript terminé avec succès.\n")