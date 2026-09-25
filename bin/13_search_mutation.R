#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    vcf = NULL,
    tables_dir = NULL,
    ref_fasta = NULL,
    genotype = "D",
    freq_min = 8.0,
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
    } else if (arg == "--ref-fasta" && i < length(args)) {
      params$ref_fasta <- args[i + 1]; i <- i + 2
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

cli_params       <- parse_args(args)
vcf_input        <- cli_params$vcf
tablemut_input   <- cli_params$tables_dir
ref_fasta_input  <- cli_params$ref_fasta
genotype         <- cli_params$genotype
filter_threshold <- cli_params$freq_min

if (!is.na(filter_threshold) && filter_threshold < 1.0 && filter_threshold > 0) {
  filter_threshold <- filter_threshold * 100.0
}

window_pos <- cli_params$window_pos
sample_id  <- cli_params$sample_id
file_out_f <- cli_params$out_variants

dir.create("filtered", showWarnings = FALSE, recursive = TRUE)

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
  codon <- toupper(codon)
  if (is.na(codon) || grepl("N", codon) || nchar(codon) != 3) return("ININT")
  if (codon %in% names(GENETIC_CODE)) return(unname(GENETIC_CODE[codon])) else return("ININT")
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

# 1. Sélection de la table de mutations
gt_letter <- toupper(substr(gsub("[^A-Za-z]", "", genotype), 1, 1))
if (gt_letter == "") gt_letter <- "D"
cat(sprintf("[13_search_mutation] sample=%s genotype_brut='%s' -> gt_letter='%s'\n", sample_id, genotype, gt_letter))

table_name <- paste0("mutation_GT", gt_letter, ".csv")
mut_table_path <- file.path(tablemut_input, table_name)

if (gt_letter %in% c("G", "H") || !file.exists(mut_table_path)) {
  cat(sprintf("[13_search_mutation] Pas de table de mutation pour GT%s, sortie vide.\n", gt_letter))
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

table_mut <- read.csv2(mut_table_path, stringsAsFactors = FALSE)
for (col in colnames(table_mut)) {
  if (is.character(table_mut[[col]])) {
    table_mut[[col]] <- trimws(table_mut[[col]])
  }
}

# 2. Extraction du contig dans le FASTA de référence
ref_seq <- ""
if (!is.null(ref_fasta_input) && file.exists(ref_fasta_input)) {
  f_lines <- readLines(ref_fasta_input)
  target_header <- paste0("GT", gt_letter)

  in_seq <- FALSE
  seq_chunks <- c()
  for (line in f_lines) {
    if (startsWith(line, ">")) {
      header_id <- gsub("^>\\s*", "", line)
      header_id <- strsplit(header_id, "\\s+")[[1]][1]
      if (header_id == target_header) {
        in_seq <- TRUE
      } else {
        if (in_seq) break
        in_seq <- FALSE
      }
    } else if (in_seq) {
      seq_chunks <- c(seq_chunks, trimws(line))
    }
  }
  ref_seq <- toupper(paste(seq_chunks, collapse = ""))
}

if (nchar(ref_seq) == 0) {
  cat(sprintf("⚠️ [13_search_mutation] ALERTE : contig '%s' introuvable dans %s pour %s.\n",
              paste0("GT", gt_letter), ref_fasta_input, sample_id))
} else {
  cat(sprintf("[13_search_mutation] Contig GT%s trouvé, longueur=%d bp.\n", gt_letter, nchar(ref_seq)))
}

# ==============================================================================
# 3. Parsing du VCF (Extraction ordonnée des allèles réels observés)
# ==============================================================================
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
if (ncol(raw_vcf) >= 7) colnames(raw_vcf)[7] <- "FILTER"
if (ncol(raw_vcf) >= 9) colnames(raw_vcf)[9] <- "FORMAT"
if (ncol(raw_vcf) >= 10) colnames(raw_vcf)[10] <- "SAMPLE"

parsed_records <- list()
rank_tags <- c("majo", "mino_1st", "mino_2nd", "mino_3d")
all_possible_bases <- c("A", "C", "G", "T")

for (i in seq_len(nrow(raw_vcf))) {
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
      if (idx <= length(fmt_vals)) {
        ad_vals <- as.numeric(unlist(strsplit(fmt_vals[idx], ",")))
      }
    }
  }

  ref_base <- toupper(substr(raw_vcf$REF[i], 1, 1))
  raw_alt  <- raw_vcf$ALT[i]
  
  base_counts <- setNames(rep(0, 4), all_possible_bases)
  
  if (ref_base %in% all_possible_bases) {
    if (length(ad_vals) >= 1) {
      base_counts[ref_base] <- ad_vals[1]
    } else {
      base_counts[ref_base] <- dp_val
    }
  }
  
  if (!is.na(raw_alt) && raw_alt != "." && raw_alt != "") {
    alt_list <- unlist(strsplit(raw_alt, ","))
    for (k in seq_along(alt_list)) {
      b_alt <- toupper(substr(alt_list[k], 1, 1))
      if (b_alt %in% all_possible_bases) {
        cnt <- 0
        if (length(ad_vals) >= (k + 1)) {
          cnt <- ad_vals[k + 1]
        }
        base_counts[b_alt] <- cnt
      }
    }
  }
  
  sum_counts <- sum(base_counts)
  actual_dp <- max(dp_val, sum_counts)
  if (actual_dp == 0) next
  
  sorted_indices <- order(base_counts, decreasing = TRUE)
  sorted_bases   <- names(base_counts)[sorted_indices]
  sorted_counts  <- base_counts[sorted_indices]
  
  for (r in 1:4) {
    b_name <- sorted_bases[r]
    b_cnt  <- sorted_counts[r]
    b_stat <- rank_tags[r]
    
    if (b_cnt == 0 && b_stat != "majo") {
      next
    }
    
    b_freq <- round((b_cnt / actual_dp) * 100.0, 2)
    
    parsed_records[[length(parsed_records) + 1]] <- data.frame(
      REFERENCE   = raw_vcf$CHROM[i],
      REF_POS     = as.numeric(raw_vcf$POS[i]),
      REF         = ref_base,
      base        = b_name,
      base_status = b_stat,
      count       = b_cnt,
      total_count = actual_dp,
      freq        = b_freq,
      stringsAsFactors = FALSE
    )
  }
}

if (length(parsed_records) == 0) {
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

table_vcf <- do.call(rbind, parsed_records)
cat(sprintf("[13_search_mutation] %d allèles retenus sur %d positions.\n",
            nrow(table_vcf), length(unique(table_vcf$REF_POS))))

# ==============================================================================
# 4. Criblage par Codon et par Nucléotide (Avec Composition Multi-Nucléotidique)
# ==============================================================================
searchMUT_CODON <- function(vcf, mut_table, pattern_region, gene_pfx) {
  mut_sub <- mut_table[grepl(pattern_region, mut_table$Region, ignore.case = TRUE), ]
  if (nrow(mut_sub) == 0 || nrow(vcf) == 0) return(data.frame())
  results <- list()

  for (i in seq_len(nrow(mut_sub))) {
    pos_eco <- mut_sub$Position_EcoR1[i]
    p1 <- as.numeric(mut_sub$NUC1_P3[i]) + window_pos
    p2 <- as.numeric(mut_sub$NUC2_P3[i]) + window_pos
    p3 <- as.numeric(mut_sub$NUC3_P3[i]) + window_pos

    codon_ref <- ""
    if (nchar(ref_seq) >= max(p1, p2, p3, na.rm = TRUE)) {
      b1 <- substr(ref_seq, p1, p1)
      b2 <- substr(ref_seq, p2, p2)
      b3 <- substr(ref_seq, p3, p3)
      codon_ref <- paste0(b1, b2, b3)
    }

    hit_vcf <- vcf[vcf$REF_POS %in% c(p1, p2, p3), ]
    if (nrow(hit_vcf) > 0) {
      ref_aa <- if (nchar(codon_ref) == 3) getAA(codon_ref) else "X"

      # Reconstitution du codon majoritaire complet en combinant toutes les positions majoritaires observées
      codon_majo_composite <- codon_ref
      if (nchar(codon_ref) == 3) {
        majo_hits <- hit_vcf[hit_vcf$base_status == "majo", ]
        if (nrow(majo_hits) > 0) {
          for (m in seq_len(nrow(majo_hits))) {
            pos_m <- majo_hits$REF_POS[m]
            if (pos_m == p1) substr(codon_majo_composite, 1, 1) <- majo_hits$base[m]
            if (pos_m == p2) substr(codon_majo_composite, 2, 2) <- majo_hits$base[m]
            if (pos_m == p3) substr(codon_majo_composite, 3, 3) <- majo_hits$base[m]
          }
        }
      }

      for (k in seq_len(nrow(hit_vcf))) {
        r_row <- hit_vcf[k, ]
        pt <- ifelse(r_row$REF_POS == p1, "NUC1_P3", ifelse(r_row$REF_POS == p2, "NUC2_P3", "NUC3_P3"))

        if (nchar(codon_ref) == 3) {
          if (r_row$base_status == "majo") {
            alt_aa <- getAA(codon_majo_composite)
          } else {
            codon_mino <- codon_majo_composite
            if (pt == "NUC1_P3") substr(codon_mino, 1, 1) <- r_row$base
            if (pt == "NUC2_P3") substr(codon_mino, 2, 2) <- r_row$base
            if (pt == "NUC3_P3") substr(codon_mino, 3, 3) <- r_row$base
            alt_aa <- getAA(codon_mino)
          }
        } else {
          ref_aa <- r_row$REF
          alt_aa <- r_row$base
        }

        mut_name <- paste0(gene_pfx, ref_aa, pos_eco, alt_aa)
        
        # Alerte si l'acide aminé est modifié ET que ce nucléotide précis a muté
        nt_has_mutated <- (r_row$base != r_row$REF)
        warn <- ifelse(ref_aa != alt_aa && nt_has_mutated, "alerte!", "OK")

        results[[length(results) + 1]] <- data.frame(
          REFERENCE     = r_row$REFERENCE,
          REF_POS       = r_row$REF_POS,
          REF           = r_row$REF,
          total_count   = r_row$total_count,
          base_status   = r_row$base_status,
          base          = r_row$base,
          count         = r_row$count,
          NUM_CODON     = pos_eco,
          POS_TYPE      = pt,
          REF_AA        = ref_aa,
          ALT_AA        = alt_aa,
          Mutation_name = mut_name,
          freq          = r_row$freq,
          warning       = warn,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) return(data.frame())
  return(do.call(rbind, results))
}

searchMUT_NT <- function(vcf, mut_table, pattern_region) {
  mut_sub <- mut_table[grepl(pattern_region, mut_table$Region, ignore.case = TRUE), ]
  if (nrow(mut_sub) == 0 || nrow(vcf) == 0) return(data.frame())

  mut_sub$REF_POS <- as.numeric(mut_sub$NUC1_P3) + window_pos
  table_vcf_mut <- merge(vcf, mut_sub, by = "REF_POS")
  if (nrow(table_vcf_mut) == 0) return(data.frame())

  table_vcf_mut$Mutation_name <- paste0(table_vcf_mut$REF, as.character(table_vcf_mut$Position_EcoR1), table_vcf_mut$base)
  table_vcf_mut$warning       <- ifelse(table_vcf_mut$REF == table_vcf_mut$base, "OK", "alerte!")
  table_vcf_mut$NUM_CODON     <- table_vcf_mut$Position_EcoR1
  table_vcf_mut$POS_TYPE      <- "NUC1_P3"
  table_vcf_mut$REF_AA        <- table_vcf_mut$REF
  table_vcf_mut$ALT_AA        <- table_vcf_mut$base

  cols_order <- c("REFERENCE", "REF_POS", "REF", "total_count", "base_status",
                  "base", "count", "NUM_CODON", "POS_TYPE", "REF_AA", "ALT_AA", "Mutation_name", "freq", "warning")
  cols_avail <- intersect(cols_order, colnames(table_vcf_mut))
  return(table_vcf_mut[, cols_avail, drop = FALSE])
}

add_meta <- function(df, gene) {
  if (nrow(df) > 0) {
    df$GENE <- gene
    df$GENOTYPE <- paste0("GT", gt_letter)
  }
  return(df)
}

regions_list <- list(
  add_meta(searchMUT_NT(table_vcf, table_mut, "BCP"), "BCP"),
  add_meta(searchMUT_NT(table_vcf, table_mut, "PreCore"), "PreCore"),
  add_meta(searchMUT_CODON(table_vcf, table_mut, "RT", "rt"), "Domaine RT"),
  add_meta(searchMUT_CODON(table_vcf, table_mut, "^Core", "c"), "Core"),
  add_meta(searchMUT_CODON(table_vcf, table_mut, "Domaine S", "s"), "Domaine S"),
  add_meta(searchMUT_CODON(table_vcf, table_mut, "PreS1", "preS1_"), "Domaine PreS1"),
  add_meta(searchMUT_CODON(table_vcf, table_mut, "PreS2", "preS2_"), "Domaine PreS2"),
  add_meta(searchMUT_CODON(table_vcf, table_mut, "HBx", "x"), "Domaine HBx")
)

valid_dfs <- Filter(function(x) nrow(x) > 0, regions_list)

if (length(valid_dfs) == 0) {
  cat("[13_search_mutation] Aucune mutation trouvée dans les régions criblées.\n")
  write_empty(file_out_f)
  quit(save = "no", status = 0)
}

all_cols <- unique(unlist(lapply(valid_dfs, colnames)))
aligned <- lapply(valid_dfs, function(df) {
  for (col in setdiff(all_cols, colnames(df))) df[[col]] <- NA
  return(df[, all_cols, drop = FALSE])
})

df_combined <- do.call(rbind, aligned)
df_filtered <- subset(df_combined, freq >= filter_threshold)

cat(sprintf("[13_search_mutation] %d mutations avant filtre freq>=%.1f%%, %d après filtre.\n",
            nrow(df_combined), filter_threshold, nrow(df_filtered)))

if (nrow(df_filtered) > 0) {
  write.table(df_filtered, file_out_f, sep = ";", row.names = FALSE, quote = FALSE)
} else {
  write_empty(file_out_f)
}