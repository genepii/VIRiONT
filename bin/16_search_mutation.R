#!/usr/bin/env Rscript

# =====================================================================
# VIRiONT_NF - Annotation & Screening des Mutations VHB (Base-R Native)
# Compatible Clair3 VCF v4.2 - Ajout de la colonne GENE/Region_Origin
# =====================================================================

argv <- commandArgs(TRUE)

if (length(argv) < 20) {
  stop("Usage: Rscript 16_search_mutation.R <vcf_input> <tablemut_dir> <filter_threshold> <window_pos> <8_raw_outputs> <8_filtered_outputs> [output_variants_global]")
}

vcf_input        <- argv[1]
tablemut_input   <- argv[2]
filter_threshold <- as.numeric(as.character(argv[3]))
window_pos       <- as.numeric(as.character(argv[4]))

file_output_PC    <- argv[5];  file_output_BCP   <- argv[6]
file_output_DS    <- argv[7];  file_output_RT    <- argv[8]
file_output_DPS1  <- argv[9];  file_output_DPS2  <- argv[10]
file_output_DHBx  <- argv[11]; file_output_C     <- argv[12]

file_output_PC_F   <- argv[13]; file_output_BCP_F  <- argv[14]
file_output_DS_F   <- argv[15]; file_output_RT_F   <- argv[16]
file_output_DPS1_F <- argv[17]; file_output_DPS2_F <- argv[18]
file_output_DHBx_F <- argv[19]; file_output_C_F    <- argv[20]

file_output_VARIANTS <- ifelse(length(argv) >= 21, argv[21], "")

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

# --- FONCTION DE RECHERCHE DE MUTATIONS SUR LES CODONS ---
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

# --- FONCTION DE RECHERCHE DE MUTATIONS NUCLÉOTIDIQUES (BCP/PreCore) ---
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

# ---------------------------------------------------------------------
# 1. NETTOYAGE VCF & LECTURE CLAIR3
# ---------------------------------------------------------------------
write_empty <- function(path) {
  if (is.na(path) || path == "") return()
  empty_df <- data.frame(REFERENCE=character(), REF_POS=numeric(), REF=character(), ALT=character(), freq=numeric())
  write.csv2(empty_df, path, row.names=FALSE, quote=FALSE)
}

if (!file.exists(vcf_input) || file.info(vcf_input)$size == 0) {
  lapply(c(file_output_PC, file_output_BCP, file_output_DS, file_output_RT, file_output_DPS1, file_output_DPS2, file_output_DHBx, file_output_C,
           file_output_PC_F, file_output_BCP_F, file_output_DS_F, file_output_RT_F, file_output_DPS1_F, file_output_DPS2_F, file_output_DHBx_F, file_output_C_F,
           file_output_VARIANTS), write_empty)
  quit(save = "no", status = 0)
}

vcf_lines <- readLines(vcf_input)
vcf_clean_lines <- vcf_lines[!startsWith(vcf_lines, "##")]

if (length(vcf_clean_lines) <= 1) {
  lapply(c(file_output_PC, file_output_BCP, file_output_DS, file_output_RT, file_output_DPS1, file_output_DPS2, file_output_DHBx, file_output_C,
           file_output_PC_F, file_output_BCP_F, file_output_DS_F, file_output_RT_F, file_output_DPS1_F, file_output_DPS2_F, file_output_DHBx_F, file_output_C_F,
           file_output_VARIANTS), write_empty)
  quit(save = "no", status = 0)
}

raw_vcf <- read.delim(text = vcf_clean_lines, sep = "\t", stringsAsFactors = FALSE)
colnames(raw_vcf)[1] <- "CHROM"

if (nrow(raw_vcf) == 0) {
  lapply(c(file_output_PC, file_output_BCP, file_output_DS, file_output_RT, file_output_DPS1, file_output_DPS2, file_output_DHBx, file_output_C,
           file_output_PC_F, file_output_BCP_F, file_output_DS_F, file_output_RT_F, file_output_DPS1_F, file_output_DPS2_F, file_output_DHBx_F, file_output_C_F,
           file_output_VARIANTS), write_empty)
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------
# 2. PARSING DE LA COLONNE SAMPLE (EXTRACTION DP & AF)
# ---------------------------------------------------------------------
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

# ---------------------------------------------------------------------
# 3. IDENTIFICATION DE LA TABLE DE GÉNOTYPE (GTA à GTI)
# ---------------------------------------------------------------------
reference <- unique(as.character(table_vcf$REFERENCE))
vcf_filename <- basename(vcf_input)

table_name <- "mutation_GTD.csv"
if (grepl("_A|GTA|A", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTA.csv"
if (grepl("_B|GTB|B", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTB.csv"
if (grepl("_C|GTC|C", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTC.csv"
if (grepl("_D|GTD|D3|D", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTD.csv"
if (grepl("_E|GTE|E", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTE.csv"
if (grepl("_F|GTF|F", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTF.csv"
if (grepl("_G|GTG|G", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTG.csv"
if (grepl("_H|GTH|H", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTH.csv"
if (grepl("_I|GTI|I", paste(reference, vcf_filename), ignore.case=T)) table_name <- "mutation_GTI.csv"

mut_table_path <- file.path(tablemut_input, table_name)

if (!file.exists(mut_table_path)) {
  lapply(c(file_output_PC, file_output_BCP, file_output_DS, file_output_RT, file_output_DPS1, file_output_DPS2, file_output_DHBx, file_output_C,
           file_output_PC_F, file_output_BCP_F, file_output_DS_F, file_output_RT_F, file_output_DPS1_F, file_output_DPS2_F, file_output_DHBx_F, file_output_C_F,
           file_output_VARIANTS), write_empty)
  quit(save = "no", status = 0)
}

table_mut <- read.csv2(mut_table_path, stringsAsFactors = FALSE)
table_mut$Region <- as.character(table_mut$Region)

# ---------------------------------------------------------------------
# 4. EXÉCUTION DU SCREENING
# ---------------------------------------------------------------------
table_vcf_PC   <- searchMUT_NT(table_vcf, table_mut, "PreCore")
table_vcf_BCP  <- searchMUT_NT(table_vcf, table_mut, "BCP")
table_vcf_RT   <- searchMUT_CODON(table_vcf, table_mut, "Domaine RT")
table_vcf_C    <- searchMUT_CODON(table_vcf, table_mut, "Core")
table_vcf_DS   <- searchMUT_CODON(table_vcf, table_mut, "Domaine S")
table_vcf_DPS1 <- searchMUT_CODON(table_vcf, table_mut, "Domaine PreS1")
table_vcf_DPS2 <- searchMUT_CODON(table_vcf, table_mut, "Domaine PreS2")
table_vcf_DHBx <- searchMUT_CODON(table_vcf, table_mut, "Domaine HBx")

filter_tbl <- function(df) if (nrow(df) > 0) subset(df, freq >= filter_threshold) else df

# Écriture des résultats bruts
write.csv2(table_vcf_BCP, file_output_BCP, row.names = F, quote = F)
write.csv2(table_vcf_PC, file_output_PC, row.names = F, quote = F)
write.csv2(table_vcf_RT, file_output_RT, row.names = F, quote = F)
write.csv2(table_vcf_DS, file_output_DS, row.names = F, quote = F)
write.csv2(table_vcf_DPS1, file_output_DPS1, row.names = F, quote = F)
write.csv2(table_vcf_DPS2, file_output_DPS2, row.names = F, quote = F)
write.csv2(table_vcf_DHBx, file_output_DHBx, row.names = F, quote = F)
write.csv2(table_vcf_C, file_output_C, row.names = F, quote = F)

# Écriture des résultats filtrés
write.csv2(filter_tbl(table_vcf_BCP), file_output_BCP_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_PC), file_output_PC_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_RT), file_output_RT_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_DS), file_output_DS_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_DPS1), file_output_DPS1_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_DPS2), file_output_DPS2_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_DHBx), file_output_DHBx_F, row.names = F, quote = F)
write.csv2(filter_tbl(table_vcf_C), file_output_C_F, row.names = F, quote = F)

# ---------------------------------------------------------------------
# 5. ÉCRITURE DU FICHIER GLOBAL *_vcf_variants.csv AVEC LA COLONNE GENE
# ---------------------------------------------------------------------
if (file_output_VARIANTS != "") {
  
  add_gene <- function(df, gene_name) {
    if (nrow(df) > 0) {
      df$GENE <- gene_name
    } else {
      df$GENE <- character()
    }
    return(df)
  }

  t_BCP  <- add_gene(table_vcf_BCP, "BCP")
  t_PC   <- add_gene(table_vcf_PC, "PreCore")
  t_RT   <- add_gene(table_vcf_RT, "Domaine RT")
  t_C    <- add_gene(table_vcf_C, "Core")
  t_DS   <- add_gene(table_vcf_DS, "Domaine S")
  t_DPS1 <- add_gene(table_vcf_DPS1, "Domaine PreS1")
  t_DPS2 <- add_gene(table_vcf_DPS2, "Domaine PreS2")
  t_DHBx <- add_gene(table_vcf_DHBx, "Domaine HBx")

  all_vars_list <- list(t_BCP, t_PC, t_RT, t_C, t_DS, t_DPS1, t_DPS2, t_DHBx)
  valid_dfs <- Filter(function(x) nrow(x) > 0, all_vars_list)
  
  if (length(valid_dfs) > 0) {
    all_cols <- unique(unlist(lapply(valid_dfs, colnames)))
    all_cols <- c("REFERENCE", "GENE", setdiff(all_cols, c("REFERENCE", "GENE")))
    
    aligned_dfs <- lapply(valid_dfs, function(df) {
      missing_cols <- setdiff(all_cols, colnames(df))
      for (col in missing_cols) {
        df[[col]] <- NA
      }
      return(df[, all_cols, drop = FALSE])
    })
    
    all_vars <- do.call(rbind, aligned_dfs)
    write.csv2(all_vars, file_output_VARIANTS, row.names = F, quote = F)
  } else {
    write_empty(file_output_VARIANTS)
  }
}

cat("✅ Screening Rscript terminé avec succès.\n")