#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(seqinr)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  library(optparse)
})

option_list = list(
  make_option(c("--input"), type="character", default=NULL, help="input aln file", metavar="character"),
  make_option(c("--matrix"), type="character", default=NULL, help="output matrix pdf", metavar="character"),
  make_option(c("--table"), type="character", default=NULL, help="output table tsv", metavar="character"),
  make_option(c("--countSNP"), type="character", default=NULL, help="count SNP as a discordance", metavar="character"),
  make_option(c("--countDEL"), type="character", default=NULL, help="count deletion stretches as a discordance", metavar="character"),
  make_option(c("--countN"), type="character", default=NULL, help="count N stretches as a discordance", metavar="character"),
  make_option(c("--countIUPAC"), type="character", default=NULL, help="count IUPAC bases as a discordance", metavar="character"),
  make_option(c("--noheader"), type="character", default=NULL, help="do not write header in TSV", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ------------------ Helpers ------------------
# Nettoyage pour l'affichage uniquement (axes + TSV)
clean_display_name <- function(x) {
  x <- str_replace_all(x, "^barcode_", "")
  x <- str_replace_all(x, "_ONT0\\.08$", "")
  x <- str_trim(x)
  return(x)
}

# Règle demandée :
# - si l'ID finit par "_ONT0.08", génotype = token juste avant (1–3 alphanum)
# - sinon, génotype = dernier token (1–3 alphanum)
extract_genotype <- function(x) {
  if (grepl("_ONT0\\.08$", x)) {
    m <- str_match(x, "_([A-Za-z0-9]{1,3})_ONT0\\.08$")
    g <- if (!is.na(m[1,2])) m[1,2] else NA_character_
  } else {
    m <- str_match(x, "_([A-Za-z0-9]{1,3})$")
    g <- if (!is.na(m[1,2])) m[1,2] else NA_character_
  }
  return(g)
}

# Pour tri naturel : convertir g en (lettre, has_num, num)
split_genotype <- function(g) {
  if (is.na(g) || length(g) == 0) {
    return(list(letter=NA_character_, has_num=NA_integer_, num=NA_real_))
  }
  letter <- toupper(substr(g, 1, 1))
  rest <- if (nchar(g) > 1) substr(g, 2, nchar(g)) else ""
  if (!grepl("^[A-Z]$", letter)) {
    # pas une lettre en tête -> on met en fin (NA)
    return(list(letter=NA_character_, has_num=NA_integer_, num=NA_real_))
  }
  if (rest == "") {
    return(list(letter=letter, has_num=0L, num=NA_real_))  # ex: "E"
  } else if (grepl("^[0-9]{1,2}$", rest)) {
    return(list(letter=letter, has_num=1L, num=as.numeric(rest)))  # ex: "A10"
  } else {
    # token atypique -> en fin
    return(list(letter=NA_character_, has_num=NA_integer_, num=NA_real_))
  }
}

# ---------------- Lecture & Préparation ----------------
allseq <- read.fasta(opt$input)
orig_samples <- names(allseq)

# mapping noms originels -> affichage + génotype
display_names <- setNames(clean_display_name(orig_samples), orig_samples)
genotypes     <- setNames(vapply(orig_samples, extract_genotype, FUN.VALUE = character(1)),
                          orig_samples)

# ---------------- Comparaisons pair-à-pair ----------------
list_result <- vector("list", length(orig_samples) * length(orig_samples))
cmpt <- 1

for (sample1 in orig_samples) {
  aln1 <- toupper(allseq[[sample1]])
  for (sample2 in orig_samples) {
    aln2 <- toupper(allseq[[sample2]])

    nb_diff <- 0
    list_diff <- character(0)
    k <- 1
    L <- length(aln1)

    while (k <= L) {
      if (aln1[k] != aln2[k]) {
        if (aln1[k] == "-") {
          newk <- k
          while ((newk + 1) <= L && aln1[newk + 1] == "-") newk <- newk + 1
          if (opt$countDEL == "TRUE") {
            list_diff <- c(list_diff, paste0("1:DEL", k, "-", newk))
            nb_diff <- nb_diff + 1
          }
          k <- newk + 1
          next
        } else if (aln2[k] == "-") {
          newk <- k
          while ((newk + 1) <= L && aln2[newk + 1] == "-") newk <- newk + 1
          if (opt$countDEL == "TRUE") {
            list_diff <- c(list_diff, paste0("2:DEL", k, "-", newk))
            nb_diff <- nb_diff + 1
          }
          k <- newk + 1
          next
        } else if (aln1[k] == "N" && aln2[k] != "N") {
          newk <- k
          while (newk <= L && aln1[newk] == "N" && aln2[newk] != "N") newk <- newk + 1
          if (opt$countN == "TRUE") {
            list_diff <- c(list_diff, paste0("1:Nstretch", k, "-", newk))
            nb_diff <- nb_diff + 1
          }
          k <- newk + 1
          next
        } else if (aln2[k] == "N" && aln1[k] != "N") {
          newk <- k
          while (newk <= L && aln2[newk] == "N" && aln1[newk] != "N") newk <- newk + 1
          if (opt$countN == "TRUE") {
            list_diff <- c(list_diff, paste0("2:Nstretch", k, "-", newk))
            nb_diff <- nb_diff + 1
          }
          k <- newk
          next
        } else if (!(aln1[k] %in% c("A","C","G","T","N"))) {
          if (opt$countIUPAC == "TRUE") {
            list_diff <- c(list_diff, paste0("IUPAC:", k, "/1:", aln1[k], "/2:", aln2[k]))
            nb_diff <- nb_diff + 1
          }
          k <- k + 1
          next
        } else if (!(aln2[k] %in% c("A","C","G","T","N"))) {
          if (opt$countIUPAC == "TRUE") {
            list_diff <- c(list_diff, paste0("IUPAC:", k, "/1:", aln1[k], "/2:", aln2[k]))
            nb_diff <- nb_diff + 1
          }
          k <- k + 1
          next
        } else {
          if (opt$countSNP == "TRUE") {
            list_diff <- c(list_diff, paste0("pos:", k, "/1:", aln1[k], "/2:", aln2[k]))
            nb_diff <- nb_diff + 1
          }
          k <- k + 1
          next
        }
      } else {
        k <- k + 1
        next
      }
    }

    # comp1/comp2 = NOMS NETTOYÉS (affichage)
    vec_result <- c(display_names[[sample1]],
                    display_names[[sample2]],
                    as.character(nb_diff),
                    paste0(list_diff, collapse = ";"))
    list_result[[cmpt]] <- vec_result
    cmpt <- cmpt + 1
  }
}

table_merge <- as.data.frame(do.call("rbind", list_result), stringsAsFactors = FALSE)
colnames(table_merge) <- c("comp1","comp2","count","description")
table_merge$count <- as.numeric(table_merge$count)

# ---------------- Ordonnancement naturel par génotype ----------------
# mapping display -> génotype (via noms originels)
display_to_genotype <- setNames(genotypes, clean_display_name(names(genotypes)))

all_display <- sort(unique(c(table_merge$comp1, table_merge$comp2)))
geno_for_display <- display_to_genotype[all_display]

# Clés de tri naturelles
keys <- lapply(as.list(geno_for_display), split_genotype)
letters_vec <- vapply(keys, function(k) if (is.null(k$letter) || is.na(k$letter)) NA_character_ else k$letter, character(1))
hasnum_vec  <- vapply(keys, function(k) if (is.null(k$has_num) || is.na(k$has_num)) NA_integer_ else k$has_num, integer(1))
num_vec     <- vapply(keys, function(k) if (is.null(k$num) || is.na(k$num)) NA_real_ else k$num, numeric(1))

na_mask <- is.na(letters_vec)

# Tri :
# 1) NA (génotype manquant) à la fin
# 2) par lettre A..Z
# 3) sans numéro (0) avant numérotés (1)
# 4) par numéro croissant
# 5) puis par nom d'affichage
ord <- order(na_mask, letters_vec, hasnum_vec, num_vec, all_display)
ordered_levels <- all_display[ord]

# ---------------- Écriture TSV ----------------
if (is.null(opt$noheader)) {
  write.table(table_merge, opt$table, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  write.table(table_merge, opt$table, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ---------------- Matrice & Heatmap ----------------
table_merge$comp1 <- factor(table_merge$comp1, levels = ordered_levels)
table_merge$comp2 <- factor(table_merge$comp2, levels = ordered_levels)

mat <- acast(table_merge, comp1 ~ comp2, value.var = "count")
mat[lower.tri(mat)] <- NA
tbl <- melt(mat, na.rm = TRUE)
colnames(tbl) <- c("comp1","comp2","count")

final_plot <- ggplot(data = tbl, aes(comp1, comp2, fill = count)) +
  xlab("") + ylab("") +
  geom_tile(color = "white", show.legend = FALSE) +
  geom_text(aes(label = count), color = "black", size = 3) +
  theme(axis.text.x = element_text(vjust = 1, size = 6, hjust = 1, angle = 45),
        axis.text.y = element_text(vjust = 0.5, size = 6, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 8)) +
  scale_fill_gradient(low = "blue", high = "red")

pdf(opt$matrix)
final_plot
dev.off()
