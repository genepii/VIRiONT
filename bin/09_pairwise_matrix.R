#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT 09: GENERATION DE LA MATRICE PAIRWISE & HEATMAP DE CONTAMINATION
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)  # <-- TRÈS IMPORTANT : Évite l'erreur "could not find function ggplot"
})

# 1. Chargement de tous les consensus renommés
fasta_files <- list.files(".", pattern = ".*_Geno_.*\\.fasta$", full.names = TRUE)

if (length(fasta_files) < 2) {
  cat("Pas assez de fichiers consensus validés pour générer une matrice pairwise.\n")
  quit(save = "no", status = 0)
}

read_fasta <- function(file) {
  lines <- readLines(file)
  header <- gsub("^>", "", lines[1])
  seq <- paste(lines[-1], collapse = "")
  return(list(header = header, seq = toupper(seq)))
}

seqs <- list()
for (f in fasta_files) {
  s <- read_fasta(f)
  seqs[[s$header]] <- s$seq
}

samples <- sort(names(seqs))
n <- length(samples)

matrix_rows <- list()
plot_data <- data.frame()

iupac_map <- list(
  'R'=c('A','G'), 'Y'=c('C','T'), 'S'=c('G','C'), 'W'=c('A','T'),
  'K'=c('G','T'), 'M'=c('A','C'), 'B'=c('C','G','T'), 'D'=c('A','G','T'),
  'H'=c('A','C','T'), 'V'=c('A','C','G')
)

for (i in 1:n) {
  for (j in 1:n) {
    s1_id <- samples[i]
    s2_id <- samples[j]
    
    seq1 <- unlist(strsplit(seqs[[s1_id]], ""))
    seq2 <- unlist(strsplit(seqs[[s2_id]], ""))
    
    len <- min(length(seq1), length(seq2))
    
    mismatches <- 0
    desc_list <- c()
    
    if (i != j) {
      for (pos in 1:len) {
        b1 <- seq1[pos]
        b2 <- seq2[pos]
        
        if (b1 != b2 && b1 != "N" && b2 != "N" && b1 != "-" && b2 != "-") {
          mismatches <- mismatches + 1
          
          is_iupac <- b1 %in% names(iupac_map) || b2 %in% names(iupac_map)
          prefix <- if (is_iupac) "IUPAC:" else "pos:"
          
          desc_list <- c(desc_list, paste0(prefix, pos, "/1:", b1, "/2:", b2))
        }
      }
    }
    
    desc_str <- paste(desc_list, collapse = ";")
    
    if (j >= i) {
      plot_data <- rbind(plot_data, data.frame(
        comp1 = s1_id,
        comp2 = s2_id,
        count = mismatches,
        stringsAsFactors = FALSE
      ))
    }
    
    matrix_rows[[length(matrix_rows) + 1]] <- data.frame(
      comp1 = s1_id,
      comp2 = s2_id,
      count = mismatches,
      description = desc_str,
      stringsAsFactors = FALSE
    )
  }
}

# Export TSV
final_matrix_df <- do.call(rbind, matrix_rows)
write.table(final_matrix_df, "matrix_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Export PDF avec gradient contrasté (Nouvelles couleurs & échelle)
plot_data$comp1 <- factor(plot_data$comp1, levels = samples)
plot_data$comp2 <- factor(plot_data$comp2, levels = rev(samples))

p <- ggplot(plot_data, aes(x = comp1, y = comp2, fill = count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = count), color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#2b83ba", "#abdda4", "#fdae61", "#d7191c"),
    trans = "sqrt",
    breaks = c(0, 10, 50, 100, 300, 700),
    name = "Mismatches"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("matrix_comp.pdf", plot = p, width = 9, height = 9)