#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
fasta_file <- ifelse(length(args) >= 1, args[1], "validated_consensus_all.fasta")

if (!file.exists(fasta_file)) {
  cat("⚠️ Aucun fichier consensus validé trouvé pour la matrice.\n")
  quit(save = "no", status = 0)
}

# Lecture multi-FASTA
lines <- readLines(fasta_file)
seqs <- list()
current_hdr <- ""
current_seq <- ""

for (line in lines) {
  if (grepl("^>", line)) {
    if (current_hdr != "") {
      seqs[[current_hdr]] <- toupper(current_seq)
    }
    current_hdr <- gsub("^>", "", line)
    current_seq <- ""
  } else {
    current_seq <- paste0(current_seq, gsub("\\s+", "", line))
  }
}
if (current_hdr != "") {
  seqs[[current_hdr]] <- toupper(current_seq)
}

samples <- sort(names(seqs))
n <- length(samples)

if (n < 2) {
  cat("Moins de 2 échantillons validés : matrice pairwise non générée.\n")
  quit(save = "no", status = 0)
}

matrix_rows <- list()
plot_data <- data.frame()

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
          desc_list <- c(desc_list, paste0("pos:", pos, "/1:", b1, "/2:", b2))
        }
      }
    }
    
    if (j >= i) {
      plot_data <- rbind(plot_data, data.frame(
        comp1 = s1_id, comp2 = s2_id, count = mismatches, stringsAsFactors = FALSE
      ))
    }
    
    matrix_rows[[length(matrix_rows) + 1]] <- data.frame(
      comp1 = s1_id, comp2 = s2_id, count = mismatches,
      description = paste(desc_list, collapse = ";"), stringsAsFactors = FALSE
    )
  }
}

# Export TSV
write.table(do.call(rbind, matrix_rows), "matrix_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot Heatmap
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
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("matrix_comp.pdf", plot = p, width = 9, height = 9)