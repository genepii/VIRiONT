#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript generate_report.R <validated_summary.tsv> <cutoff_percent> <output_report.pdf>")
}

data_file <- args[1]
cutoff    <- as.numeric(args[2])
pdf_out   <- args[3]

if (!file.exists(data_file)) {
  quit(save = "no", status = 0)
}

df <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (nrow(df) == 0) {
  quit(save = "no", status = 0)
}

# Extraction propre du ratio numérique
df$ratio_num <- as.numeric(gsub("%", "", df$ratio_percent))
samplename <- df$sample[1]

# Graphique 1 : Distribution des Reads
p1 <- ggplot(df, aes(x = factor(genotype), y = total_reads, fill = factor(status))) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  scale_fill_manual(values = c("VALIDATED" = "#2b5c8f", "REJECTED" = "#e74c3c")) +
  geom_text(aes(label = total_reads), vjust = -0.5, fontface = "bold", size = 3.5) +
  labs(title = paste("Échantillon :", samplename, "\nDistribution des Reads par Génotype"),
       x = "Génotype Identifié", y = "Nombre de Reads", fill = "Statut") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Graphique 2 : Ratio % vs Seuil Clinique
p2 <- ggplot(df, aes(x = factor(genotype), y = ratio_num, fill = factor(status))) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  scale_fill_manual(values = c("VALIDATED" = "#41b6c4", "REJECTED" = "#e74c3c")) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red", linewidth = 1.2) +
  geom_text(aes(label = paste0(ratio_num, "%")), vjust = -0.5, fontface = "bold", size = 3.5) +
  labs(title = paste("Ratio d'Infection Relative (Seuil Clinique =", cutoff, "%)"),
       x = "Génotype Identifié", y = "% / Génotype Majeur", fill = "Statut") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Génération du PDF
pdf(pdf_out, width = 9, height = 9)
grid.arrange(p1, p2, nrow = 2)
dev.off()

cat(paste0("   📊 Rapport PDF généré avec succès : ", basename(pdf_out), "\n"))