#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 5_generate_report.R <validated_summary.tsv> <cutoff_percent> <pdf_out>")
}

data_file <- args[1]
cutoff    <- as.numeric(args[2])
pdf_out   <- args[3]

if (!file.exists(data_file)) {
  stop(paste("❌ Fichier de données introuvable :", data_file))
}

# Lecture des données de génotypage
df <- tryCatch({
  read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) NULL)

if (is.null(df) || nrow(df) == 0) {
  cat("⚠️ Aucune donnée dans le fichier TSV, rapport non généré.\n")
  quit(save = "no", status = 0)
}

df$status    <- factor(df$status, levels = c("VALIDATED", "REJECTED"))
df$ratio_num <- as.numeric(gsub("%", "", df$ratio_percent))
samplename   <- df$sample[1]

status_colors <- c("VALIDATED" = "#2b5c8f", "REJECTED" = "#e74c3c")
ratio_colors  <- c("VALIDATED" = "#41b6c4", "REJECTED" = "#e74c3c")

# Graphique 1 : Distribution des Reads
p1 <- ggplot(df, aes(x = factor(genotype), y = total_reads, fill = status)) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  geom_text(aes(label = total_reads), vjust = -0.5, fontface = "bold", size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = paste("Échantillon :", samplename, "\nDistribution des Reads par Génotype"),
       x = "Génotype Identifié", y = "Nombre de Reads", fill = "Statut") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Graphique 2 : Ratio % vs Seuil Clinique
p2 <- ggplot(df, aes(x = factor(genotype), y = ratio_num, fill = status)) +
  geom_bar(stat = "identity", color = "black", width = 0.4) +
  scale_fill_manual(values = ratio_colors, drop = FALSE) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red", linewidth = 1.2) +
  geom_text(aes(label = paste0(ratio_num, "%")), vjust = -0.5, fontface = "bold", size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = paste("Ratio d'Infection Relative (Seuil Clinique =", cutoff, "%)"),
       x = "Génotype Identifié", y = "% / Génotype Majeur", fill = "Statut") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Tableau récapitulatif
table_df <- df[, c("genotype", "best_cluster", "total_reads", "ratio_percent", "pident", "length", "strand", "status")]
colnames(table_df) <- c("Génotype", "Best Cluster", "Reads", "Ratio", "%Ident", "Taille", "Brin", "Statut")

table_theme <- ttheme_default(
  core = list(bg_params = list(fill = c("#f1f4f9", "#ffffff"), col = NA)),
  colhead = list(bg_params = list(fill = "#2b5c8f"), fg_params = list(col = "white", fontface = "bold"))
)

p_table <- tableGrob(table_df, rows = NULL, theme = table_theme)

# Ouverture du PDF multi-pages
pdf(pdf_out, width = 8.5, height = 11)

# Page 1 : Vue d'ensemble graphique & tableau
grid.arrange(p1, p2, p_table, heights = c(2.5, 2.5, 1.5))

# Page 2 : Détail complet de la distribution des filtres et co-infections
p_detail <- ggplot(df, aes(x = factor(genotype), y = ratio_num, fill = genotype)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  facet_wrap(~ status, scales = "free_y") +
  labs(title = paste("Détail par Statut de Validation -", samplename),
       x = "Génotype", y = "Ratio (%)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

grid.arrange(p_detail, nrow = 1)

dev.off()

cat(paste0("📊 Rapport PDF multi-pages créé avec succès : ", pdf_out, "\n"))