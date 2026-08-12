#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 06_generate_report.R <summary_file.tsv> <cutoff_percent> <pdf_out>")
}

data_file <- args[1]
cutoff    <- as.numeric(args[2])
pdf_out   <- args[3]

if (!file.exists(data_file)) {
  stop(paste("❌ Fichier de données introuvable :", data_file))
}

# Lecture du fichier TSV multi-échantillons / résumé
df_all <- tryCatch({
  read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) NULL)

if (is.null(df_all) || nrow(df_all) == 0) {
  cat("⚠️ Aucune donnée dans le fichier TSV, rapport non généré.\n")
  pdf(pdf_out, width = 8.5, height = 11)
  grid.newpage()
  grid.text("Aucune donnée disponible", x = 0.5, y = 0.5, gp = gpar(fontsize = 14, col = "red"))
  dev.off()
  quit(save = "no", status = 0)
}

# Harmonisation des colonnes numériques
if (!"ratio_num" %in% colnames(df_all)) {
  df_all$ratio_num <- as.numeric(gsub("%", "", as.character(df_all$ratio_percent)))
}

if (!"status" %in% colnames(df_all)) {
  df_all$status <- ifelse(df_all$ratio_num >= cutoff, "VALIDATED", "REJECTED")
}

# Fixation stricte des niveaux du facteur status
df_all$status <- factor(df_all$status, levels = c("VALIDATED", "REJECTED"))

# Palettes de couleurs
status_colors <- c("VALIDATED" = "#2b5c8f", "REJECTED" = "#e74c3c")
ratio_colors  <- c("VALIDATED" = "#41b6c4", "REJECTED" = "#e74c3c")

# Ouverture du PDF
pdf(pdf_out, width = 8.5, height = 11)

samples <- unique(df_all$sample)

for (samplename in samples) {
  df_sample <- df_all[df_all$sample == samplename, ]
  
  max_reads <- max(df_sample$total_reads, na.rm = TRUE)
  
  df_plot <- df_sample[(df_sample$total_reads >= (max_reads * 0.0001) & df_sample$total_reads >= 2) | df_sample$status == "VALIDATED", ]
  
  if (nrow(df_plot) == 0) {
    df_plot <- head(df_sample[order(-df_sample$total_reads), ], 3)
  }

  df_plot <- df_plot[order(-df_plot$total_reads), ]
  df_plot$genotype <- factor(df_plot$genotype, levels = unique(df_plot$genotype))

  # ---------------- 1. Distribution des Reads ----------------
  p1 <- ggplot(df_plot, aes(x = genotype, y = total_reads, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 0.35) +
    scale_fill_manual(values = status_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    geom_text(aes(label = total_reads), vjust = -0.5, fontface = "bold", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    guides(fill = guide_legend(override.aes = list(fill = c("#2b5c8f", "#e74c3c")))) + # <-- Force la couleur dans la légende
    labs(title = paste("Échantillon :", samplename, "\n1. Distribution des Reads par Génotype"),
         x = "Génotype Identifié", y = "Nombre de Reads", fill = "Statut") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, face = "bold", size = 10)
    )

  # ---------------- 2. Ratio % vs Seuil Clinique ----------------
  p2 <- ggplot(df_plot, aes(x = genotype, y = ratio_num, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 0.35) +
    scale_fill_manual(values = ratio_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "red", linewidth = 1) +
    geom_text(aes(label = paste0(ratio_num, "%")), vjust = -0.5, fontface = "bold", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    guides(fill = guide_legend(override.aes = list(fill = c("#41b6c4", "#e74c3c")))) + # <-- Force la couleur dans la légende
    labs(title = paste("2. Ratio d'Infection Relative (Seuil Clinique =", cutoff, "%)"),
         x = "Génotype Identifié", y = "% / Génotype Majeur", fill = "Statut") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, face = "bold", size = 10)
    )

  # ---------------- 3. Tableau Récapitulatif ----------------
  cols_req <- c("genotype", "best_cluster", "total_reads", "ratio_percent", "pident", "length", "strand", "status")
  cols_present <- intersect(cols_req, colnames(df_plot))
  table_df <- head(df_plot[, cols_present, drop = FALSE], 6)
  
  if (ncol(table_df) == 8) {
    colnames(table_df) <- c("Génotype", "Best Cluster", "Reads", "Ratio", "%Ident", "Taille", "Brin", "Statut")
  }

  table_theme <- ttheme_default(
    core = list(fg_params = list(fontsize = 8), bg_params = list(fill = c("#f1f4f9", "#ffffff"), col = NA)),
    colhead = list(bg_params = list(fill = "#2b5c8f"), fg_params = list(col = "white", fontface = "bold", fontsize = 8))
  )

  p_table <- tableGrob(table_df, rows = NULL, theme = table_theme)

  # Assemblage
  grid.arrange(p1, p2, p_table, heights = c(2.2, 2.2, 1.2))
}

dev.off()
cat(paste0("📊 Rapport PDF optimisé créé avec succès : ", pdf_out, "\n"))