#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 08_generate_report.R <summary_file.tsv> <cutoff_percent> <pdf_out>")
}

data_file <- args[1]
cutoff    <- as.numeric(args[2])
pdf_out   <- args[3]

if (!file.exists(data_file)) {
  stop(paste("❌ Fichier introuvable :", data_file))
}

df_all <- tryCatch({
  read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) NULL)

if (is.null(df_all) || nrow(df_all) == 0) {
  pdf(pdf_out, width = 8.5, height = 11)
  grid.newpage()
  grid.text("Aucune donnée disponible", x = 0.5, y = 0.5, gp = gpar(fontsize = 14, col = "red"))
  dev.off()
  quit(save = "no", status = 0)
}

# 1. FILTRE STRICT : Suppression définitive des Unknown et non-assignés
df_all <- df_all[!toupper(df_all$genotype) %in% c("UNKNOWN", "NO_MATCH", "NONE", ""), ]

if (nrow(df_all) == 0) {
  pdf(pdf_out, width = 8.5, height = 11)
  grid.newpage()
  grid.text("Aucun génotype valide identifié pour ce run", x = 0.5, y = 0.5, gp = gpar(fontsize = 14, col = "black", fontface = "bold"))
  dev.off()
  quit(save = "no", status = 0)
}

if (!"ratio_num" %in% colnames(df_all)) {
  df_all$ratio_num <- as.numeric(gsub("%", "", as.character(df_all$ratio_percent)))
}

if (!"status" %in% colnames(df_all)) {
  df_all$status <- ifelse(df_all$ratio_num >= cutoff, "VALIDATED", "REJECTED")
}

df_all$status <- factor(df_all$status, levels = c("VALIDATED", "REJECTED"))

status_colors <- c("VALIDATED" = "#2b5c8f", "REJECTED" = "#e74c3c")
ratio_colors  <- c("VALIDATED" = "#41b6c4", "REJECTED" = "#e74c3c")

pdf(pdf_out, width = 8.5, height = 11)

samples <- unique(df_all$sample)

for (samplename in samples) {
  df_sample <- df_all[df_all$sample == samplename, ]
  if (nrow(df_sample) == 0) next
  
  df_plot <- df_sample[order(-df_sample$total_reads), ]
  
  # Ordonnancement sur le génotype unique
  df_plot$genotype <- factor(df_plot$genotype, levels = unique(df_plot$genotype))

  # 1. Distribution des Lectures
  p1 <- ggplot(df_plot, aes(x = genotype, y = total_reads, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 0.4) +
    scale_fill_manual(values = status_colors, drop = FALSE) +
    geom_text(aes(label = total_reads), vjust = -0.5, fontface = "bold", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = paste("Échantillon :", samplename, "\n1. Distribution des Reads réels par Génotype"),
         x = "Génotype Identifié", y = "Nombre de Reads", fill = "Statut") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  # 2. Ratio % vs Seuil Clinique
  p2 <- ggplot(df_plot, aes(x = genotype, y = ratio_num, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 0.4) +
    scale_fill_manual(values = ratio_colors, drop = FALSE) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "red", linewidth = 0.9) +
    geom_text(aes(label = paste0(round(ratio_num, 1), "%")), vjust = -0.5, fontface = "bold", size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)), limits = c(0, max(105, max(df_plot$ratio_num) * 1.15))) +
    labs(title = paste("2. Ratio d'Infection Relative (Seuil Clinique =", cutoff, "%)"),
         x = "Génotype Identifié", y = "% / Génotype Majeur", fill = "Statut") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(face = "bold", size = 10),
      legend.position = "right"
    )

  # 3. Tableau Récapitulatif
  cols_base <- c("genotype", "total_reads", "ratio_percent", "pident", "length", "strand", "status")
  cols_to_use <- intersect(cols_base, colnames(df_plot))
  table_df <- df_plot[, cols_to_use, drop = FALSE]
  
  if (ncol(table_df) == 7) {
    colnames(table_df) <- c("Génotype", "Reads", "Ratio (%)", "%Ident", "Taille (bp)", "Brin", "Statut")
  }

  table_theme <- ttheme_default(
    core = list(fg_params = list(fontsize = 8), bg_params = list(fill = c("#f1f4f9", "#ffffff"), col = NA)),
    colhead = list(bg_params = list(fill = "#2b5c8f"), fg_params = list(col = "white", fontface = "bold", fontsize = 8))
  )

  p_table <- tableGrob(head(table_df, 6), rows = NULL, theme = table_theme)

  grid.arrange(p1, p2, p_table, heights = c(2.2, 2.2, 1.2))
}

dev.off()
cat(paste0("📊 Rapport PDF généré avec succès : ", pdf_out, "\n"))