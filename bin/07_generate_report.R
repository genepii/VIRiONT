#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 07_generate_report.R <dir_counts_04> <cutoff_percent> <pdf_out>")
}

counts_dir <- args[1]
cutoff     <- as.numeric(args[2])
pdf_out    <- args[3]

files <- list.files(counts_dir, pattern = "_real_counts\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) {
  stop(paste("Aucun fichier *_real_counts.tsv trouvé dans :", counts_dir))
}

status_colors <- c("VALIDATED" = "#2b5c8f", "REJECTED" = "#e74c3c")

pdf(pdf_out, width = 8.5, height = 11)

for (f in sort(files)) {
  sample_id <- gsub("_real_counts\\.tsv$", "", basename(f))
  
  df <- tryCatch({
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(df) || nrow(df) == 0) next
  
  colnames(df)[1:2] <- c("reference", "raw_count")
  df$raw_count <- as.numeric(df$raw_count)
  
  max_count <- max(df$raw_count, na.rm = TRUE)
  df$ratio  <- if (max_count > 0) (df$raw_count / max_count) * 100 else 0
  
  # Statut clinique dynamique
  df$status <- ifelse(df$ratio >= cutoff, "VALIDATED", "REJECTED")
  df$status <- factor(df$status, levels = c("VALIDATED", "REJECTED"))
  
  # Ordonnancement alphabétique des références
  df <- df[order(df$reference), ]
  df$reference <- factor(df$reference, levels = unique(df$reference))

  # 1. Graphe haut : Lectures brutes
  p1 <- ggplot(df, aes(x = reference, y = raw_count, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 0.65) +
    scale_fill_manual(values = status_colors, drop = FALSE) +
    geom_text(aes(label = ifelse(raw_count > 0, raw_count, "")), 
              vjust = -0.5, fontface = "bold", size = 2.8) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    labs(
      title = paste("Échantillon :", sample_id, "\n1. Distribution globale des lectures par référence"),
      x = "Référence / Sous-type",
      y = "Nombre de reads",
      fill = "Statut clinique"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8.5),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  # 2. Graphe bas : Matching ratio vs Seuil clinique
  p2 <- ggplot(df, aes(x = reference, y = ratio, fill = status)) +
    geom_bar(stat = "identity", color = "black", width = 0.65) +
    scale_fill_manual(values = status_colors, drop = FALSE) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "#c0392b", linewidth = 0.9) +
    geom_text(aes(label = ifelse(raw_count > 0, sprintf("%.1f%%", ratio), "")), 
              vjust = -0.5, fontface = "bold", size = 2.8) +
    scale_y_continuous(limits = c(0, max(110, max(df$ratio) * 1.15)), expand = c(0, 0)) +
    labs(
      title = paste0("2. Ratio d'infection relative (Seuil de validation = ", cutoff, "%)"),
      x = "Référence / Sous-type",
      y = "% / Sous-type majoritaire",
      fill = "Statut clinique"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8.5),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  # 3. Tableau Récapitulatif simple : Filtrage strict des bruits de fond (< 5 reads)
  df_table <- df[df$raw_count >= 5, c("reference", "raw_count", "ratio", "status")]
  df_table <- df_table[order(-df_table$raw_count), ]
  
  if (nrow(df_table) > 0) {
    df_table$ratio <- paste0(round(df_table$ratio, 1), " %")
    colnames(df_table) <- c("Sous-type", "Reads réels", "Ratio relatif", "Statut")
    
    bg_rows <- rep(c("#f8f9fa", "#ffffff"), length.out = nrow(df_table))
    
    table_theme <- ttheme_default(
      core = list(
        fg_params = list(fontsize = 8.5, fontface = "plain"),
        bg_params = list(fill = bg_rows, col = "#dcdde1")
      ),
      colhead = list(
        bg_params = list(fill = "#2b5c8f", col = "#dcdde1"),
        fg_params = list(col = "white", fontface = "bold", fontsize = 9)
      )
    )
    p_table <- tableGrob(head(df_table, 8), rows = NULL, theme = table_theme)
  } else {
    p_table <- textGrob("Aucun sous-type avec au moins 5 reads détecté.", 
                        gp = gpar(fontsize = 9, fontface = "italic", col = "red"))
  }

  grid.arrange(p1, p2, p_table, heights = c(2.3, 2.3, 1.4))
}

dev.off()
cat(paste("Rapport clinique PDF généré avec succès :", pdf_out, "\n"))