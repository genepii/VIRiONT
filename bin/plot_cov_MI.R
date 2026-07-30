#!/usr/bin/env Rscript

# =====================================================================
# VIRiONT V2 - Rendu des Profils de Couverture Génomique (10_COVERAGE)
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

argv <- commandArgs(TRUE)

if (length(argv) < 2) {
  stop("Usage: Rscript plot_cov_MI.R <cov_sum_file> <output_pdf>")
}

cov_file   <- as.character(argv[1])
output_pdf <- as.character(argv[2])

if (!file.exists(cov_file) || file.info(cov_file)$size == 0) {
  stop(paste("❌ Fichier de couverture introuvable ou vide :", cov_file))
}

covdata <- read.table(cov_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(covdata)[1:5] <- c("REF", "POS", "COV", "SAMPLE", "METHOD")

covdata$POS <- as.numeric(covdata$POS)
covdata$COV <- as.numeric(covdata$COV)

# ---------------------------------------------------------------------
# NETTOYAGE DU TITRE DE LA RÉFÉRENCE / GÉNOTYPE
# Remplace ex: "barcode_01_26104456601_B4" par juste "B4"
# ---------------------------------------------------------------------
covdata$REF_CLEAN <- sub(".*_([^_]+)$", "\\1", covdata$REF)

covplot <- ggplot(covdata, aes(x = POS, y = COV)) +
  geom_line(aes(color = REF_CLEAN), show.legend = FALSE, alpha = 0.75, linewidth = 0.6) +
  scale_color_manual(values = c("#D9534F", "#2B6CB0", "#38A169", "#D69E2E", "#805AD5")) +
  labs(
    title = "VIRiONT V2 - Profils de Couverture Génomique",
    x = "Genomic position",
    y = "Coverage"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#E2E8F0", color = NA),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  facet_wrap(vars(SAMPLE, REF_CLEAN), scales = "free", ncol = 2)

num_samples <- length(unique(covdata$SAMPLE))
calc_height <- max(4, 2.5 * num_samples)

ggsave(
  filename = output_pdf,
  plot = covplot,
  path = NULL,
  width = 10,
  height = calc_height,
  limitsize = FALSE
)

cat("✅ Graphique de couverture généré avec succès (Titres nettoyés) :", output_pdf, "\n")