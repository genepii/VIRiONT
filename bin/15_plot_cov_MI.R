#!/usr/bin/env Rscript

# =====================================================================
# VIRiONT_NF - Rendu des Profils de Couverture Génomique (10_COVERAGE)
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

argv <- commandArgs(TRUE)

if (length(argv) < 2) {
  stop("Usage: Rscript 15_plot_cov_MI.R <cov_sum_file> <output_pdf>")
}

cov_file   <- as.character(argv[1])
output_pdf <- as.character(argv[2])

if (!file.exists(cov_file) || file.info(cov_file)$size == 0) {
  cat("⚠️ Fichier de couverture introuvable ou vide. Graphique non généré.\n")
  quit(save = "no", status = 0)
}

covdata <- tryCatch({
  read.table(cov_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) NULL)

if (is.null(covdata) || nrow(covdata) == 0) {
  cat("⚠️ Aucune donnée de couverture valide à lire.\n")
  quit(save = "no", status = 0)
}

colnames(covdata)[1:5] <- c("REF", "POS", "COV", "SAMPLE", "METHOD")

covdata$POS <- as.numeric(covdata$POS)
covdata$COV <- as.numeric(covdata$COV)

# Filtrage des lignes valides
covdata <- covdata[!is.na(covdata$POS) & !is.na(covdata$COV), ]

if (nrow(covdata) == 0) {
  cat("⚠️ Aucune coordonnée valide après nettoyage.\n")
  quit(save = "no", status = 0)
}

# Nettoyage du titre de la référence / génotype
covdata$REF_CLEAN <- sub(".*_([^_]+)$", "\\1", covdata$REF)

# Génération du graphique ggplot
covplot <- ggplot(covdata, aes(x = POS, y = COV)) +
  geom_line(aes(color = REF_CLEAN), show.legend = FALSE, alpha = 0.8, linewidth = 0.6) +
  scale_color_discrete() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "VIRiONT_NF - Profils de Couverture Génomique",
    x = "Position Génomique (pb)",
    y = "Profondeur de Couverture (X)"
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

# Calcul dynamique de la hauteur du fichier PDF
num_facets  <- length(unique(paste(covdata$SAMPLE, covdata$REF_CLEAN)))
num_rows    <- ceiling(num_facets / 2)
calc_height <- max(4, 3 * num_rows)

ggsave(
  filename  = output_pdf,
  plot      = covplot,
  width     = 10,
  height    = calc_height,
  limitsize = FALSE
)

cat("✅ Graphique de couverture généré avec succès :", output_pdf, "\n")