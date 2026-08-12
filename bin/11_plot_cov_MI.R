#!/usr/bin/env Rscript

# =====================================================================
# VIRiONT_NF - Rendu des Profils de Couverture Génomique (11_COVERAGE)
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

argv <- commandArgs(TRUE)

if (length(argv) < 2) {
  stop("Usage: Rscript 15_plot_cov_MI.R <cov_sum_file> <output_pdf> [summary_tsv]")
}

cov_file    <- as.character(argv[1])
output_pdf  <- as.character(argv[2])
summary_tsv <- ifelse(length(argv) >= 3, as.character(argv[3]), "")

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

# Chargement optionnel des génotypes depuis SUMMARY_Multi_Infection.tsv
geno_map <- list()
if (summary_tsv != "" && file.exists(summary_tsv)) {
  sum_data <- tryCatch({
    read.table(summary_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(sum_data) && "best_cluster" %in% colnames(sum_data) && "genotype" %in% colnames(sum_data)) {
    for (i in seq_len(nrow(sum_data))) {
      geno_map[[ as.character(sum_data$best_cluster[i]) ]] <- as.character(sum_data$genotype[i])
    }
  }
}

# Association propre du génotype
covdata$GENOTYPE <- sapply(covdata$REF, function(r) {
  if (r %in% names(geno_map)) {
    return(geno_map[[r]])
  } else {
    sub(".*_([^_]+)$", "\\1", r)
  }
})

# Formatage propre sans répétitions
covdata$FACET_TITLE <- paste0(covdata$SAMPLE, "\nGénotype : ", covdata$GENOTYPE)

# Génération du graphique ggplot
covplot <- ggplot(covdata, aes(x = POS, y = COV)) +
  geom_line(aes(color = GENOTYPE), show.legend = FALSE, alpha = 0.8, linewidth = 0.6) +
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
    strip.background = element_rect(fill = "#2b5c8f", color = NA),
    strip.text = element_text(face = "bold", color = "white", size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  facet_wrap(~ FACET_TITLE, scales = "free", ncol = 2)

# Calcul dynamique de la hauteur du fichier PDF
num_facets  <- length(unique(covdata$FACET_TITLE))
num_rows    <- ceiling(num_facets / 2)
calc_height <- max(4, 3.2 * num_rows)

ggsave(
  filename  = output_pdf,
  plot      = covplot,
  width     = 10,
  height    = calc_height,
  limitsize = FALSE
)

cat("✅ Graphique propre généré avec succès :", output_pdf, "\n")