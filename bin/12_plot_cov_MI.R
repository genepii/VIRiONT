#!/usr/bin/env Rscript
# =====================================================================
# VIRiONT_NF - Rendu Dynamique des Profils de Couverture (12_COVERAGE)
# Support : VHB WG (3,2 kb), VHB POL (~1 kb), VHD WG (1,7 kb), VHD R0 (400 pb)
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
})

argv <- commandArgs(TRUE)
if (length(argv) < 2) {
  stop("Usage: Rscript 12_plot_cov_MI.R <cov_sum_file> <output_pdf> [virus_name] [tech_name]")
}

cov_file   <- as.character(argv[1])
output_pdf <- as.character(argv[2])
virus_name <- if (length(argv) >= 3) as.character(argv[3]) else "VIRUS"
tech_name  <- if (length(argv) >= 4) as.character(argv[4]) else "WG"

if (!file.exists(cov_file) || file.info(cov_file)$size == 0) {
  cat("Fichier de couverture introuvable ou vide. Graphique non généré.\n")
  quit(save = "no", status = 0)
}

covdata <- tryCatch({
  read.table(cov_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) NULL)

if (is.null(covdata) || nrow(covdata) == 0) {
  cat("Aucune donnée de couverture valide à lire.\n")
  quit(save = "no", status = 0)
}

colnames(covdata)[1:6] <- c("REF", "POS", "COV", "SAMPLE", "GENOTYPE", "METHOD")
covdata$POS <- as.numeric(covdata$POS)
covdata$COV <- as.numeric(covdata$COV)

covdata <- covdata[!is.na(covdata$POS) & !is.na(covdata$COV), ]
if (nrow(covdata) == 0) {
  cat("Aucune coordonnée valide après nettoyage.\n")
  quit(save = "no", status = 0)
}

# 1. En-tête simple : Échantillon et Sous-type uniquement (médiane supprimée)
covdata$LABEL <- paste0(covdata$SAMPLE, "\nSous-type : ", covdata$GENOTYPE)

# 2. Recadrage automatique de l'axe X pour amplicons ciblés (POL / R0)
pos_covered <- covdata$POS[covdata$COV >= 10]
total_ref_len <- max(covdata$POS)

if (tech_name == "POL" || (length(pos_covered) > 0 && (max(pos_covered) - min(pos_covered)) < (0.6 * total_ref_len))) {
  if (length(pos_covered) > 0) {
    min_x <- max(1, min(pos_covered) - 60)
    max_x <- min(total_ref_len, max(pos_covered) + 60)
    covdata <- covdata[covdata$POS >= min_x & covdata$POS <= max_x, ]
  }
}

# 3. Libellé de l'axe X
x_label_title <- if (tech_name == "POL") {
  "Position Génomique EcoRI (pb) [Région RT / POL ciblée]"
} else if (tech_name == "R0") {
  "Position Génomique (pb) [Amplicon R0 (~400 pb)]"
} else {
  "Position Génomique (pb)"
}

plot_title <- sprintf("VIRiONT_NF – Profils de Couverture Génomique [%s %s]", virus_name, tech_name)

# 4. Tracé graphique
covplot <- ggplot(covdata, aes(x = POS, y = COV)) +
  geom_area(fill = "#2b5c8f", alpha = 0.35) +
  geom_line(color = "#1b3a5a", linewidth = 0.7) +
  scale_y_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = " "),
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE, big.mark = " ")
  ) +
  labs(
    title = plot_title,
    x = x_label_title,
    y = "Profondeur de Couverture (X)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#2b5c8f", color = NA),
    strip.text = element_text(face = "bold", color = "white", size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
  ) +
  facet_wrap(~ LABEL, scales = "free", ncol = 2)

num_facets  <- length(unique(covdata$LABEL))
num_rows    <- ceiling(num_facets / 2)
calc_height <- max(4.5, 3.2 * num_rows)

ggsave(
  filename  = output_pdf,
  plot      = covplot,
  width     = 11,
  height    = calc_height,
  limitsize = FALSE
)
cat("Graphique de couverture généré avec succès :", output_pdf, "\n")