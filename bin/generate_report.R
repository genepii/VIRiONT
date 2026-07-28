#!/usr/bin/env Rscript

# Capture des arguments transmis par le Bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript generate_report.R <data_file.tsv> <cutoff_percentage> <output_report.pdf>")
}

data_file <- args[1]
cutoff    <- as.numeric(args[2])
pdf_out   <- args[3]

# Vérification de l'existence du fichier de données
if (!file.exists(data_file)) {
  quit(save = "no", status = 0)
}

df <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (nrow(df) == 0) {
  quit(save = "no", status = 0)
}

# Nettoyage et conversion du ratio en valeur numérique
df$ratio_num <- as.numeric(gsub("%", "", df$ratio_percent))
samplename <- df$sample[1]

# Génération du rapport PDF
pdf(pdf_out, width = 9, height = 10)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))

# Graphique 1 : Nombre de Reads par Génotype
colors_status <- ifelse(df$status == "VALIDATED", "#2b5c8f", "#e74c3c")
b1 <- barplot(df$total_reads, names.arg = df$genotype, col = colors_status,
              main = paste("Échantillon :", samplename, "\nDistribution des Reads par Génotype"),
              ylab = "Nombre de Reads", las = 1, cex.names = 0.9)
text(b1, df$total_reads / 2, labels = df$total_reads, col = "white", font = 2)
legend("topright", legend = c("VALIDATED", "REJECTED"), fill = c("#2b5c8f", "#e74c3c"), bty = "n")

# Graphique 2 : Ratio % par rapport au génotype majeur vs Seuil Clinique
b2 <- barplot(df$ratio_num, names.arg = df$genotype, col = colors_status,
              main = paste("Ratio d'Infection Relative (Seuil Clinique =", cutoff, "%)"),
              ylab = "% / Génotype Majeur", ylim = c(0, max(df$ratio_num) * 1.25), las = 1, cex.names = 0.9)
abline(h = cutoff, col = "red", lty = 2, lwd = 2)
text(b2, df$ratio_num + 4, labels = paste0(df$ratio_num, "%"), font = 2, cex = 0.9)
legend("topright", legend = c(paste("Seuil Co-inf (", cutoff, "%)", sep="")), lty = 2, col = "red", lwd = 2, bty = "n")

dev.off()