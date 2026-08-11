#!/usr/bin/env Rscript

# =====================================================================
# VIRiONT V2 - Rendu Arbre Rectangulaire, Midpoint Rooting & Contamination
# =====================================================================

suppressPackageStartupMessages({
  library(ape)
})

argv <- commandArgs(TRUE)

if (length(argv) < 2) {
  stop("Usage: Rscript 13_plot_tree.R <treefile> <output_pdf>")
}

tree_file  <- as.character(argv[1])
output_pdf <- as.character(argv[2])

if (!file.exists(tree_file)) {
  stop(paste("❌ Fichier arbre introuvable :", tree_file))
}

# 1. Lecture de l'arbre au format Newick
tree <- read.tree(tree_file)

if (is.null(tree) || length(tree$tip.label) == 0) {
  cat("⚠️ Fichier d'arbre vide ou illisible.\n")
  quit(save = "no", status = 0)
}

num_tips <- length(tree$tip.label)

# 2. ENRACINEMENT AUTOMATIQUE PAR POINT MÉDIAN (MIDPOINT ROOTING FIX)
cat("🌱 Enracinement automatique par Point Médian...\n")

tree <- tryCatch({
  if (requireNamespace("phangorn", quietly = TRUE)) {
    phangorn::midpoint(tree)
  } else {
    node_depths <- node.depth.edgelength(tree)
    internal_node_depths <- node_depths[(num_tips + 1):length(node_depths)]
    target_node <- num_tips + which.max(internal_node_depths)
    root(tree, node = target_node, resolve.root = TRUE)
  }
}, error = function(e) {
  cat("⚠️ Avertissement lors de l'enracinement point médian, conservation de l'arbre brut.\n")
  return(tree)
})

tree <- ladderize(tree)

# 3. Identification des échantillons du run vs Références
tip_labels <- tree$tip.label
is_sample  <- grepl("^barcode_", tip_labels)

# 4. Métrologie : Détection des contaminations inter-échantillons (distance cophenetique == 0)
dist_matrix <- cophenetic(tree)
sample_indices <- which(is_sample)

cat("=====================================================================\n")
cat(" 🔍 DÉTECTION DES CONTAMINATIONS INTER-ÉCHANTILLONS (Distance = 0)\n")
cat("=====================================================================\n")

contamination_found <- FALSE
if (length(sample_indices) > 1) {
  for (i in 1:(length(sample_indices) - 1)) {
    for (j in (i + 1):length(sample_indices)) {
      idx1 <- sample_indices[i]
      idx2 <- sample_indices[j]
      s1   <- tip_labels[idx1]
      s2   <- tip_labels[idx2]
      d    <- dist_matrix[idx1, idx2]
      
      if (d == 0) {
        cat(sprintf("⚠️ ALERTE CONTAMINATION SUSPECTE : %s ET %s SONT 100%% IDENTIQUES !\n", s1, s2))
        contamination_found <- TRUE
      }
    }
  }
}

if (!contamination_found) {
  cat("✅ Aucune identité à 100% (distance = 0) détectée entre les échantillons du run.\n")
}
cat("=====================================================================\n")

# 5. Stylisation visuelle dynamique
tip_colors <- ifelse(is_sample, "#D9534F", "#2B6CB0") # Rouge pour le run, Bleu pour les refs
tip_cex    <- ifelse(is_sample, 0.85, 0.65)
font_type  <- ifelse(is_sample, 2, 1)

# 6. Dimensionnement dynamique de la hauteur du PDF
pdf_height <- max(10, num_tips * 0.25)

pdf(output_pdf, width = 11, height = pdf_height)
par(mar = c(5, 2, 4, 2))

max_depth <- max(node.depth.edgelength(tree))

plot(
  tree,
  type = "phylogram",
  show.tip.label = TRUE,
  tip.color = tip_colors,
  cex = tip_cex,
  font = font_type,
  edge.width = 1.2,
  no.margin = FALSE,
  x.lim = c(0, max_depth * 1.4),
  main = "VIRiONT_NF - Arbre Phylogénétique Global"
)

# Positionnement de la barre d'échelle
add.scale.bar(x = 0, y = 1, cex = 0.8, lwd = 1.5)

# Légende en haut à droite
legend(
  "topright",
  legend = c("Consensus du Run", "Séquences de Référence"),
  col = c("#D9534F", "#2B6CB0"),
  pch = 19,
  pt.cex = 1.2,
  bty = "n",
  cex = 0.85
)

dev.off()
cat("✅ Arbre rectangulaire PDF généré avec succès :", output_pdf, "\n")