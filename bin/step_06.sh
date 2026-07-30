#!/bin/bash

# Ne pas couper brutalement sur les erreurs mineures de boucle
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

chmod +x "$SCRIPT_DIR"/*.sh "$SCRIPT_DIR"/*.R 2>/dev/null || true

RESULTS_DIR="results_test"
COVERAGE_DIR="${RESULTS_DIR}/10_COVERAGE"
mkdir -p "$COVERAGE_DIR"

# Détection du dossier BAM source
if [ -d "${RESULTS_DIR}/07_BAM" ] && [ "$(find "${RESULTS_DIR}/07_BAM" -name "*.bam" 2>/dev/null | head -n 1)" ]; then
    BAM_DIR="${RESULTS_DIR}/07_BAM"
elif [ -d "${RESULTS_DIR}/06_PRECONSENSUS/BAM" ]; then
    BAM_DIR="${RESULTS_DIR}/06_PRECONSENSUS/BAM"
else
    echo "❌ Erreur : Dossier BAM introuvable."
    exit 1
fi

COV_SUM="${COVERAGE_DIR}/cov_sum.cov"
COV_PLOT_PDF="${COVERAGE_DIR}/cov_plot.pdf"
R_PLOT_COV="${SCRIPT_DIR}/plot_cov_MI.R"

echo "====================================================================="
echo "  ÉTAPES 10_COVERAGE : GÉNÉRATION DU COV_PLOT.PDF"
echo "  Dossier BAM source : $BAM_DIR"
echo "  Dossier de sortie  : $COVERAGE_DIR"
echo "====================================================================="

# Reset du fichier global de couverture
> "$COV_SUM"

echo "---------------------------------------------------------------------"
echo "--> 1. Recherche et calcul de la couverture pour TOUS les échantillons..."

count=0

# Recherche robuste de TOUS les fichiers .bam dans tous les sous-dossiers
while IFS= read -r bam; do
    # Ignorer les index .bai s'ils remontent
    if [[ "$bam" == *.bai ]]; then
        continue
    fi
    
    sample_id=$(basename "$(dirname "$bam")")
    filename=$(basename "$bam")
    ref_name=$(echo "$filename" | sed 's/\.bam$//; s/_sorted$//; s/_ampliconclip$//')
    
    sample_cov_dir="${COVERAGE_DIR}/${sample_id}"
    mkdir -p "$sample_cov_dir"
    sample_cov_file="${sample_cov_dir}/${ref_name}.cov"
    
    echo "   📊 Traitement du fichier : $filename (Barcode : $sample_id)"

    # Génération de la couverture avec bedtools
    bedtools genomecov -ibam "$bam" -d -split | sed "s/$/\t${sample_id}\tMINION/" > "$sample_cov_file"
    
    cat "$sample_cov_file" >> "$COV_SUM"
    ((count++))
done < <(find "$BAM_DIR" -type f -name "*.bam")

echo "   ✅ $count profils de couverture extraits dans $COV_SUM."

if [ "$count" -eq 0 ]; then
    echo "❌ Erreur : Aucun BAM valide n'a été traité."
    exit 1
fi

# =====================================================================
# 2. GÉNÉRATION DU PDF (Rscript)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 2. Génération du graphique PDF avec Rscript..."

if [ ! -f "$R_PLOT_COV" ]; then
    echo "❌ Erreur : Script R introuvable à l'emplacement $R_PLOT_COV"
    exit 1
fi

# Exécution explicite de R
Rscript "$R_PLOT_COV" "$COV_SUM" "$COV_PLOT_PDF"

echo "====================================================================="
echo " 🎉 Étape 10_COVERAGE terminée avec succès !"
echo " 📄 Fichier résumé : $COV_SUM"
echo " 📄 Fichier PDF    : $COV_PLOT_PDF"
echo "====================================================================="