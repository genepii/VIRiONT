#!/bin/bash

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS VIRiONT V2
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

chmod +x "$SCRIPT_DIR"/*.sh "$SCRIPT_DIR"/*.R 2>/dev/null || true

CSV_FILE="fastq_pass/260630_VHB_WG_ABS.csv"
RESULTS_DIR="results_test"
PRECONS_VCF_DIR="${RESULTS_DIR}/06_PRECONSENSUS/VCF"
MUTATION_DIR="${RESULTS_DIR}/11_MUTATION_SCREENING"
MUTATION_TABLES_DIR="${BASE_DIR}/mutation_table"

R_SEARCH_MUT="${SCRIPT_DIR}/search_mutation.R"

FREQ_MIN=5          # Seuil minimal de fréquence (%)
WINDOW_POS=0

# -----------------------------------------------------------------
# 1. DÉTECTION DU VIRUS
# -----------------------------------------------------------------
virus_name="VHB"
if [ -f "$CSV_FILE" ]; then
    csv_upper=$(echo "$CSV_FILE" | tr '[:lower:]' '[:upper:]')
    if [[ "$csv_upper" == *"VHD"* ]]; then virus_name="VHD"; fi
fi

echo "====================================================================="
echo "  ÉTAPES 11_MUTATION_SCREENING : RECHERCHE DES MUTATIONS CLINIQUES VHB"
echo "  Virus détecté      : $virus_name"
echo "  Dossier VCF source : $PRECONS_VCF_DIR"
echo "  Dossier des tables : $MUTATION_TABLES_DIR"
echo "  Dossier de sortie  : $MUTATION_DIR"
echo "====================================================================="

if [ "$virus_name" == "VHD" ]; then
    echo "ℹ️  L'échantillon analysé est du VHD (HDV)."
    echo "⏭️  Le screening des mutations cliniques est réservé au VHB. Étape ignorée."
    exit 0
fi

if [ ! -d "$PRECONS_VCF_DIR" ]; then
    echo "⚠️ Dossier VCF introuvable ($PRECONS_VCF_DIR). Exécutez step_03.sh d'abord."
    exit 0
fi

if [ ! -d "$MUTATION_TABLES_DIR" ]; then
    echo "❌ Erreur : Dossier des tables de mutations VHB introuvable ($MUTATION_TABLES_DIR)."
    exit 1
fi

mkdir -p "$MUTATION_DIR"

# =====================================================================
# 2. PARCOURS SYSTEMATIQUE DE TOUS LES BARCODES DE 06_PRECONSENSUS/VCF
# =====================================================================
total_vcfs=0
total_barcodes=0

for sample_dir in "$PRECONS_VCF_DIR"/*; do
    
    [ ! -d "$sample_dir" ] && continue
    
    sample_id=$(basename "$sample_dir")
    ((total_barcodes++))

    echo "---------------------------------------------------------------------"
    echo "📂 [$total_barcodes] Screening des mutations pour : $sample_id"
    echo "---------------------------------------------------------------------"

    vcf_files=$(find "$sample_dir" -maxdepth 1 -name "*.vcf" ! -name "*.idx" 2>/dev/null || true)

    if [ -z "$vcf_files" ]; then
        echo "   ⚠️ Aucun fichier .vcf dans $sample_id. Sauté."
        continue
    fi

    for vcf in $vcf_files; do
        filename=$(basename "$vcf")
        ref_name="${filename%.vcf}"

        echo "   --> 🧬 Analyse du VCF : $filename"

        sample_mut_dir="${MUTATION_DIR}/${sample_id}"
        all_res_dir="${sample_mut_dir}/all_results"
        filt_res_dir="${sample_mut_dir}/filtered"

        mkdir -p "$all_res_dir" "$filt_res_dir"

        vcf_copy="${sample_mut_dir}/${filename}"
        cp "$vcf" "$vcf_copy"

        # Fichier d'annotation globale (ex: GTD_vcf_variants) déposé à la racine de l'échantillon
        f_VARIANTS="${sample_mut_dir}/${ref_name}_vcf_variants.csv"

        # Chemins des résultats bruts
        f_PC="${all_res_dir}/${ref_name}_PreCore.csv"
        f_BCP="${all_res_dir}/${ref_name}_BCP.csv"
        f_DS="${all_res_dir}/${ref_name}_DomaineS.csv"
        f_RT="${all_res_dir}/${ref_name}_DomaineRT.csv"
        f_DPS1="${all_res_dir}/${ref_name}_DomainePreS1.csv"
        f_DPS2="${all_res_dir}/${ref_name}_DomainePreS2.csv"
        f_DHBx="${all_res_dir}/${ref_name}_DomaineHBx.csv"
        f_C="${all_res_dir}/${ref_name}_Core.csv"

        # Chemins des résultats filtrés
        f_PC_F="${filt_res_dir}/${ref_name}_PreCore.csv"
        f_BCP_F="${filt_res_dir}/${ref_name}_BCP.csv"
        f_DS_F="${filt_res_dir}/${ref_name}_DomaineS.csv"
        f_RT_F="${filt_res_dir}/${ref_name}_DomaineRT.csv"
        f_DPS1_F="${filt_res_dir}/${ref_name}_DomainePreS1.csv"
        f_DPS2_F="${filt_res_dir}/${ref_name}_DomainePreS2.csv"
        f_DHBx_F="${filt_res_dir}/${ref_name}_DomaineHBx.csv"
        f_C_F="${filt_res_dir}/${ref_name}_Core.csv"

        # Exécution du script R (avec ajout de f_VARIANTS comme 21ème argument)
        Rscript "$R_SEARCH_MUT" \
            "$vcf_copy" \
            "$MUTATION_TABLES_DIR" \
            "$FREQ_MIN" \
            "$WINDOW_POS" \
            "$f_PC" "$f_BCP" "$f_DS" "$f_RT" "$f_DPS1" "$f_DPS2" "$f_DHBx" "$f_C" \
            "$f_PC_F" "$f_BCP_F" "$f_DS_F" "$f_RT_F" "$f_DPS1_F" "$f_DPS2_F" "$f_DHBx_F" "$f_C_F" \
            "$f_VARIANTS" || echo "⚠️ Erreur mineure Rscript sur $filename (Poursuite de la boucle)"

        ((total_vcfs++))
    done

done

echo "====================================================================="
echo " 🎉 Étape 11_MUTATION_SCREENING terminée avec succès !"
echo " 📊 Total : $total_barcodes barcode(s) et $total_vcfs fichier(s) VCF analysé(s) sous $MUTATION_DIR"
echo "====================================================================="