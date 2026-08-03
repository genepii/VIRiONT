#!/bin/bash

# Interruption en cas d'erreur globale bloquante
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS VIRiONT V2
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="fastq_pass"

chmod +x "$SCRIPT_DIR"/*.sh "$SCRIPT_DIR"/*.R 2>/dev/null || true

RESULTS_DIR="results_test"
PRECONS_VCF_DIR="${RESULTS_DIR}/06_PRECONSENSUS/VCF"
MUTATION_DIR="${RESULTS_DIR}/11_MUTATION_SCREENING"
MUTATION_TABLES_DIR="${BASE_DIR}/mutation_table"

R_SEARCH_MUT="${SCRIPT_DIR}/search_mutation.R"

FREQ_MIN=5          # Seuil minimal de fréquence (%)
WINDOW_POS=0
PARALLEL_JOBS=4     # Traitement parallèle (4 VCFs en simultané)

# -----------------------------------------------------------------
# 1. DÉTECTION DU VIRUS VIA SAMPLE SHEET
# -----------------------------------------------------------------
CSV_FILE=$(find "$DATA_DIR" -maxdepth 1 \( -name "*.csv" -o -name "*Sample_Sheet*" -o -name "*VHB*" -o -name "*VHD*" \) -type f | head -n 1 || true)

virus_name="VHB"
if [ -n "$CSV_FILE" ] && [ -f "$CSV_FILE" ]; then
    csv_upper=$(echo "$CSV_FILE" | tr '[:lower:]' '[:upper:]')
    if [[ "$csv_upper" == *"VHD"* ]]; then virus_name="VHD"; fi
fi

echo "====================================================================="
echo "   ÉTAPE 11_MUTATION_SCREENING : RECHERCHE DES MUTATIONS CLINIQUES VHB"
echo "   Virus détecté      : $virus_name"
echo "   Sample Sheet       : ${CSV_FILE:-"Non trouvé (Défaut VHB)"}"
echo "   Dossier VCF source : $PRECONS_VCF_DIR"
echo "   Dossier des tables : $MUTATION_TABLES_DIR"
echo "   Dossier de sortie  : $MUTATION_DIR"
echo "   Jobs en parallèle  : $PARALLEL_JOBS"
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
# FONCTION UNIFIÉE DE SCREENING R (EXPORTÉE POUR XARGS)
# =====================================================================
process_single_vcf() {
    local vcf_path="$1"
    local MUTATION_DIR="$2"
    local MUTATION_TABLES_DIR="$3"
    local R_SEARCH_MUT="$4"
    local FREQ_MIN="$5"
    local WINDOW_POS="$6"

    [ ! -f "$vcf_path" ] && return 0

    local sample_dir="$(basename "$(dirname "$vcf_path")")"
    local filename="$(basename "$vcf_path")"
    local ref_name="${filename%.vcf}"

    echo "---------------------------------------------------------------------"
    echo "▶️ [MUTATION SCREENING] $sample_dir --> $filename"
    echo "---------------------------------------------------------------------"

    local sample_mut_dir="${MUTATION_DIR}/${sample_dir}"
    local all_res_dir="${sample_mut_dir}/all_results"
    local filt_res_dir="${sample_mut_dir}/filtered"

    mkdir -p "$all_res_dir" "$filt_res_dir"

    local vcf_copy="${sample_mut_dir}/${filename}"
    cp "$vcf_path" "$vcf_copy"

    # GARANTIE DU SAUT DE LIGNE FINAL (Supprime le warning R readLines)
    [ -s "$vcf_copy" ] && sed -i -e '$a\' "$vcf_copy" 2>/dev/null || true

    # Chemins des fichiers de sortie
    local f_VARIANTS="${sample_mut_dir}/${ref_name}_vcf_variants.csv"

    local f_PC="${all_res_dir}/${ref_name}_PreCore.csv"
    local f_BCP="${all_res_dir}/${ref_name}_BCP.csv"
    local f_DS="${all_res_dir}/${ref_name}_DomaineS.csv"
    local f_RT="${all_res_dir}/${ref_name}_DomaineRT.csv"
    local f_DPS1="${all_res_dir}/${ref_name}_DomainePreS1.csv"
    local f_DPS2="${all_res_dir}/${ref_name}_DomainePreS2.csv"
    local f_DHBx="${all_res_dir}/${ref_name}_DomaineHBx.csv"
    local f_C="${all_res_dir}/${ref_name}_Core.csv"

    local f_PC_F="${filt_res_dir}/${ref_name}_PreCore.csv"
    local f_BCP_F="${filt_res_dir}/${ref_name}_BCP.csv"
    local f_DS_F="${filt_res_dir}/${ref_name}_DomaineS.csv"
    local f_RT_F="${filt_res_dir}/${ref_name}_DomaineRT.csv"
    local f_DPS1_F="${filt_res_dir}/${ref_name}_DomainePreS1.csv"
    local f_DPS2_F="${filt_res_dir}/${ref_name}_DomainePreS2.csv"
    local f_DHBx_F="${filt_res_dir}/${ref_name}_DomaineHBx.csv"
    local f_C_F="${filt_res_dir}/${ref_name}_Core.csv"

    # Appel Rscript strict avec les 21 arguments
    Rscript "$R_SEARCH_MUT" \
        "$vcf_copy" \
        "$MUTATION_TABLES_DIR" \
        "$FREQ_MIN" \
        "$WINDOW_POS" \
        "$f_PC" "$f_BCP" "$f_DS" "$f_RT" "$f_DPS1" "$f_DPS2" "$f_DHBx" "$f_C" \
        "$f_PC_F" "$f_BCP_F" "$f_DS_F" "$f_RT_F" "$f_DPS1_F" "$f_DPS2_F" "$f_DHBx_F" "$f_C_F" \
        "$f_VARIANTS" > /dev/null 2>&1 || true

    echo "✅ [MUTATION SCREENING] Terminé pour $sample_dir ($ref_name)"
}

export -f process_single_vcf

# =====================================================================
# LISTAGE ET EXÉCUTION PARALLÈLE
# =====================================================================
VCF_LIST=$(mktemp)

find "$PRECONS_VCF_DIR" -type f -name "*.vcf" ! -name "*.idx" > "$VCF_LIST"

total_vcfs=$(wc -l < "$VCF_LIST" || echo 0)

if [ "$total_vcfs" -eq 0 ]; then
    echo "⚠️ Aucun fichier .vcf n'a été trouvé dans $PRECONS_VCF_DIR."
    rm -f "$VCF_LIST"
    exit 0
fi

# Traitement parallèle xargs
cat "$VCF_LIST" | xargs -I {} -P "$PARALLEL_JOBS" bash -c 'process_single_vcf "$@"' _ {} "$MUTATION_DIR" "$MUTATION_TABLES_DIR" "$R_SEARCH_MUT" "$FREQ_MIN" "$WINDOW_POS"

rm -f "$VCF_LIST"

echo "====================================================================="
echo " 🎉 Étape 11_MUTATION_SCREENING terminée avec succès !"
echo " 📊 Total : $total_vcfs fichier(s) VCF analysé(s) sous $MUTATION_DIR"
echo "====================================================================="