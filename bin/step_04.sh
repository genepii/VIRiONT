#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS VIRiONT V2
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="fastq_pass"

RESULTS_DIR="results_test"
MERGED_DIR="${RESULTS_DIR}/01_MERGED"
DEHOST_DIR="${RESULTS_DIR}/02_DEHOSTING"
TRIMMED_DIR="${RESULTS_DIR}/03_FILTERED_TRIMMED"
GENOTYPING_DIR="${RESULTS_DIR}/05_GENOTYPING"
PRECONS_DIR="${RESULTS_DIR}/06_PRECONSENSUS"
QC_DIR="${RESULTS_DIR}/08_QC_ANALYSIS"

R_COMPUTE_METRICS="${SCRIPT_DIR}/compute_metrics.R"
SUMMARY_TABLE="${QC_DIR}/METRIC_summary_table.csv"

# --- RESSOURCES PARALLÈLES ---
PARALLEL_JOBS=4

# --- DÉTECTION AUTOMATIQUE DU SAMPLE SHEET DANS FASTQ_PASS ---
CSV_FILE=$(find "$DATA_DIR" -maxdepth 1 \( -name "*.csv" -o -name "*Sample_Sheet*" -o -name "*VHB*" -o -name "*VHD*" \) -type f | head -n 1 || true)

if [ -z "$CSV_FILE" ] || [ ! -f "$CSV_FILE" ]; then
    echo "❌ Erreur : Aucun fichier CSV / Sample Sheet trouvé dans $DATA_DIR/."
    exit 1
fi

echo "📄 Sample Sheet détecté : $CSV_FILE"

# Détection du séparateur (virgule ou point-virgule)
SEP=";"
if head -n 5 "$CSV_FILE" | grep -q ","; then
    SEP=","
fi

mkdir -p "$QC_DIR"

echo "====================================================================="
echo "   ÉTAPE 04 : MÉTROLOGIE QC COMPLÈTE (PARALLÉLISÉE AVEC XARGS : $PARALLEL_JOBS JOBS)"
echo "   Tableau récapitulatif final : $SUMMARY_TABLE"
echo "====================================================================="

# =====================================================================
# FONCTION UNIFIÉE DE TRAITEMENT EXPORTÉE POUR XARGS
# =====================================================================
process_single_sample_qc() {
    local sample_id="$1"
    local MERGED_DIR="$2"
    local DEHOST_DIR="$3"
    local TRIMMED_DIR="$4"
    local GENOTYPING_DIR="$5"
    local PRECONS_DIR="$6"
    local QC_DIR="$7"

    local sample_raw="${QC_DIR}/${sample_id}_raw.csv"
    local sample_dehost="${QC_DIR}/${sample_id}_dehost.csv"
    local sample_trimm="${QC_DIR}/${sample_id}_trimm.csv"
    local sample_geno="${QC_DIR}/${sample_id}_geno.csv"

    > "$sample_raw"
    > "$sample_dehost"
    > "$sample_trimm"
    > "$sample_geno"

    echo "---------------------------------------------------------------------"
    echo "▶️ [QC] Extraction des données : $sample_id"
    echo "---------------------------------------------------------------------"

    # --- 1. Étape 01_RAW ---
    local raw_fastq="${MERGED_DIR}/${sample_id}_merged.fastq.gz"
    if [ -f "$raw_fastq" ] && [ -s "$raw_fastq" ]; then
        gzip -dc "$raw_fastq" | awk -v sid="$sample_id" '(NR%4==2){print length($0)";"sid";01_RAW;NONE;NA;NA;NA;NA"}' >> "$sample_raw" || true
    fi

    # --- 2. Étape 02_DEHOSTING ---
    local dehost_fastq="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"
    [ ! -f "$dehost_fastq" ] && dehost_fastq="${DEHOST_DIR}/${sample_id}_meta.fastq.gz"
    if [ -f "$dehost_fastq" ] && [ -s "$dehost_fastq" ]; then
        gzip -dc "$dehost_fastq" | awk -v sid="$sample_id" '(NR%4==2){print length($0)";"sid";02_DEHOSTING;NONE;NA;NA;NA;NA"}' >> "$sample_dehost" || true
    fi

    # --- 3. Étape 03_FILTERED_TRIMMED ---
    local trim_fastq="${TRIMMED_DIR}/${sample_id}_trimmed.fastq.gz"
    if [ -f "$trim_fastq" ] && [ -s "$trim_fastq" ]; then
        gzip -dc "$trim_fastq" | awk -v sid="$sample_id" '(NR%4==2){print length($0)";"sid";03_FILTERED_TRIMMED;NONE;NA;NA;NA;NA"}' >> "$sample_trimm" || true
    fi

    # --- Métriques de 06_PRECONSENSUS ---
    local mean_depth="NA"
    local bam_file=$(find "${PRECONS_DIR}/BAM/${sample_id}" -name "*.bam" 2>/dev/null | head -n 1 || true)
    if [ -n "$bam_file" ] && [ -f "$bam_file" ]; then
        mean_depth=$(samtools depth "$bam_file" 2>/dev/null | awk '{sum+=$3; cnt++} END {if(cnt>0) printf "%.1f", sum/cnt; else print "NA"}' || echo "NA")
    fi

    local clair3_vars="NA"
    local vcf_file=$(find "${PRECONS_DIR}/VCF/${sample_id}" -name "*.vcf" 2>/dev/null | head -n 1 || true)
    if [ -n "$vcf_file" ] && [ -f "$vcf_file" ]; then
        clair3_vars=$(grep -v "^#" "$vcf_file" | grep -cw "PASS" || grep -v -c "^#" || echo "NA")
    fi

    # --- 4. Étape 05_GENOTYPING ---
    local tsv_file="${GENOTYPING_DIR}/${sample_id}/${sample_id}_validated_genotypes.tsv"
    if [ -f "$tsv_file" ]; then
        grep -w "VALIDATED" "$tsv_file" | while IFS=$'\t' read -r s_id geno best_cluster total_reads ratio pident length status; do
            # NETTOYAGE STRICT DES RETOURS À LA LIGNE POUR ÉVITER LA CORRUPTION DU FICHIER CSV
            geno=$(echo "$geno" | tr -d '\r\n ')
            pident=$(echo "$pident" | tr -d '\r\n ')
            total_reads=$(echo "$total_reads" | tr -d '\r\n ')

            if [ -f "$trim_fastq" ] && [ -s "$trim_fastq" ]; then
                gzip -dc "$trim_fastq" | awk -v sid="$sample_id" -v g="$geno" -v p="$pident" -v d="$mean_depth" -v v="$clair3_vars" -v a="$total_reads" \
                    '(NR%4==2){print length($0)";"sid";05_GENOTYPING;"g";"p";"d";"v";"a}' >> "$sample_geno" || true
            fi
        done
    fi

    echo "    ✅ [QC] Complété avec succès pour $sample_id"
}

export -f process_single_sample_qc

# =====================================================================
# EXTRACTION PARALLÈLE DE LA LISTE DES SAMPLES (XARGS -P 4)
# =====================================================================
SAMPLE_LIST=$(mktemp)

tail -n +2 "$CSV_FILE" | while IFS="$SEP" read -r col1 col2 col3 col4 col5 col6 col7 col8 || [ -n "$col1" ]; do
    sample=$(echo "$col2" | tr -d '\r\n ')
    component=$(echo "$col3" | tr -d '\r\n ')
    col1=$(echo "$col1" | tr -d '\r\n ')

    [ -z "$sample" ] && continue
    [[ "$col1" == *"Plate"* ]] && continue
    [[ "$sample" == *"Sample"* ]] && continue

    num_barcode=$(echo "$component" | grep -o '[0-9]\+' | head -n 1 || true)
    [ -z "$num_barcode" ] && continue

    formatted_barcode=$(printf "%02d" "$num_barcode")
    echo "barcode_${formatted_barcode}_${sample}" >> "$SAMPLE_LIST"
done

# Exécution parallèle avec xargs
cat "$SAMPLE_LIST" | xargs -I {} -P "$PARALLEL_JOBS" bash -c 'process_single_sample_qc "$@"' _ {} "$MERGED_DIR" "$DEHOST_DIR" "$TRIMMED_DIR" "$GENOTYPING_DIR" "$PRECONS_DIR" "$QC_DIR"

rm -f "$SAMPLE_LIST"

# =====================================================================
# CONCATÉNATION DES RÉSULTATS INTERMÉDIAIRES ET CALCUL R
# =====================================================================
ALL_RAW="${QC_DIR}/all_rawcount.csv"
ALL_DEHOST="${QC_DIR}/all_dehostcount.csv"
ALL_TRIMM="${QC_DIR}/all_trimcount.csv"
ALL_GENO="${QC_DIR}/all_genocount.csv"

rm -f "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO"

cat "${QC_DIR}"/*_raw.csv > "$ALL_RAW" 2>/dev/null || touch "$ALL_RAW"
cat "${QC_DIR}"/*_dehost.csv > "$ALL_DEHOST" 2>/dev/null || touch "$ALL_DEHOST"
cat "${QC_DIR}"/*_trimm.csv > "$ALL_TRIMM" 2>/dev/null || touch "$ALL_TRIMM"
cat "${QC_DIR}"/*_geno.csv > "$ALL_GENO" 2>/dev/null || touch "$ALL_GENO"

# Nettoyage des petits fichiers CSV temporaires
rm -f "${QC_DIR}"/*_raw.csv "${QC_DIR}"/*_dehost.csv "${QC_DIR}"/*_trimm.csv "${QC_DIR}"/*_geno.csv

# --- CALCULATION VIA R ---
echo "---------------------------------------------------------------------"
echo "--> [R] Calcul du tableau métrologique final..."

Rscript "$R_COMPUTE_METRICS" "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO" "$SUMMARY_TABLE"

rm -f "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO"

echo "====================================================================="
echo " 🎉 Étape 04 terminée avec succès en mode parallèle !"
echo " 📁 Fichier récapitulatif généré dans : $SUMMARY_TABLE"
echo "====================================================================="