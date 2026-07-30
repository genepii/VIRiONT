#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS VIRiONT V2
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CSV_FILE="fastq_pass/260630_VHB_WG_ABS.csv"

RESULTS_DIR="results_test"
MERGED_DIR="${RESULTS_DIR}/01_MERGED"
DEHOST_DIR="${RESULTS_DIR}/02_DEHOSTING"
TRIMMED_DIR="${RESULTS_DIR}/03_FILTERED_TRIMMED"
GENOTYPING_DIR="${RESULTS_DIR}/05_GENOTYPING"
PRECONS_DIR="${RESULTS_DIR}/06_PRECONSENSUS"
QC_DIR="${RESULTS_DIR}/08_QC_ANALYSIS"

R_COMPUTE_METRICS="${SCRIPT_DIR}/compute_metrics.R"

ALL_RAW="${QC_DIR}/all_rawcount.csv"
ALL_DEHOST="${QC_DIR}/all_dehostcount.csv"
ALL_TRIMM="${QC_DIR}/all_trimcount.csv"
ALL_GENO="${QC_DIR}/all_genocount.csv"

SUMMARY_TABLE="${QC_DIR}/METRIC_summary_table.csv"

mkdir -p "$QC_DIR"

if [ ! -f "$CSV_FILE" ]; then
    echo "❌ Erreur : Fichier CSV introuvable ($CSV_FILE)."
    exit 1
fi

echo "====================================================================="
echo "  ÉTAPE 04 : MÉTROLOGIE QC COMPLÈTE (4 ÉTAPES STRATÉGIQUES)"
echo "  Tableau récapitulatif final : $SUMMARY_TABLE"
echo "====================================================================="

rm -f "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO"
touch "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO"

# =====================================================================
# EXTRACTION DES 4 ÉTAPES CLÉS
# =====================================================================
tail -n +2 "$CSV_FILE" | while IFS=';' read -r plate_pos sample component forward reverse flowcell kit_id primers; do

    [ -z "$sample" ] && continue

    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')

    num_barcode=$(echo "$component" | sed 's/[^0-9]//g')
    formatted_barcode=$(printf "%02d" "$num_barcode")
    sample_id="barcode_${formatted_barcode}_${sample}"

    echo "---------------------------------------------------------------------"
    echo "Extraction des données : $sample_id"
    echo "---------------------------------------------------------------------"

    # --- 1. Étape 01_RAW ---
    raw_fastq="${MERGED_DIR}/${sample_id}_merged.fastq.gz"
    if [ -f "$raw_fastq" ] && [ -s "$raw_fastq" ]; then
        zcat "$raw_fastq" | awk -v sid="$sample_id" '(NR%4==2){print length($0)";"sid";01_RAW;NONE;NA;NA;NA;NA"}' >> "$ALL_RAW"
    fi

    # --- 2. Étape 02_DEHOSTING ---
    dehost_fastq="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"
    [ ! -f "$dehost_fastq" ] && dehost_fastq="${DEHOST_DIR}/${sample_id}_meta.fastq.gz"
    if [ -f "$dehost_fastq" ] && [ -s "$dehost_fastq" ]; then
        zcat "$dehost_fastq" | awk -v sid="$sample_id" '(NR%4==2){print length($0)";"sid";02_DEHOSTING;NONE;NA;NA;NA;NA"}' >> "$ALL_DEHOST"
    fi

    # --- 3. Étape 03_FILTERED_TRIMMED ---
    trim_fastq="${TRIMMED_DIR}/${sample_id}_trimmed.fastq.gz"
    if [ -f "$trim_fastq" ] && [ -s "$trim_fastq" ]; then
        zcat "$trim_fastq" | awk -v sid="$sample_id" '(NR%4==2){print length($0)";"sid";03_FILTERED_TRIMMED;NONE;NA;NA;NA;NA"}' >> "$ALL_TRIMM"
    fi

    # --- Métriques de 06_PRECONSENSUS ---
    mean_depth="NA"
    bam_file=$(find "${PRECONS_DIR}/BAM/${sample_id}" -name "*.bam" 2>/dev/null | head -n 1 || true)
    if [ -n "$bam_file" ] && [ -f "$bam_file" ]; then
        mean_depth=$(samtools depth "$bam_file" 2>/dev/null | awk '{sum+=$3; cnt++} END {if(cnt>0) printf "%.1f", sum/cnt; else print "NA"}')
    fi

    clair3_vars="NA"
    vcf_file=$(find "${PRECONS_DIR}/VCF/${sample_id}" -name "*.vcf" 2>/dev/null | head -n 1 || true)
    if [ -n "$vcf_file" ] && [ -f "$vcf_file" ]; then
        clair3_vars=$(grep -v "^#" "$vcf_file" | grep -cw "PASS" || grep -v -c "^#")
    fi

    # --- 4. Étape 05_GENOTYPING ---
    tsv_file="${GENOTYPING_DIR}/${sample_id}/${sample_id}_validated_genotypes.tsv"
    if [ -f "$tsv_file" ]; then
        grep -w "VALIDATED" "$tsv_file" | while IFS=$'\t' read -r s_id geno best_cluster total_reads ratio pident length status; do
            if [ -f "$trim_fastq" ] && [ -s "$trim_fastq" ]; then
                zcat "$trim_fastq" | awk -v sid="$sample_id" -v g="$geno" -v p="$pident" -v d="$mean_depth" -v v="$clair3_vars" -v a="$total_reads" \
                    '(NR%4==2){print length($0)";"sid";05_GENOTYPING;"g";"p";"d";"v";"a}' >> "$ALL_GENO"
            fi
        done
    fi

    echo "   ✅ Métriques extraites avec succès."

done

# =====================================================================
# CALCULATION VIA R
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> [R] Calcul du tableau final..."

Rscript "$R_COMPUTE_METRICS" "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO" "$SUMMARY_TABLE"

rm -f "$ALL_RAW" "$ALL_DEHOST" "$ALL_TRIMM" "$ALL_GENO"

echo "====================================================================="
echo " 🎉 Étape 04 terminée avec succès !"
echo " 📁 Fichier récapitulatif généré dans : $SUMMARY_TABLE"
echo "====================================================================="