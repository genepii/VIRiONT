#!/bin/bash

# Interruption en cas d'erreur globale bloquante
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="results_test"

TRIMMED_DIR="${RESULTS_DIR}/03_FILTERED_TRIMMED"
GENOTYPING_DIR="${RESULTS_DIR}/05_GENOTYPING"
AMPLICON_DIR="${RESULTS_DIR}/04_AMPLICON_SORTER"

PRECONSENSUS_DIR="${RESULTS_DIR}/06_PRECONSENSUS"
PRECONS_BAM_DIR="${PRECONSENSUS_DIR}/BAM"
PRECONS_SEQ_DIR="${PRECONSENSUS_DIR}/SEQUENCES"
PRECONS_VCF_DIR="${PRECONSENSUS_DIR}/VCF"

# Renommé de 09_CONSENSUS -> 07_CONSENSUS pour suivre l'ordre logique
FINAL_CONSENSUS_DIR="${RESULTS_DIR}/07_CONSENSUS"
LOCAL_MODELS_DIR="${BASE_DIR}/models"

# --- CONFIGURATION DES MODÈLES D'IA ---
MEDAKA_MODEL="r1041_e82_400bps_sup_g615"
CLAIR3_MODEL_NAME="r1041_e82_400bps_sup"

# =====================================================================
# DÉTECTION DU MODÈLE CLAIR3 (EMBARQUÉ OU LOCAL)
# =====================================================================
SEARCH_PATHS=(
    "/opt/models"
    "$LOCAL_MODELS_DIR"
    "${CONDA_PREFIX}/bin/models"
)

CLAIR3_MODEL_PATH=""
for path in "${SEARCH_PATHS[@]}"; do
    if [ -d "$path" ]; then
        # Exclut les modèles avec move-table (*_with_mv) incompatibles
        found=$(find "$path" -maxdepth 2 -type d -name "${CLAIR3_MODEL_NAME}*" ! -name "*_with_mv" 2>/dev/null | head -n 1 || true)
        if [ -n "$found" ] && [ -d "$found" ]; then
            CLAIR3_MODEL_PATH="$found"
            break
        fi
    fi
done

# Création des dossiers principaux de sortie
mkdir -p "$PRECONS_BAM_DIR" "$PRECONS_SEQ_DIR" "$PRECONS_VCF_DIR" "$FINAL_CONSENSUS_DIR"

echo "====================================================================="
echo "  ÉTAPE 03 : DUAL-ENGINE (Polissage Medaka & Variant Calling Clair3)"
echo "  Arborescence pré-consensus : $PRECONSENSUS_DIR"
echo "  Séquence consensus finale  : $FINAL_CONSENSUS_DIR"
echo "  Modèle Medaka utilisé      : $MEDAKA_MODEL"
echo "  Dossier modèle Clair3      : ${CLAIR3_MODEL_PATH:-"NON DISPONIBLE (Repli FreeBayes)"}"
echo "====================================================================="

# Boucle de traitement récursive sur les sous-dossiers de 05_GENOTYPING
find "$GENOTYPING_DIR" -name "*_validated_genotypes.tsv" | while read -r tsv_file; do

    grep -w "VALIDATED" "$tsv_file" | while IFS=$'\t' read -r sample_id geno best_cluster total_reads ratio pident length status; do

        [ -z "$geno" ] && continue

        # --- CRÉATION DES SOUS-DOSSIERS PAR ÉCHANTILLON ---
        SAMPLE_SEQ_DIR="${PRECONS_SEQ_DIR}/${sample_id}"
        SAMPLE_BAM_DIR="${PRECONS_BAM_DIR}/${sample_id}"
        SAMPLE_VCF_DIR="${PRECONS_VCF_DIR}/${sample_id}"

        mkdir -p "$SAMPLE_SEQ_DIR" "$SAMPLE_BAM_DIR" "$SAMPLE_VCF_DIR"

        sample_geno_id="${sample_id}_${geno}"
        TRIMMED_FASTQ="${TRIMMED_DIR}/${sample_id}_trimmed.fastq.gz"
        
        PRECONS_FASTA="${SAMPLE_SEQ_DIR}/${sample_geno_id}_preconsensus.fasta"
        OUT_BAM="${SAMPLE_BAM_DIR}/${sample_geno_id}.bam"
        OUT_VCF="${SAMPLE_VCF_DIR}/${sample_geno_id}.vcf"
        FINAL_FASTA="${FINAL_CONSENSUS_DIR}/${sample_geno_id}_consensus.fasta"
        
        TEMP_DIR="${FINAL_CONSENSUS_DIR}/temp_${sample_geno_id}"

        echo "---------------------------------------------------------------------"
        echo "Traitement : $sample_geno_id (Cluster source : $best_cluster)"
        echo "---------------------------------------------------------------------"

        if [ ! -f "$TRIMMED_FASTQ" ]; then
            echo "⚠️ FASTQ filtré introuvable pour $sample_id. Sauté."
            continue
        fi

        # --- 1. Extraction Séquence Draft ---
        AMPLICON_FASTA=$(find "${AMPLICON_DIR}/${sample_id}" -name "*consensusfile.fasta" 2>/dev/null | head -n 1 || true)

        if [ -n "$AMPLICON_FASTA" ] && [ -f "$AMPLICON_FASTA" ]; then
            awk -v target=">$best_cluster" -v new_id=">$sample_geno_id" '
                $0 ~ target {flag=1; print new_id; next}
                /^>/ {flag=0}
                flag {print}
            ' "$AMPLICON_FASTA" > "$PRECONS_FASTA"
        fi

        # Recherche du master consensus dans le sous-dossier dédié de l'échantillon
        if [ ! -s "$PRECONS_FASTA" ]; then
            sample_dir="$(dirname "$tsv_file")"
            master_fasta="${sample_dir}/${sample_id}_${geno}_master_consensus.fasta"
            if [ -f "$master_fasta" ]; then
                awk -v new_id=">$sample_geno_id" '
                    NR==1 {print new_id; next}
                    /^>/ {exit}
                    {print}
                ' "$master_fasta" > "$PRECONS_FASTA"
            fi
        fi

        if [ ! -s "$PRECONS_FASTA" ]; then
            echo "⚠️ Impossible d'extraire la séquence unique pour $best_cluster."
            continue
        fi

        echo "--> [06_PRECONSENSUS/SEQUENCES/${sample_id}] Patron pré-consensus prêt : $(basename "$PRECONS_FASTA")"

        # --- 2. Alignement BAM ---
        echo "--> [06_PRECONSENSUS/BAM/${sample_id}] Alignement Minimap2 et indexation BAM..."
        minimap2 -ax map-ont --secondary=no -L -t 4 "$PRECONS_FASTA" "$TRIMMED_FASTQ" 2>/dev/null | \
        samtools sort -o "$OUT_BAM"
        samtools index "$OUT_BAM"

        # --- 3. Polissage Medaka ---
        echo "--> [MEDAKA] Polissage de la séquence FASTA..."
        rm -rf "${TEMP_DIR}_medaka"

        medaka_consensus \
            -i "$TRIMMED_FASTQ" \
            -d "$PRECONS_FASTA" \
            -o "${TEMP_DIR}_medaka" \
            -m "$MEDAKA_MODEL" \
            -t 4 > /dev/null 2>&1 || true

        if [ -f "${TEMP_DIR}_medaka/consensus.fasta" ] && [ -s "${TEMP_DIR}_medaka/consensus.fasta" ]; then
            cp "${TEMP_DIR}_medaka/consensus.fasta" "$FINAL_FASTA"
            rm -rf "${TEMP_DIR}_medaka"
        else
            echo "⚠️ Polissage Medaka non concluant pour $sample_geno_id. Pré-consensus conservé en secours."
            cp "$PRECONS_FASTA" "$FINAL_FASTA"
            rm -rf "${TEMP_DIR}_medaka"
        fi

        # --- 4. Variant Calling Clair3 ---
        echo "--> [CLAIR3] Génération du VCF des variants..."
        rm -rf "${TEMP_DIR}_clair3"

        if [ -n "$CLAIR3_MODEL_PATH" ] && [ -d "$CLAIR3_MODEL_PATH" ]; then
            run_clair3.sh \
                --bam_fn="$OUT_BAM" \
                --ref_fn="$PRECONS_FASTA" \
                --threads=4 \
                --platform="ont" \
                --model_path="$CLAIR3_MODEL_PATH" \
                --output="${TEMP_DIR}_clair3" \
                --include_all_ctgs \
                --haploid_precise \
                --no_phasing_for_fa \
                --chunk_size=5000 > /dev/null 2>&1 || true

            if [ -f "${TEMP_DIR}_clair3/merge_output.vcf.gz" ]; then
                gunzip -c "${TEMP_DIR}_clair3/merge_output.vcf.gz" > "$OUT_VCF"
            fi
            rm -rf "${TEMP_DIR}_clair3"
        else
            echo "⚠️ Modèle Clair3 non disponible. Repli sur FreeBayes..."
            freebayes -f "$PRECONS_FASTA" "$OUT_BAM" > "$OUT_VCF" 2>/dev/null || true
        fi

        if [ ! -s "$OUT_VCF" ]; then
            echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$OUT_VCF"
        fi

        echo "   ✅ SEQUENCES : ${sample_id}/$(basename "$PRECONS_FASTA")"
        echo "   ✅ BAM       : ${sample_id}/$(basename "$OUT_BAM")"
        echo "   ✅ VCF       : ${sample_id}/$(basename "$OUT_VCF")"
        echo "   ✅ CONSENSUS : $(basename "$FINAL_FASTA")"

    done
done

# =====================================================================
# 5. CONCATÉNATION GLOBALE DE TOUS LES CONSENSUS DANS 07_CONSENSUS
# =====================================================================
ALL_CONSENSUS_FILE="${FINAL_CONSENSUS_DIR}/all_samples_consensus.fasta"

echo "---------------------------------------------------------------------"
echo "--> [07_CONSENSUS] Concaténation globale des consensus..."

rm -f "$ALL_CONSENSUS_FILE"
find "$FINAL_CONSENSUS_DIR" -maxdepth 1 -name "*_consensus.fasta" ! -name "all_samples_consensus.fasta" -exec cat {} + > "$ALL_CONSENSUS_FILE" 2>/dev/null || true

if [ -s "$ALL_CONSENSUS_FILE" ]; then
    echo "✅ Fichier multi-FASTA concaténé créé : $(basename "$ALL_CONSENSUS_FILE")"
fi

echo "====================================================================="
echo " 🎉 Étape 03 terminée avec succès !"
echo "====================================================================="