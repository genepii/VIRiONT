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

FINAL_CONSENSUS_DIR="${RESULTS_DIR}/07_CONSENSUS"
LOCAL_MODELS_DIR="${BASE_DIR}/models"

# --- CONFIGURATION DES RESSOURCES CPU PARALLÈLES ---
PARALLEL_JOBS=4
THREADS_PER_JOB=16

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
echo "   ÉTAPE 03 : DUAL-ENGINE (PARALLÉLISÉ AVEC XARGS : $PARALLEL_JOBS ÉCHANTILLONS EN SIMULTANÉ)"
echo "   Arborescence pré-consensus : $PRECONSENSUS_DIR"
echo "   Séquence consensus finale  : $FINAL_CONSENSUS_DIR"
echo "   Modèle Medaka utilisé      : $MEDAKA_MODEL"
echo "   Dossier modèle Clair3      : ${CLAIR3_MODEL_PATH:-"NON DISPONIBLE (Repli FreeBayes)"}"
echo "   Threads / Job              : $THREADS_PER_JOB"
echo "====================================================================="

# =====================================================================
# FONCTION UNIFIÉE DE TRAITEMENT EXPORTÉE POUR XARGS
# =====================================================================
process_single_tsv() {
    local tsv_file="$1"
    local TRIMMED_DIR="$2"
    local AMPLICON_DIR="$3"
    local PRECONS_SEQ_DIR="$4"
    local PRECONS_BAM_DIR="$5"
    local PRECONS_VCF_DIR="$6"
    local FINAL_CONSENSUS_DIR="$7"
    local MEDAKA_MODEL="$8"
    local CLAIR3_MODEL_PATH="$9"
    local THREADS_PER_JOB="${10}"

    [ ! -f "$tsv_file" ] && return 0

    grep -w "VALIDATED" "$tsv_file" | while IFS=$'\t' read -r sample_id geno best_cluster total_reads ratio pident length status; do

        sample_id=$(echo "$sample_id" | tr -d '\r\n ')
        geno=$(echo "$geno" | tr -d '\r\n ')
        best_cluster=$(echo "$best_cluster" | tr -d '\r\n ')

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
        echo "▶️ [DEBUT] Traitement : $sample_geno_id (Cluster source : $best_cluster)"
        echo "---------------------------------------------------------------------"

        if [ ! -f "$TRIMMED_FASTQ" ]; then
            echo "⚠️ FASTQ filtré introuvable pour $sample_id dans $TRIMMED_DIR. Sauté."
            continue
        fi

        # --- 1. EXTRACTION BLINDÉE ET SECOURUE DE LA SÉQUENCE DRAFT ---
        > "$PRECONS_FASTA"
        
        # Secours 1 : Fichier Master FASTA créé par l'Étape 02
        sample_dir="$(dirname "$tsv_file")"
        master_fasta="${sample_dir}/${sample_id}_${geno}_master_consensus.fasta"
        if [ -f "$master_fasta" ] && [ -s "$master_fasta" ]; then
            awk -v new_id=">$sample_geno_id" '
                NR==1 {print new_id; next}
                /^>/ {exit}
                {print}
            ' "$master_fasta" > "$PRECONS_FASTA"
        fi

        # Secours 2 : Recherche dans le dossier amplicon_sorter par nom de fichier
        if [ ! -s "$PRECONS_FASTA" ]; then
            cluster_file=$(find "${AMPLICON_DIR}/${sample_id}" -type f \( -name "*${best_cluster}*.fasta" -o -name "*${best_cluster}*.fa" \) ! -name "*unique*" 2>/dev/null | head -n 1 || true)
            if [ -n "$cluster_file" ] && [ -s "$cluster_file" ]; then
                awk -v new_id=">$sample_geno_id" '
                    NR==1 {print new_id; next}
                    /^>/ {exit}
                    {print}
                ' "$cluster_file" > "$PRECONS_FASTA"
            fi
        fi

        # Secours 3 : Prendre n'importe quel fichier FASTA dans le dossier amplicon_sorter
        if [ ! -s "$PRECONS_FASTA" ]; then
            any_fasta=$(find "${AMPLICON_DIR}/${sample_id}" -type f \( -name "*.fasta" -o -name "*.fa" \) ! -name "*unique*" 2>/dev/null | head -n 1 || true)
            if [ -n "$any_fasta" ] && [ -s "$any_fasta" ]; then
                awk -v new_id=">$sample_geno_id" '
                    NR==1 {print new_id; next}
                    /^>/ {exit}
                    {print}
                ' "$any_fasta" > "$PRECONS_FASTA"
            fi
        fi

        if [ ! -s "$PRECONS_FASTA" ]; then
            echo "❌ Erreur : Impossible de récupérer la séquence draft pour $sample_id ($geno). Sauté."
            continue
        fi

        echo "--> [06_PRECONSENSUS/SEQUENCES/${sample_id}] Patron pré-consensus prêt : $(basename "$PRECONS_FASTA")"

        # --- 2. Alignement BAM ---
        echo "--> [06_PRECONSENSUS/BAM/${sample_id}] Alignement Minimap2 et indexation BAM..."
        minimap2 -ax map-ont --secondary=no -L -t "$THREADS_PER_JOB" "$PRECONS_FASTA" "$TRIMMED_FASTQ" 2>/dev/null | \
        samtools sort -@ "$THREADS_PER_JOB" -o "$OUT_BAM"
        samtools index "$OUT_BAM"

        # --- 3. Polissage Medaka ---
        echo "--> [MEDAKA] Polissage de la séquence FASTA pour $sample_geno_id..."
        rm -rf "${TEMP_DIR}_medaka"

        medaka_consensus \
            -i "$TRIMMED_FASTQ" \
            -d "$PRECONS_FASTA" \
            -o "${TEMP_DIR}_medaka" \
            -m "$MEDAKA_MODEL" \
            -t "$THREADS_PER_JOB" > /dev/null 2>&1 || true

        if [ -f "${TEMP_DIR}_medaka/consensus.fasta" ] && [ -s "${TEMP_DIR}_medaka/consensus.fasta" ]; then
            cp "${TEMP_DIR}_medaka/consensus.fasta" "$FINAL_FASTA"
            rm -rf "${TEMP_DIR}_medaka"
        else
            echo "⚠️ Polissage Medaka non concluant pour $sample_geno_id. Pré-consensus conservé en secours."
            cp "$PRECONS_FASTA" "$FINAL_FASTA"
            rm -rf "${TEMP_DIR}_medaka"
        fi

        # --- 4. Variant Calling Clair3 ---
        echo "--> [CLAIR3] Génération du VCF des variants pour $sample_geno_id..."
        rm -rf "${TEMP_DIR}_clair3"

        if [ -n "$CLAIR3_MODEL_PATH" ] && [ -d "$CLAIR3_MODEL_PATH" ]; then
            run_clair3.sh \
                --bam_fn="$OUT_BAM" \
                --ref_fn="$PRECONS_FASTA" \
                --threads="$THREADS_PER_JOB" \
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

        echo "✅ [FIN] Traitement terminé pour $sample_geno_id"

    done
}

export -f process_single_tsv

# =====================================================================
# LISTAGE DES FICHIERS TSV ET PARALLÉLISATEUR XARGS
# =====================================================================
TSV_LIST=$(mktemp)

find "$GENOTYPING_DIR" -maxdepth 2 -name "*_validated_genotypes.tsv" > "$TSV_LIST"

if [ ! -s "$TSV_LIST" ]; then
    echo "⚠️ Aucun fichier de génotypage n'a été trouvé dans $GENOTYPING_DIR."
    rm -f "$TSV_LIST"
    exit 0
fi

# Lancement parallèle sur les fichiers TSV
cat "$TSV_LIST" | xargs -I {} -P "$PARALLEL_JOBS" bash -c 'process_single_tsv "$@"' _ {} "$TRIMMED_DIR" "$AMPLICON_DIR" "$PRECONS_SEQ_DIR" "$PRECONS_BAM_DIR" "$PRECONS_VCF_DIR" "$FINAL_CONSENSUS_DIR" "$MEDAKA_MODEL" "$CLAIR3_MODEL_PATH" "$THREADS_PER_JOB"

rm -f "$TSV_LIST"

# =====================================================================
# CONCATÉNATION GLOBALE DE TOUS LES CONSENSUS DANS 07_CONSENSUS
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
echo " 🎉 Étape 03 terminée avec succès en mode parallèle !"
echo "====================================================================="