#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS
# =====================================================================
RESULTS_DIR="results_test"
REF_GLOBAL="${RESULTS_DIR}/00_SUPDATA/REFSEQ/HBV_subtype.fasta"
REFILTERED_DIR="${RESULTS_DIR}/05_REFILTERED_FASTQ"

MEDAKA_DIR="${RESULTS_DIR}/06_MEDAKA_CONSENSUS"
BAM_DIR="${RESULTS_DIR}/07_BAM"
VCF_DIR="${RESULTS_DIR}/08_VCF"
CONSENSUS_DIR="${RESULTS_DIR}/09_CONSENSUS"

# Paramètres Medaka (Chimie R10.4.1 / SUP)
MEDAKA_MODEL="r1041_e82_400bps_sup_g615"
THREADS=4

# Paramètres Samtools Mpileup
MPILEUP_DEPTH=10000
MPILEUP_QUAL=0

# Création des dossiers de sortie
mkdir -p "$MEDAKA_DIR" "$BAM_DIR" "$VCF_DIR" "$CONSENSUS_DIR"

echo "====================================================================="
echo "  EXECUTION ÉTAPE 03 (Consensus Hybride : Medaka + Samtools Mpileup)"
echo "  Référence globale : $REF_GLOBAL"
echo "  Modèle Medaka     : $MEDAKA_MODEL"
echo "====================================================================="

# Indexation de la référence globale si nécessaire
if [ -f "$REF_GLOBAL" ] && [ ! -f "${REF_GLOBAL}.fai" ]; then
    samtools faidx "$REF_GLOBAL"
fi

# =====================================================================
# BOUCLE SUR TOUS LES ÉCHANTILLONS ET LEURS GÉNOTYPES
# =====================================================================
for sample_path in "${REFILTERED_DIR}"/barcode_*; do
    [ -d "$sample_path" ] || continue
    
    sample_id=$(basename "$sample_path")

    # Parcours des fichiers FASTQ filtrés générés à l'étape 02
    for fastq_file in "${sample_path}"/*_filtered.fastq.gz; do
        [ -f "$fastq_file" ] || continue

        # Extraction du nom du génotype (ex: B4, D3, etc.)
        genotype=$(basename "$fastq_file" | sed 's/_filtered.fastq.gz//')

        echo "---------------------------------------------------------------------"
        echo "Traitement : $sample_id | Génotype : $genotype"
        echo "---------------------------------------------------------------------"

        # Dossiers temporaires et définitifs pour cet échantillon
        SAMPLE_MEDAKA_DIR="${MEDAKA_DIR}/${sample_id}/${genotype}"
        SAMPLE_BAM_DIR="${BAM_DIR}/${sample_id}"
        SAMPLE_VCF_DIR="${VCF_DIR}/${sample_id}"
        SAMPLE_CONS_DIR="${CONSENSUS_DIR}/${sample_id}"

        mkdir -p "$SAMPLE_MEDAKA_DIR" "$SAMPLE_BAM_DIR" "$SAMPLE_VCF_DIR" "$SAMPLE_CONS_DIR"

        # Fichier temporaire de référence isolé juste pour cet édit Medaka
        TMP_REF="${SAMPLE_MEDAKA_DIR}/ref_${genotype}.fasta"
        samtools faidx "$REF_GLOBAL" "$genotype" > "$TMP_REF" 2>/dev/null || true

        if [ ! -s "$TMP_REF" ]; then
            echo "❌ Erreur : Génotype $genotype non trouvé dans $REF_GLOBAL. Sauté."
            rm -rf "$SAMPLE_MEDAKA_DIR"
            continue
        fi

        # -----------------------------------------------------------------
        # 1. POLISSAGE ET CONSENSUS AVEC MEDAKA (.fasta + .fai)
        # -----------------------------------------------------------------
        echo "--> [1/4] Génération du consensus Medaka..."

        if ! medaka_consensus -i "$fastq_file" -d "$TMP_REF" -o "$SAMPLE_MEDAKA_DIR" -m "$MEDAKA_MODEL" -t "$THREADS" -f > /dev/null 2>&1; then
            echo "⚠️ Échec avec $MEDAKA_MODEL. Tentative avec le modèle fallback..."
            medaka_consensus -i "$fastq_file" -d "$TMP_REF" -o "$SAMPLE_MEDAKA_DIR" -m "r941_min_fast_g303" -t "$THREADS" -f > /dev/null 2>&1 || true
        fi

        MEDAKA_OUTPUT="${SAMPLE_MEDAKA_DIR}/consensus.fasta"
        CONS_FINAL="${SAMPLE_CONS_DIR}/${genotype}_cons.fasta"

        if [ -s "$MEDAKA_OUTPUT" ]; then
            # Header unique propre
            echo ">${sample_id}_${genotype}_medaka" > "$CONS_FINAL"
            grep -v "^>" "$MEDAKA_OUTPUT" >> "$CONS_FINAL"
            
            # Indexation du FASTA (.fai)
            samtools faidx "$CONS_FINAL"
            echo "    ✅ FASTA & .fai créés : $CONS_FINAL"
        else
            echo "    ❌ Erreur : Medaka n'a pas produit de fichier consensus pour $sample_id ($genotype)."
            continue
        fi

        # Nettoyage du fichier temporaire de référence
        rm -f "$TMP_REF" "${TMP_REF}.fai"

        # -----------------------------------------------------------------
        # 2. ALIGNEMENT ET INDEXATION BAM (.bam + .bai)
        # -----------------------------------------------------------------
        FINAL_BAM="${SAMPLE_BAM_DIR}/${genotype}_sorted.bam"

        echo "--> [2/4] Re-alignement des reads sur le consensus Medaka (minimap2)..."
        minimap2 -ax splice "$CONS_FINAL" "$fastq_file" 2>/dev/null | samtools sort -o "$FINAL_BAM"
        samtools index "$FINAL_BAM"
        echo "    ✅ BAM & .bai créés : $FINAL_BAM"

        # -----------------------------------------------------------------
        # 3. COMPTAGE EXHAUSTIF POSITION PAR POSITION (.vcf)
        # -----------------------------------------------------------------
        FINAL_VCF="${SAMPLE_VCF_DIR}/${genotype}_sammpileup.vcf"

        echo "--> [3/4] Génération du VCF exhaustif position par position (mpileup)..."
        samtools mpileup -d "$MPILEUP_DEPTH" -Q "$MPILEUP_QUAL" -f "$CONS_FINAL" "$FINAL_BAM" > "$FINAL_VCF" 2>/dev/null
        echo "    ✅ VCF exhaustif créé : $FINAL_VCF"

        # -----------------------------------------------------------------
        # 4. EXTRACTION DU TABLEAU DE VARIANTS (_variants.txt)
        # -----------------------------------------------------------------
        VARIANTS_TXT="${SAMPLE_VCF_DIR}/${genotype}_sammpileup.vcf_variants.txt"

        echo "--> [4/4] Extraction du tableau des variants (_variants.txt)..."
        awk 'BEGIN {FS="\t"; OFS="\t"} 
             !/^#/ && $5 != "." && $5 != "" {
                 print $1, $2, $4, $5, $6, $7
             }' "$FINAL_VCF" > "$VARIANTS_TXT"
        echo "    ✅ Fichier variants créé : $VARIANTS_TXT"

    done
done

# Concaténation globale de tous les consensus du run
cat "${CONSENSUS_DIR}"/*/*_cons.fasta > "${CONSENSUS_DIR}/all_cons.fasta" 2>/dev/null || true

echo "====================================================================="
echo "  Étape 03 terminée avec succès !"
echo "====================================================================="