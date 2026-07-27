#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS
# =====================================================================
RESULTS_DIR="results_test"
REF_DIR="${RESULTS_DIR}/00_SUPDATA/REFSEQ"
AMPLICON_DIR="${RESULTS_DIR}/04_AMPLICON_SORTER"
GENOTYPING_DIR="${RESULTS_DIR}/05_FINAL_GENOTYPING"

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
echo "  EXECUTION ÉTAPE 03 (Polissage Medaka & Calling Variants)"
echo "  Modèle Medaka : $MEDAKA_MODEL"
echo "====================================================================="

# Localisation du fichier FASTA de référence dans 00_SUPDATA
REF_GLOBAL=$(find "$REF_DIR" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) | head -n 1)

if [ -z "$REF_GLOBAL" ] || [ ! -f "$REF_GLOBAL" ]; then
    echo "❌ Erreur : Fichier de référence introuvable dans $REF_DIR"
    exit 1
fi

# Indexation de la référence globale si nécessaire
if [ ! -f "${REF_GLOBAL}.fai" ]; then
    samtools faidx "$REF_GLOBAL"
fi

# =====================================================================
# BOUCLE SUR TOUS LES FIFIERS DE BLAST (ISSUS DE L'ÉTAPE 02)
# =====================================================================
for blast_file in "${GENOTYPING_DIR}"/*_blast.tsv; do
    [ -f "$blast_file" ] || continue

    # Extraction du nom du sample et du cluster depuis le nom du fichier BLAST
    # Format attendu: barcode_SAMPLE-X_cluster_N_blast.tsv
    filename=$(basename "$blast_file")
    sample_id=$(echo "$filename" | sed -E 's/_(gene_|species_|cluster_)[0-9]+_.*_blast\.tsv//')
    cluster_id=$(echo "$filename" | sed 's/_blast.tsv//' | sed "s/${sample_id}_//")

    # Génotype identifié par BLAST à l'étape 02
    genotype=$(awk '{print $2}' "$blast_file" | head -n 1)

    if [ -z "$genotype" ]; then
        echo "⚠️ Aucun génotype trouvé dans $blast_file. Sauté."
        continue
    fi

    echo "---------------------------------------------------------------------"
    echo "Traitement : $sample_id | Cluster : $cluster_id | Génotype : $genotype"
    echo "---------------------------------------------------------------------"

    # Recherche du FASTQ/FASTA des reads de ce cluster produit par Amplicon_sorter
    READS_FILE=$(find "${AMPLICON_DIR}/${sample_id}" -name "*${cluster_id}*" \( -name "*.fastq*" -o -name "*.fasta*" \) ! -name "*consensus*" ! -name "*unique*" | head -n 1 || true)

    if [ -z "$READS_FILE" ] || [ ! -f "$READS_FILE" ]; then
        echo "⚠️ Fichier de reads introuvable pour $sample_id ($cluster_id) dans $AMPLICON_DIR. Sauté."
        continue
    fi

    # Dossiers temporaires et définitifs pour cet échantillon
    SAMPLE_MEDAKA_DIR="${MEDAKA_DIR}/${sample_id}/${cluster_id}_${genotype}"
    SAMPLE_BAM_DIR="${BAM_DIR}/${sample_id}"
    SAMPLE_VCF_DIR="${VCF_DIR}/${sample_id}"
    SAMPLE_CONS_DIR="${CONSENSUS_DIR}/${sample_id}"

    mkdir -p "$SAMPLE_MEDAKA_DIR" "$SAMPLE_BAM_DIR" "$SAMPLE_VCF_DIR" "$SAMPLE_CONS_DIR"

    # Isolation de la séquence de référence du génotype sélectionné
    TMP_REF="${SAMPLE_MEDAKA_DIR}/ref_${genotype}.fasta"
    samtools faidx "$REF_GLOBAL" "$genotype" > "$TMP_REF" 2>/dev/null || true

    if [ ! -s "$TMP_REF" ]; then
        echo "⚠️ Génotype [$genotype] non trouvé dans $REF_GLOBAL. Utilisation du consensus de-novo comme référence."
        # Fallback: Utiliser le consensus de-novo d'Amplicon_sorter comme référence
        DE_NOVO_FASTA=$(find "${AMPLICON_DIR}/${sample_id}" -name "*${cluster_id}*.fasta" | head -n 1)
        cp "$DE_NOVO_FASTA" "$TMP_REF"
    fi

    # -----------------------------------------------------------------
    # 1. POLISSAGE ET CONSENSUS AVEC MEDAKA
    # -----------------------------------------------------------------
    echo "--> [1/4] Polissage Medaka..."

    if ! medaka_consensus -i "$READS_FILE" -d "$TMP_REF" -o "$SAMPLE_MEDAKA_DIR" -m "$MEDAKA_MODEL" -t "$THREADS" -f > /dev/null 2>&1; then
        echo "⚠️ Échec avec $MEDAKA_MODEL. Tentative avec le modèle fallback (r941_min_fast_g303)..."
        medaka_consensus -i "$READS_FILE" -d "$TMP_REF" -o "$SAMPLE_MEDAKA_DIR" -m "r941_min_fast_g303" -t "$THREADS" -f > /dev/null 2>&1 || true
    fi

    MEDAKA_OUTPUT="${SAMPLE_MEDAKA_DIR}/consensus.fasta"
    CONS_FINAL="${SAMPLE_CONS_DIR}/${cluster_id}_${genotype}_cons.fasta"

    if [ -s "$MEDAKA_OUTPUT" ]; then
        # Header unique propre
        echo ">${sample_id}_${genotype}_medaka" > "$CONS_FINAL"
        grep -v "^>" "$MEDAKA_OUTPUT" >> "$CONS_FINAL"
        
        # Indexation du FASTA (.fai)
        samtools faidx "$CONS_FINAL"
        echo "    ✅ FASTA & .fai créés : $(basename "$CONS_FINAL")"
    else
        echo "    ❌ Erreur : Medaka n'a pas produit de consensus pour $sample_id ($cluster_id)."
        rm -rf "$SAMPLE_MEDAKA_DIR"
        continue
    fi

    # Nettoyage de la référence temporaire
    rm -f "$TMP_REF" "${TMP_REF}.fai"

    # -----------------------------------------------------------------
    # 2. ALIGNEMENT ET INDEXATION BAM (.bam + .bai)
    # -----------------------------------------------------------------
    FINAL_BAM="${SAMPLE_BAM_DIR}/${cluster_id}_${genotype}_sorted.bam"

    echo "--> [2/4] Re-alignement des reads sur le consensus Medaka (minimap2)..."
    minimap2 -ax splice "$CONS_FINAL" "$READS_FILE" 2>/dev/null | samtools sort -o "$FINAL_BAM"
    samtools index "$FINAL_BAM"
    echo "    ✅ BAM & .bai créés : $(basename "$FINAL_BAM")"

    # -----------------------------------------------------------------
    # 3. COMPTAGE EXHAUSTIF POSITION PAR POSITION (.vcf)
    # -----------------------------------------------------------------
    FINAL_VCF="${SAMPLE_VCF_DIR}/${cluster_id}_${genotype}_sammpileup.vcf"

    echo "--> [3/4] Génération du VCF exhaustif (samtools mpileup)..."
    samtools mpileup -d "$MPILEUP_DEPTH" -Q "$MPILEUP_QUAL" -f "$CONS_FINAL" "$FINAL_BAM" > "$FINAL_VCF" 2>/dev/null
    echo "    ✅ VCF exhaustif créé : $(basename "$FINAL_VCF")"

    # -----------------------------------------------------------------
    # 4. EXTRACTION DU TABLEAU DE VARIANTS (_variants.txt)
    # -----------------------------------------------------------------
    VARIANTS_TXT="${SAMPLE_VCF_DIR}/${cluster_id}_${genotype}_variants.txt"

    echo "--> [4/4] Extraction du tableau de variants (_variants.txt)..."
    awk 'BEGIN {FS="\t"; OFS="\t"} 
         !/^#/ && $5 != "." && $5 != "" {
             print $1, $2, $4, $5, $6, $7
         }' "$FINAL_VCF" > "$VARIANTS_TXT"
    echo "    ✅ Fichier variants créé : $(basename "$VARIANTS_TXT")"

done

# Concaténation globale de tous les consensus du run
cat "${CONSENSUS_DIR}"/*/*_cons.fasta > "${CONSENSUS_DIR}/all_cons.fasta" 2>/dev/null || true

echo "====================================================================="
echo " 🎉 Étape 03 terminée avec succès !"
echo "====================================================================="