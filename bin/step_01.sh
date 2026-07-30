#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# 1. DÉFINITION DES CHEMINS ET DOSSIERS
# =====================================================================
CSV_FILE="fastq_pass/260630_VHB_WG_ABS.csv"
DATA_DIR="fastq_pass"
REF_DIR="ref_GRCh38"
MERGED_DIR="results_test/01_MERGED"
DEHOST_DIR="results_test/02_DEHOSTING"

REF_HUMAN="${REF_DIR}/GRCh38_no_alt.fa"
REF_INDEX="${REF_DIR}/GRCh38_no_alt.mmi"

# Création des dossiers de sortie
mkdir -p "$MERGED_DIR" "$DEHOST_DIR" "$REF_DIR"

if [ ! -f "$CSV_FILE" ]; then
    echo "❌ Erreur : Fichier CSV introuvable ($CSV_FILE)."
    exit 1
fi

# =====================================================================
# 2. TÉLÉCHARGEMENT AUTOMATIQUE ET INDEXATION DU GÉNOME HUMAIN
# =====================================================================
if [ ! -f "$REF_HUMAN" ]; then
    echo "=== [VIRiONT V2] Téléchargement de GRCh38 ==="
    wget -O "$REF_DIR/GRCh38_no_alt.fa.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    echo "=== [VIRiONT V2] Décompression de GRCh38 ==="
    gunzip "$REF_DIR/GRCh38_no_alt.fa.gz"
fi

if [ ! -f "$REF_INDEX" ]; then
    echo "=== [VIRiONT V2] Indexation de GRCh38 avec minimap2 ==="
    minimap2 -d "$REF_INDEX" "$REF_HUMAN"
fi

# =====================================================================
# 3. BOUCLE DE TRAITEMENT SUR LES ÉCHANTILLONS (DEHOSTING)
# =====================================================================
tail -n +2 "$CSV_FILE" | while IFS=';' read -r plate_pos sample component forward reverse flowcell kit_id primers; do
    
    [ -z "$sample" ] && continue

    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')

    num_barcode=$(echo "$component" | sed 's/[^0-9]//g')
    
    # ---------------------------------------------------------------------
    # NOMENCLATURE : barcode_<NUMERO_BARCODE_2DIGITS>_<NUMERO_GLIMS>
    # Exemple : barcode_01_26104456601
    # ---------------------------------------------------------------------
    formatted_barcode=$(printf "%02d" "$num_barcode")
    folder_name="barcode${num_barcode}"
    sample_id="barcode_${formatted_barcode}_${sample}"

    if [ ! -d "$DATA_DIR/$folder_name" ]; then
        echo "⚠️ Dossier $DATA_DIR/$folder_name introuvable. Sauté."
        continue
    fi

    echo "--------------------------------------------------"
    echo "Traitement $sample_id ($folder_name)"
    echo "--------------------------------------------------"
    
    # --- ÉTAPE 1 : Fusion des FASTQ par Barcode ---
    merged_output="${MERGED_DIR}/${sample_id}_merged.fastq.gz"
    echo "--> [01_MERGED] Fusion des fichiers FASTQ..."
    zcat "$DATA_DIR/$folder_name"/* | gzip -c > "$merged_output"

    # --- ÉTAPE 2 : Dehosting (Alignement & Extraction Non-Humain) ---
    HUMAN_BAM="${DEHOST_DIR}/${sample_id}_human.bam"
    FINAL_OUTPUT="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"

    echo "--> [02_DEHOSTING] Alignement GRCh38 (-ax splice)..."
    minimap2 -t 4 -ax splice "$REF_INDEX" "$merged_output" | samtools view -b > "$HUMAN_BAM"

    echo "--> [02_DEHOSTING] Extraction directe des reads non-humains (-f 4 -> FASTQ.GZ)..."
    samtools fastq -f 4 "$HUMAN_BAM" | gzip -c > "$FINAL_OUTPUT"

    # Nettoyage du fichier BAM intermédiaire lourd
    rm -f "$HUMAN_BAM"

    echo "✅ Créé avec succès : $(basename "$FINAL_OUTPUT")"

done

echo "====================================================================="
echo " 🎉 Déhosting terminé avec succès !"
echo "====================================================================="