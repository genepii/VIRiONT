#!/bin/bash

# 1. Définition des chemins
CSV_FILE="fastq_pass/260630_VHB_WG_ABS.csv"
DATA_DIR="fastq_pass"
REF_DIR="ref_GRCh38"
MERGED_DIR="results_test/01_MERGED"
DEHOST_DIR="results_test/02_DEHOSTING"

REF_HUMAN="${REF_DIR}/GRCh38_no_alt.fa"
REF_INDEX="${REF_DIR}/GRCh38_no_alt.mmi"

# Création des dossiers de sortie
mkdir -p "$MERGED_DIR" "$DEHOST_DIR" "$REF_DIR"

# 2. TÉLÉCHARGEMENT AUTOMATIQUE DU GÉROME HUMAIN
if [ ! -f "$REF_HUMAN" ]; then
    echo "=== [VIRiONT V2] Téléchargement de GRCh38 ==="
    wget -O "$REF_DIR/GRCh38_no_alt.fa.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    echo "=== [VIRiONT V2] Décompression de GRCh38 ==="
    gunzip "$REF_DIR/GRCh38_no_alt.fa.gz"
fi

# 3. CRÉATION DE L'INDEX
if [ ! -f "$REF_INDEX" ]; then
    echo "=== [VIRiONT V2] Indexation de GRCh38 avec minimap2 ==="
    minimap2 -d "$REF_INDEX" "$REF_HUMAN"
fi

# 4. Lecture du fichier CSV ligne par ligne
tail -n +2 "$CSV_FILE" | while IFS=';' read -r plate_pos sample component forward reverse flowcell kit_id primers; do
    
    [ -z "$sample" ] && continue

    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')

    num_barcode=$(echo "$component" | sed 's/[^0-9]//g')
    folder_name="barcode${num_barcode}"
    sample_id="barcode_${sample}-${num_barcode}"

    if [ ! -d "$DATA_DIR/$folder_name" ]; then
        echo "⚠️ Dossier $DATA_DIR/$folder_name introuvable. Sauté."
        continue
    fi

    echo "--------------------------------------------------"
    echo "Étape 1 : Fusion pour $folder_name"
    echo "--------------------------------------------------"
    merged_output="${MERGED_DIR}/${folder_name}_merged.fastq.gz"
    
    zcat "$DATA_DIR/$folder_name"/* | gzip -c > "$merged_output"
    
    echo "--> Créé : $merged_output"
    echo "--------------------------------------------------"
    echo "Étape 2 : Dehosting, Filtrage Chimères & Extraction pour $sample_id"
    echo "--------------------------------------------------"

    # Définition des fichiers intermédiaires
    HUMAN_BAM="${DEHOST_DIR}/${sample_id}_human.bam"
    NONHUMAN_BAM="${DEHOST_DIR}/${sample_id}_meta.bam"
    NONHUMAN_FASTQ="${DEHOST_DIR}/${sample_id}_meta.fastq"
    CHIMERA_TRIMMED_FASTQ="${DEHOST_DIR}/${sample_id}_meta_trimmed.fastq"
    FINAL_OUTPUT="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"

    # Alignement avec option -ax splice
    echo "--> Alignement GRCh38 (-ax splice)..."
    minimap2 -t 4 -ax splice "$REF_INDEX" "$merged_output" | samtools view -b > "$HUMAN_BAM"

    # Extraction des reads non-mappés (-f 4)
    echo "--> Extraction des reads non-humains (-f 4)..."
    samtools view -b -f 4 "$HUMAN_BAM" > "$NONHUMAN_BAM"

    # Conversion BAM vers FASTQ via bedtools
    echo "--> Conversion BAM vers FASTQ via bedtools..."
    bedtools bamtofastq -i "$NONHUMAN_BAM" -fq "$NONHUMAN_FASTQ"
    
    # Nettoyage des chimères au milieu du read (Fast-Chimera mode)
    echo "--> Nettoyage des chimères via Porechop..."
    porechop -i "$NONHUMAN_FASTQ" -o "$CHIMERA_TRIMMED_FASTQ" --end_size 0 --middle_threshold 85
    
    # Compression finale du fichier nettoyé
    echo "--> Compression finale..."
    gzip -c "$CHIMERA_TRIMMED_FASTQ" > "$FINAL_OUTPUT"

    # Nettoyage de l'ensemble des fichiers intermédiaires
    rm -f "$HUMAN_BAM" "$NONHUMAN_BAM" "$NONHUMAN_FASTQ" "$CHIMERA_TRIMMED_FASTQ"

    echo "--> Créé avec succès : $FINAL_OUTPUT"

done

echo "=== Déhosting et Filtrage Chimères V2 terminés avec succès ==="