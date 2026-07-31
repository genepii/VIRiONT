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
tail -n +2 "$CSV_FILE" | while IFS=';' read -r plate_pos sample component forward reverse flowcell kit_id primers || [ -n "$sample" ]; do

    # Nettoyage direct des caractères invisibles (\r, \n) et des espaces
    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')

    # Passer les lignes vides
    [ -z "$sample" ] && continue

    # Extraction sécurisée du numéro de barcode
    num_barcode=$(echo "$component" | grep -o '[0-9]\+' || true)

    if [ -z "$num_barcode" ]; then
        echo "⚠️ Barcode non valide trouvé dans le composant '$component'. Ligne sautée."
        continue
    fi

    # ---------------------------------------------------------------------
    # NOMENCLATURE : barcode_<NUMERO_BARCODE_2DIGITS>_<NUMERO_GLIMS>
    # Exemple : barcode_01_26104456601
    # ---------------------------------------------------------------------
    formatted_barcode=$(printf "%02d" "$num_barcode")
    folder_name="barcode${formatted_barcode}"
    sample_id="barcode_${formatted_barcode}_${sample}"

    if [ ! -d "$DATA_DIR/$folder_name" ]; then
        # Essai alternatif si le dossier est nommé sans zéro initial (ex: barcode1 au lieu de barcode01)
        folder_name="barcode${num_barcode}"
        if [ ! -d "$DATA_DIR/$folder_name" ]; then
            echo "⚠️ Dossier $DATA_DIR/$folder_name introuvable. Sauté."
            continue
        fi
    fi

    echo "--------------------------------------------------"
    echo "Traitement $sample_id ($folder_name)"
    echo "--------------------------------------------------"

    # --- ÉTAPE 1 : Fusion des FASTQ par Barcode ---
    merged_output="${MERGED_DIR}/${sample_id}_merged.fastq.gz"
    echo "--> [01_MERGED] Fusion des fichiers FASTQ..."
    
    # Gestion automatique selon si les fichiers sources sont compressés (.gz) ou non
    gz_count=$(find "$DATA_DIR/$folder_name" -type f -name "*.gz" | wc -l)
    if [ "$gz_count" -gt 0 ]; then
        zcat "$DATA_DIR/$folder_name"/* | gzip -c > "$merged_output"
    else
        cat "$DATA_DIR/$folder_name"/* | gzip -c > "$merged_output"
    fi

    # --- ÉTAPE 2 : Dehosting tolérant (Alignement -ax splice & Filtre MAPQ) ---
    FINAL_OUTPUT="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"

    echo "--> [02_DEHOSTING] Alignement GRCh38 (-ax splice) & extraction tolérante des reads (non-alignés + MAPQ < 10)..."
    
    # Explication du pipeline ci-dessous :
    # 1. minimap2 aligne avec la sensibilité 'splice'
    # 2. samtools view ne garde QUE les reads totalement non-alignés OR avec un score de qualité d'alignement humain très faible (MAPQ < 10)
    # 3. samtools fastq exclut les alignements secondaires/supplémentaires (-F 0x900) pour éviter les doublons tout en sauvant le read d'origine
    minimap2 -t 4 -ax splice "$REF_INDEX" "$merged_output" | \
    samtools view -e 'flag.unmap || mapq < 10' -u - | \
    samtools fastq -@ 4 -F 0x900 - | gzip -c > "$FINAL_OUTPUT"

    echo "✅ Créé avec succès : $(basename "$FINAL_OUTPUT")"

done

echo "====================================================================="
echo " 🎉 Déhosting terminé avec succès !"
echo "====================================================================="