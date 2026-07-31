#!/bin/bash

# Interruption si une commande échoue
set -eo pipefail

# =====================================================================
# 1. DÉFINITION DES CHEMINS, DOSSIERS ET CACHE
# =====================================================================
CSV_FILE="fastq_pass/260630_VHB_WG_ABS.csv"
DATA_DIR="fastq_pass"
MERGED_DIR="results_test/01_MERGED"
DEHOST_DIR="results_test/02_DEHOSTING"

# Dossiers de cache locaux pour éviter tout conflit /tmp ou permissions
export HOSTILE_CACHE_DIR="/srv/scratch/chu-lyon.fr/alamiso/VIRiONT_V2/ref_hostile"
export TMPDIR="/srv/scratch/chu-lyon.fr/alamiso/VIRiONT_V2/tmp"

# Création des dossiers requis
mkdir -p "$MERGED_DIR" "$DEHOST_DIR" "$HOSTILE_CACHE_DIR" "$TMPDIR"

if [ ! -f "$CSV_FILE" ]; then
    echo "❌ Erreur : Fichier CSV introuvable ($CSV_FILE)."
    exit 1
fi

if ! command -v hostile &> /dev/null; then
    echo "❌ Erreur : 'hostile' n'est pas disponible dans cet environnement."
    exit 1
fi

# =====================================================================
# 2. BOUCLE DE TRAITEMENT SUR LES ÉCHANTILLONS (DEHOSTING HOSTILE)
# =====================================================================
tail -n +2 "$CSV_FILE" | while IFS=';' read -r plate_pos sample component forward reverse flowcell kit_id primers || [ -n "$sample" ]; do

    # Nettoyage des espaces et retours chariot
    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')

    [ -z "$sample" ] && continue

    # Extraction du numéro de barcode
    num_barcode=$(echo "$component" | grep -o '[0-9]\+' || true)

    if [ -z "$num_barcode" ]; then
        echo "⚠️ Barcode non valide pour '$component'. Ligne sautée."
        continue
    fi

    formatted_barcode=$(printf "%02d" "$num_barcode")
    folder_name="barcode${formatted_barcode}"
    sample_id="barcode_${formatted_barcode}_${sample}"

    if [ ! -d "$DATA_DIR/$folder_name" ]; then
        folder_name="barcode${num_barcode}"
        if [ ! -d "$DATA_DIR/$folder_name" ]; then
            echo "⚠️ Dossier $DATA_DIR/$folder_name introuvable. Sauté."
            continue
        fi
    fi

    echo "--------------------------------------------------"
    echo "Traitement $sample_id ($folder_name)"
    echo "--------------------------------------------------"

    # --- ÉTAPE 1 : Fusion FASTQ ---
    merged_output="${MERGED_DIR}/${sample_id}_merged.fastq.gz"
    echo "--> [01_MERGED] Fusion des fichiers FASTQ..."
    
    gz_count=$(find "$DATA_DIR/$folder_name" -type f -name "*.gz" | wc -l)
    if [ "$gz_count" -gt 0 ]; then
        zcat "$DATA_DIR/$folder_name"/* | gzip -c > "$merged_output"
    else
        cat "$DATA_DIR/$folder_name"/* | gzip -c > "$merged_output"
    fi

    # --- ÉTAPE 2 : Dehosting Hostile ---
    FINAL_OUTPUT="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"

    echo "--> [02_DEHOSTING] Décontamination avec Hostile..."
    
    hostile clean \
        --fastq1 "$merged_output" \
        --index human-t2t-hla.rs-viral-202401_ml-phage-202401 \
        --aligner minimap2 \
        --threads 32 \
        -o "$DEHOST_DIR"

    # Nom exact généré par Hostile (ex: barcode_01_26104456601_merged.clean.fastq.gz)
    hostile_default_file="${DEHOST_DIR}/${sample_id}_merged.clean.fastq.gz"
    
    if [ -f "$hostile_default_file" ]; then
        mv "$hostile_default_file" "$FINAL_OUTPUT"
        echo "✅ Créé avec succès : $(basename "$FINAL_OUTPUT")"
    else
        # Sécurité supplémentaire en cas de variante de nom
        alt_file=$(find "$DEHOST_DIR" -maxdepth 1 -name "*${sample_id}*.clean.fastq.gz" | head -n 1)
        if [ -n "$alt_file" ] && [ -f "$alt_file" ]; then
            mv "$alt_file" "$FINAL_OUTPUT"
            echo "✅ Créé avec succès : $(basename "$FINAL_OUTPUT")"
        else
            echo "❌ Erreur : Le fichier Hostile n'a pas été trouvé."
            exit 1
        fi
    fi

done

echo "====================================================================="
echo " 🎉 Déhosting avec Hostile terminé pour tous les échantillons !"
echo "====================================================================="