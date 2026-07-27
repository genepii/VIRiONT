#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS
# =====================================================================
CSV_FILE="fastq_pass/260630_VHB_WG_ABS.csv"
REF_SOURCE_DIR="/srv/scratch/chu-lyon.fr/alamiso/VIRiONT_V2/ref"

RESULTS_DIR="results_test"
SUPDATA_REFSEQ="${RESULTS_DIR}/00_SUPDATA/REFSEQ"
SUPDATA_DB="${RESULTS_DIR}/00_SUPDATA/DB"
DEHOST_DIR="${RESULTS_DIR}/02_DEHOSTING"
TRIMMED_DIR="${RESULTS_DIR}/03_FILTERED_TRIMMED"
AMPLICON_DIR="${RESULTS_DIR}/04_AMPLICON_SORTER"
GENOTYPING_DIR="${RESULTS_DIR}/05_FINAL_GENOTYPING"

# Création des dossiers de sortie
mkdir -p "$SUPDATA_REFSEQ" "$SUPDATA_DB" "$TRIMMED_DIR" "$AMPLICON_DIR" "$GENOTYPING_DIR"

if [ ! -f "$CSV_FILE" ]; then
    echo "❌ Erreur : Fichier CSV introuvable ($CSV_FILE)."
    exit 1
fi

# -----------------------------------------------------------------
# 1. DÉTECTION GLOBALE DU VIRUS ET DE LA TECHNO DEPUIS LE CSV
# -----------------------------------------------------------------
virus_name=""
tech_name=""

csv_upper=$(echo "$CSV_FILE" | tr '[:lower:]' '[:upper:]')

if [[ "$csv_upper" == *"VHD"* ]]; then virus_name="VHD"; fi
if [[ "$csv_upper" == *"VHB"* ]]; then virus_name="VHB"; fi

if [[ "$csv_upper" == *"R0"* ]]; then tech_name="R0"; fi
if [[ "$csv_upper" == *"WG"* ]]; then tech_name="WG"; fi
if [[ "$csv_upper" == *"POL"* ]]; then tech_name="POL"; fi

if [ "$virus_name" == "VHD" ] && [ "$tech_name" == "R0" ]; then
    ref_filename="Ref_HDV_subtype_R0.fasta"
    min_length=300
    max_length=800
elif [ "$virus_name" == "VHD" ] && [ "$tech_name" == "WG" ]; then
    ref_filename="HDV_subtype_V2.fasta"
    min_length=1000
    max_length=2000
elif [ "$virus_name" == "VHB" ] && [ "$tech_name" == "POL" ]; then
    ref_filename="HBV_subtype_R0.fasta"
    min_length=800
    max_length=5000
elif [ "$virus_name" == "VHB" ] && [ "$tech_name" == "WG" ]; then
    ref_filename="HBV_subtype.fasta"
    min_length=1000
    max_length=5000
else
    echo "❌ Erreur : Impossible de déterminer le Virus/Techno depuis le nom du CSV ($CSV_FILE)."
    exit 1
fi

ORIGINAL_REF="${REF_SOURCE_DIR}/${ref_filename}"
TARGET_REFSEQ="${SUPDATA_REFSEQ}/${ref_filename}"
DB_PREFIX="${SUPDATA_DB}/${ref_filename%.fasta}"

# -----------------------------------------------------------------
# 2. PRÉPARATION GLOBALE DE 00_SUPDATA (INDEXATION BASE GÉNOTYPES)
# -----------------------------------------------------------------
if [ ! -f "$ORIGINAL_REF" ]; then
    echo "❌ Fichier source $ORIGINAL_REF introuvable dans $REF_SOURCE_DIR !"
    exit 1
fi

if [ ! -f "$TARGET_REFSEQ" ]; then
    echo "--> [00_SUPDATA/REFSEQ] Copie de la référence : $ref_filename"
    cp "$ORIGINAL_REF" "$TARGET_REFSEQ"
fi

if [ ! -f "${DB_PREFIX}.nhr" ]; then
    echo "--> [00_SUPDATA/DB] Indexation BLAST de la base interne des génotypes..."
    makeblastdb -in "$TARGET_REFSEQ" -out "$DB_PREFIX" -input_type fasta -dbtype nucl > /dev/null
fi

echo "====================================================================="
echo "  EXECUTION ÉTAPE 02 (Chopper -> Amplicon_sorter -> BLAST Genotyping)"
echo "  Fichier CSV    : $(basename "$CSV_FILE")"
echo "  Virus / Tech   : $virus_name | $tech_name"
echo "  Base Génotypes : $ref_filename"
echo "  Filtre Taille  : Min = $min_length pb | Max = $max_length pb"
echo "====================================================================="

# =====================================================================
# 3. BOUCLE DE TRAITEMENT SUR LES ÉCHANTILLONS
# =====================================================================
tail -n +2 "$CSV_FILE" | while IFS=';' read -r plate_pos sample component forward reverse flowcell kit_id primers; do
    
    [ -z "$sample" ] && continue

    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')

    num_barcode=$(echo "$component" | sed 's/[^0-9]//g')
    sample_id="barcode_${sample}-${num_barcode}"

    DEHOSTED_FASTQ="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"
    if [ ! -f "$DEHOSTED_FASTQ" ]; then
        DEHOSTED_FASTQ="${DEHOST_DIR}/${sample_id}_meta.fastq.gz"
    fi

    if [ ! -f "$DEHOSTED_FASTQ" ]; then
        echo "⚠️ Output de step_01 introuvable pour $sample_id. Sauté."
        continue
    fi

    TRIMMED_FASTQ="${TRIMMED_DIR}/${sample_id}_trimmed.fastq.gz"
    SAMPLE_AMPLICON_DIR="${AMPLICON_DIR}/${sample_id}"

    echo "---------------------------------------------------------------------"
    echo "Traitement : $sample_id"
    echo "---------------------------------------------------------------------"

    # --- ÉTAPE 03 : Chopper (Filtrage par longueur) ---
    echo "--> [03_FILTERED_TRIMMED] Filtration par taille avec chopper..."
    gunzip -c "$DEHOSTED_FASTQ" | chopper -l "$min_length" --maxlength "$max_length" | gzip -c > "$TRIMMED_FASTQ"

    if [ ! -s "$TRIMMED_FASTQ" ]; then
        echo "⚠️ Aucun read conservé par chopper pour $sample_id."
        rm -f "$TRIMMED_FASTQ"
        continue
    fi

    # --- ÉTAPE 04 : Amplicon_sorter (Clustering De Novo WG sur TOUS les reads) ---
    echo "--> [04_AMPLICON_SORTER] Clustering & Génération du consensus De Novo..."
    mkdir -p "$SAMPLE_AMPLICON_DIR"

    amplicon_sorter.py \
        -i "$TRIMMED_FASTQ" \
        -o "$SAMPLE_AMPLICON_DIR" \
        -min "$min_length" \
        -max "$max_length" \
        -maxr 1000000 \
        -ar \
        -ssg 85 \
        -ss 85 \
        -sc 92 \
        -np 4

    CONSENSUS_FILES=$(find "$SAMPLE_AMPLICON_DIR" -maxdepth 2 \( -name "*.fasta" -o -name "*.fa" \) ! -name "*unique*" || true)

    if [ -z "$CONSENSUS_FILES" ]; then
        echo "⚠️ Aucun consensus produit par amplicon_sorter pour $sample_id."
        rm -f "$TRIMMED_FASTQ"
        continue
    fi

    # --- ÉTAPE 05 : BLASTn & Sélection du Cluster Maître ---
    echo "--> [05_FINAL_GENOTYPING] Identification des génotypes via BLASTn..."
    
    COMBINED_BLAST="${GENOTYPING_DIR}/${sample_id}_all_consensus_blast.tsv"
    rm -f "$COMBINED_BLAST"

    # Exécution de BLASTn sur chaque consensus et capture sécurisée du 1er hit
    for cons_file in $CONSENSUS_FILES; do
        blastn -query "$cons_file" \
               -db "$DB_PREFIX" \
               -outfmt "6 qseqid sseqid pident length evalue bitscore" \
               -max_hsps 1 2>/dev/null | head -n 1 >> "$COMBINED_BLAST" || true
    done

    BEST_HIT_FILE="${GENOTYPING_DIR}/${sample_id}_best_genotype.tsv"
    
    if [ -s "$COMBINED_BLAST" ]; then
        sort -k6,6nr -k4,4nr "$COMBINED_BLAST" | head -n 1 > "$BEST_HIT_FILE"
    fi

    if [ -s "$BEST_HIT_FILE" ]; then
        best_cluster=$(awk '{print $1}' "$BEST_HIT_FILE")
        genotype=$(awk '{print $2}' "$BEST_HIT_FILE")
        pident=$(awk '{print $3}' "$BEST_HIT_FILE")
        length=$(awk '{print $4}' "$BEST_HIT_FILE")
        
        echo "   🎯 Cluster Maître Retenu : $best_cluster"
        echo "   🏆 Génotype Identifié   : $genotype"
        echo "   📊 Qualité Alignement  : Identité = ${pident}% | Longueur = ${length} pb"
        
        MASTER_FASTA="${GENOTYPING_DIR}/${sample_id}_${genotype}_master_consensus.fasta"
        
        # Extraction du FASTA consensus maître depuis le fichier global
        awk -v target=">$best_cluster" '
            $0 ~ target {flag=1; print; next}
            /^>/ {flag=0}
            flag {print}
        ' "$SAMPLE_AMPLICON_DIR"/*consensusfile.fasta > "$MASTER_FASTA" 2>/dev/null || true

        # Fallback de sécurité si l'extraction par AWK échoue
        if [ ! -s "$MASTER_FASTA" ]; then
            cp $(echo "$CONSENSUS_FILES" | head -n 1) "$MASTER_FASTA"
        fi
        
        echo "   ✅ Fichier Maître créé  : $(basename "$MASTER_FASTA")"
    else
        echo "⚠️ Aucun résultat BLAST valide obtenu pour $sample_id."
    fi

    # Nettoyage du fichier FASTQ temporaire
    rm -f "$TRIMMED_FASTQ"

done

echo "====================================================================="
echo " 🎉 Étape 02 terminée avec succès !"
echo "====================================================================="