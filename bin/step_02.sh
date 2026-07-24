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
BLAST_DIR="${RESULTS_DIR}/04_BLASTN_ANALYSIS"
REFILTERED_DIR="${RESULTS_DIR}/05_REFILTERED_FASTQ"

# Paramètres de VIRiONT
BITSCORE_MIN=200
MI_CUTOFF=10 # Seuil Multi-Infection (%)

# Création des dossiers de sortie
mkdir -p "$SUPDATA_REFSEQ" "$SUPDATA_DB" "$TRIMMED_DIR" "$BLAST_DIR" "$REFILTERED_DIR"

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
# 2. PRÉPARATION GLOBALE DE 00_SUPDATA (REFSEQ ET DB BLAST)
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
    echo "--> [00_SUPDATA/DB] Indexation BLAST (makeblastdb)..."
    makeblastdb -in "$TARGET_REFSEQ" -out "$DB_PREFIX" -input_type fasta -dbtype nucl > /dev/null
fi

echo "====================================================================="
echo "  EXECUTION ÉTAPE 02 (Filtrage & Sélection Génotype + Multi-Infection)"
echo "  Fichier CSV   : $(basename "$CSV_FILE")"
echo "  Virus / Tech  : $virus_name | $tech_name"
echo "  Réf BLAST     : $ref_filename"
echo "  Filtre Taille : -l $min_length --maxlength $max_length"
echo "  Seuil MI      : >= ${MI_CUTOFF}%"
echo "====================================================================="

# =====================================================================
# 3. BOUCLE DE TRAITEMENT SUR LES SAMPLES
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
    FASTA_FILE="${BLAST_DIR}/${sample_id}_trimmed.fasta"
    BLAST_RAW="${BLAST_DIR}/${sample_id}_fmt.txt"
    BEST_HITS_FILE="${BLAST_DIR}/${sample_id}_best_hits.tsv"

    SAMPLE_REFILTER_DIR="${REFILTERED_DIR}/${sample_id}"
    mkdir -p "$SAMPLE_REFILTER_DIR"

    echo "---------------------------------------------------------------------"
    echo "Traitement $sample_id"
    echo "---------------------------------------------------------------------"

    # --- ÉTAPE 03 : Chopper ---
    echo "--> [03_FILTERED_TRIMMED] Filtration par taille (chopper)..."
    gunzip -c "$DEHOSTED_FASTQ" | chopper -l "$min_length" --maxlength "$max_length" | gzip -c > "$TRIMMED_FASTQ"

    if [ ! -s "$TRIMMED_FASTQ" ]; then
        echo "⚠️ Aucun read conservé par chopper pour $sample_id."
        continue
    fi

    # --- ÉTAPE 04 : BLASTn ---
    echo "--> [04_BLASTN_ANALYSIS] Conversion FASTQ -> FASTA (seqkit)..."
    seqkit fq2fa "$TRIMMED_FASTQ" -o "$FASTA_FILE"

    echo "--> [04_BLASTN_ANALYSIS] Alignement BLASTn..."
    blastn -query "$FASTA_FILE" \
           -db "$DB_PREFIX" \
           -outfmt "6 qseqid sseqid bitscore slen qlen length pident" \
           -num_threads 4 \
           -out "$BLAST_RAW"

    if [ ! -s "$BLAST_RAW" ]; then
        echo "⚠️ Aucun match BLAST pour $sample_id."
        rm -f "$FASTA_FILE"
        continue
    fi

    # --- ÉTAPE 05 : Filtrage Bitscore MAX & Analyse Multi-Infection ---
    echo "--> [05_REFILTERED_FASTQ] Analyse des génotypes (Règle Multi-Infection)..."
    
    # 1. Filtre Bitscore >= 200 + Meilleur Hit Unique par Read
    awk -v min_score="$BITSCORE_MIN" '$3 >= min_score' "$BLAST_RAW" | \
    sort -k1,1 -k3,3nr | \
    sort -u -k1,1 > "$BEST_HITS_FILE"

    if [ ! -s "$BEST_HITS_FILE" ]; then
        echo "⚠️ Aucun match au-dessus du bitscore min ($BITSCORE_MIN) pour $sample_id."
        rm -f "$FASTA_FILE" "$BEST_HITS_FILE"
        continue
    fi

    # 2. Total des reads valides
    TOTAL_READS=$(wc -l < "$BEST_HITS_FILE")
    echo "    📊 Total reads alignés : $TOTAL_READS"

    # 3. Traitement de TOUS les génotypes >= MI_CUTOFF (%)
    cut -f2 "$BEST_HITS_FILE" | sort | uniq -c | sort -nr | while read -r count gen; do
        
        # Calcul du pourcentage
        pct=$(awk -v c="$count" -v t="$TOTAL_READS" 'BEGIN { printf "%.2f", (c/t)*100 }')
        is_above=$(awk -v p="$pct" -v cutoff="$MI_CUTOFF" 'BEGIN { print (p >= cutoff) ? "YES" : "NO" }')

        if [ "$is_above" == "YES" ]; then
            echo "    🏆 Génotype sélectionné : $gen ($count reads, $pct%)"

            READ_IDS_FILE="${BLAST_DIR}/${sample_id}_${gen}_readlist.txt"
            awk -v g="$gen" '$2 == g {print $1}' "$BEST_HITS_FILE" > "$READ_IDS_FILE"

            FINAL_GENOTYPE_FASTQ="${SAMPLE_REFILTER_DIR}/${gen}_filtered.fastq.gz"
            seqkit grep -f "$READ_IDS_FILE" "$TRIMMED_FASTQ" -o "$FINAL_GENOTYPE_FASTQ"
            
            rm -f "$READ_IDS_FILE"
            echo "       ✅ FASTQ généré : $(basename "$FINAL_GENOTYPE_FASTQ")"
        else
            echo "    ⏩ Génotype ignoré (sous le seuil de ${MI_CUTOFF}%) : $gen ($count reads, $pct%)"
        fi
    done

    # Nettoyage des fichiers temporaires
    rm -f "$FASTA_FILE" "$BEST_HITS_FILE"

done

echo "====================================================================="
echo "  Étape 02 terminée avec succès !"
echo "====================================================================="