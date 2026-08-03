#!/bin/bash

# Interruption immédiate si une commande du pipeline échoue
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS ET DOSSIERS (VIRiONT V2)
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="fastq_pass"
REF_SOURCE_DIR="/srv/scratch/chu-lyon.fr/alamiso/VIRiONT_V2/ref"

RESULTS_DIR="results_test"
SUPDATA_REFSEQ="${RESULTS_DIR}/00_SUPDATA/REFSEQ"
SUPDATA_DB="${RESULTS_DIR}/00_SUPDATA/DB"
DEHOST_DIR="${RESULTS_DIR}/02_DEHOSTING"
TRIMMED_DIR="${RESULTS_DIR}/03_FILTERED_TRIMMED"
AMPLICON_DIR="${RESULTS_DIR}/04_AMPLICON_SORTER"
GENOTYPING_DIR="${RESULTS_DIR}/05_GENOTYPING"

# --- CONFIGURATION DES RESSOURCES CPU PARALLÈLES ---
# Ex: 4 échantillons traités en même temps x 16 threads = 64 cœurs occupés
PARALLEL_JOBS=4
THREADS_PER_JOB=16

# --- DÉTECTION AUTOMATIQUE DU SAMPLE SHEET DANS FASTQ_PASS ---
CSV_FILE=$(find "$DATA_DIR" -maxdepth 1 \( -name "*.csv" -o -name "*Sample_Sheet*" -o -name "*VHB*" -o -name "*VHD*" \) -type f | head -n 1 || true)

if [ -z "$CSV_FILE" ] || [ ! -f "$CSV_FILE" ]; then
    echo "❌ Erreur : Aucun fichier CSV / Sample Sheet trouvé dans $DATA_DIR/."
    exit 1
fi

echo "📄 Sample Sheet détecté : $CSV_FILE"

# Détection du séparateur (virgule ou point-virgule)
SEP=";"
if head -n 5 "$CSV_FILE" | grep -q ","; then
    SEP=","
fi

# --- PARAMÈTRES VIRiONT (SEUILS DE CO-INFECTION) ---
MI_CUTOFF=30.0      # Seuil de co-infection clinique (30%)
MIN_READS=50        # Nombre minimal de reads absolu

# =====================================================================
# AUTOMATISATION DE LA GESTION DU SCRIPT R
# =====================================================================
R_SCRIPT="${SCRIPT_DIR}/generate_report.R"

if [ ! -f "$R_SCRIPT" ]; then
    echo "❌ Erreur : Le script R de rapport est introuvable ($R_SCRIPT)."
    exit 1
fi

if [ ! -x "$R_SCRIPT" ]; then
    chmod +x "$R_SCRIPT"
fi

# Création des dossiers principaux de sortie
mkdir -p "$SUPDATA_REFSEQ" "$SUPDATA_DB" "$TRIMMED_DIR" "$AMPLICON_DIR" "$GENOTYPING_DIR"

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

# -----------------------------------------------------------------
# 3. INITIALISATION DES FICHIERS DE TRAÇABILITÉ VIRiONT V2
# -----------------------------------------------------------------
SUMMARY_MULTIINF="${GENOTYPING_DIR}/SUMMARY_Multi_Infection.tsv"
FASTQ_CONTENT="${RESULTS_DIR}/fastq_content.txt"
PARAM_FILE="${RESULTS_DIR}/param_file"

# Fichier des paramètres du run
cat << EOF > "$PARAM_FILE"
######################
#### PARAMS USED #####
######################
data repository: $(dirname "$CSV_FILE")/
result repository: ${RESULTS_DIR}/
database used: ${TARGET_REFSEQ}
read minlength: ${min_length}
read maxlength: ${max_length}
quality filtering: 10
min coverage for consensus generation: 20
multi-infection cutoff: ${MI_CUTOFF}
min reads threshold: ${MIN_READS}
EOF

# Fichier d'analyse du contenu FASTQ
cat << EOF > "$FASTQ_CONTENT"
##########################
##### FASTQ ANALYSIS #####
##########################
barcode repository containing fastq/gz files and used for analysis:
EOF

tail -n +2 "$CSV_FILE" | while IFS="$SEP" read -r col1 col2 col3 col4 col5 col6 col7 col8 || [ -n "$col1" ]; do
    [ -z "$col2" ] && continue
    [[ "$col1" == *"Plate"* ]] && continue
    
    component=$(echo "$col3" | tr -d '\r\n ')
    [ -z "$component" ] && component="$col2"
    
    num_barcode=$(echo "$component" | grep -o '[0-9]\+' | head -n 1 || true)
    if [ -n "$num_barcode" ]; then
        formatted_barcode=$(printf "%02d" "$num_barcode")
        sample_name=$(echo "$col2" | tr -d '\r\n ')
        echo "barcode_${formatted_barcode}_${sample_name}" >> "$FASTQ_CONTENT"
    fi
done

cat << EOF >> "$FASTQ_CONTENT"
##########################
barcode repository containing other files than fastq/gz and ignored for analysis:
##########################
list of problematic files:
##########################
empty barcode repositories and ignored for analysis:
EOF

# Réinitialisation du résumé global
> "$SUMMARY_MULTIINF"

echo "====================================================================="
echo "   EXECUTION ÉTAPE 02 (PARALLÉLISÉ AVEC XARGS : $PARALLEL_JOBS ÉCHANTILLONS EN SIMULTANÉ)"
echo "   Fichier CSV    : $(basename "$CSV_FILE")"
echo "   Virus / Tech   : $virus_name | $tech_name"
echo "   Base Génotypes : $ref_filename"
echo "   Filtre Taille  : Min = $min_length pb | Max = $max_length pb"
echo "   Threads / Job  : $THREADS_PER_JOB"
echo "====================================================================="

# =====================================================================
# 4. FONCTION UNIFIÉE DE TRAITEMENT EXPORTÉE POUR XARGS
# =====================================================================
process_single_sample() {
    local sample_id="$1"
    local min_length="$2"
    local max_length="$3"
    local DEHOST_DIR="$4"
    local TRIMMED_DIR="$5"
    local AMPLICON_DIR="$6"
    local GENOTYPING_DIR="$7"
    local DB_PREFIX="$8"
    local MI_CUTOFF="$9"
    local MIN_READS="${10}"
    local R_SCRIPT="${11}"
    local THREADS_PER_JOB="${12}"

    DEHOSTED_FASTQ="${DEHOST_DIR}/${sample_id}_dehosted.fastq.gz"
    if [ ! -f "$DEHOSTED_FASTQ" ]; then
        DEHOSTED_FASTQ="${DEHOST_DIR}/${sample_id}_meta.fastq.gz"
    fi

    if [ ! -f "$DEHOSTED_FASTQ" ]; then
        echo "⚠️ Output de step_01 introuvable pour $sample_id. Sauté."
        return 0
    fi

    TRIMMED_FASTQ="${TRIMMED_DIR}/${sample_id}_trimmed.fastq.gz"
    SAMPLE_AMPLICON_DIR="${AMPLICON_DIR}/${sample_id}"
    SAMPLE_GENO_DIR="${GENOTYPING_DIR}/${sample_id}"
    mkdir -p "$SAMPLE_GENO_DIR"

    echo "---------------------------------------------------------------------"
    echo "▶️ [DEBUT] Traitement : $sample_id"
    echo "---------------------------------------------------------------------"

    # --- ÉTAPE 03 : Chopper ---
    if [ ! -s "$TRIMMED_FASTQ" ]; then
        gunzip -c "$DEHOSTED_FASTQ" | chopper -l "$min_length" --maxlength "$max_length" | gzip -c > "$TRIMMED_FASTQ"
    fi

    if [ ! -s "$TRIMMED_FASTQ" ]; then
        echo "⚠️ Aucun read conservé par chopper pour $sample_id."
        rm -f "$TRIMMED_FASTQ"
        return 0
    fi

    # --- ÉTAPE 04 : Amplicon_sorter ---
    mkdir -p "$SAMPLE_AMPLICON_DIR"
    CONSENSUS_FILES=$(find "$SAMPLE_AMPLICON_DIR" -maxdepth 2 \( -name "*.fasta" -o -name "*.fa" \) ! -name "*unique*" 2>/dev/null || true)

    if [ -z "$CONSENSUS_FILES" ]; then
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
            -np "$THREADS_PER_JOB" > /dev/null

        CONSENSUS_FILES=$(find "$SAMPLE_AMPLICON_DIR" -maxdepth 2 \( -name "*.fasta" -o -name "*.fa" \) ! -name "*unique*" 2>/dev/null || true)
    fi

    if [ -z "$CONSENSUS_FILES" ]; then
        echo "⚠️ Aucun consensus produit par amplicon_sorter pour $sample_id."
        rm -f "$TRIMMED_FASTQ"
        return 0
    fi

    # --- ÉTAPE 05 : BLASTn ---
    RAW_COUNT_TABLE="${SAMPLE_GENO_DIR}/${sample_id}_raw_cluster_counts.tsv"
    echo -e "sample\tcluster\tread_count\tgenotype\tpident\tlength\tbitscore" > "$RAW_COUNT_TABLE"

    for cons_file in $CONSENSUS_FILES; do
        cons_name=$(basename "$cons_file" | sed 's/\.[^.]*$//')
        header=$(head -n 1 "$cons_file")
        read_count=$(echo "$header" | grep -oP '\(\K[0-9]+(?=\))' || true)
        [ -z "$read_count" ] && read_count=$(grep -c "^>" "$cons_file" || echo "1")

        blast_line=$(blastn -query "$cons_file" \
                            -db "$DB_PREFIX" \
                            -outfmt "6 qseqid sseqid pident length evalue bitscore" \
                            -max_hsps 1 2>/dev/null | head -n 1 || true)

        if [ -n "$blast_line" ]; then
            genotype=$(echo "$blast_line" | awk '{print $2}')
            pident=$(echo "$blast_line" | awk '{print $3}')
            length=$(echo "$blast_line" | awk '{print $4}')
            bitscore=$(echo "$blast_line" | awk '{print $6}')

            echo -e "${sample_id}\t${cons_name}\t${read_count}\t${genotype}\t${pident}\t${length}\t${bitscore}" >> "$RAW_COUNT_TABLE"
        fi
    done

    # --- ÉTAPE 06 : Agrégation & Filtre Co-Infections ---
    VALIDATED_SUMMARY="${SAMPLE_GENO_DIR}/${sample_id}_validated_genotypes.tsv"
    echo -e "sample\tgenotype\tbest_cluster\ttotal_reads\tratio_percent\tpident\tlength\tstatus" > "$VALIDATED_SUMMARY"

    awk -F'\t' -v cutoff="$MI_CUTOFF" -v min_r="$MIN_READS" '
        NR > 1 {
            geno = $4
            reads = $3 + 0
            bit = $7 + 0
            sum_reads[geno] += reads
            if (bit > max_bit[geno]) {
                max_bit[geno] = bit
                best_cluster[geno] = $2
                best_pident[geno] = $5
                best_len[geno] = $6
            }
        }
        END {
            max_g_reads = 0
            for (g in sum_reads) {
                if (sum_reads[g] > max_g_reads) max_g_reads = sum_reads[g]
            }
            for (g in sum_reads) {
                ratio = (sum_reads[g] / max_g_reads) * 100
                status = "REJECTED"
                if (ratio >= cutoff && sum_reads[g] >= min_r) status = "VALIDATED"
                printf "%s\t%s\t%s\t%d\t%.1f%%\t%s\t%s\t%s\n", 
                       "'"$sample_id"'", g, best_cluster[g], sum_reads[g], ratio, best_pident[g], best_len[g], status
            }
        }
    ' "$RAW_COUNT_TABLE" | sort -k4,4nr >> "$VALIDATED_SUMMARY"

    # --- ÉTAPE 07 : Script R ---
    PDF_REPORT="${SAMPLE_GENO_DIR}/${sample_id}_genotype_report.pdf"
    Rscript "$R_SCRIPT" "$VALIDATED_SUMMARY" "$MI_CUTOFF" "$PDF_REPORT" &>/dev/null || true

    # --- ÉTAPE 08 : Export FASTA ---
    grep -w "VALIDATED" "$VALIDATED_SUMMARY" | while IFS=$'\t' read -r s_id geno cluster rcount ratio pident length status; do
        MASTER_FASTA="${SAMPLE_GENO_DIR}/${sample_id}_${geno}_master_consensus.fasta"
        target_file=$(find "$SAMPLE_AMPLICON_DIR" -type f \( -name "${cluster}*.fasta" -o -name "${cluster}*.fa" -o -name "*${cluster}*" \) ! -name "*unique*" | head -n 1 || true)
        if [ -n "$target_file" ] && [ -s "$target_file" ]; then
            cp "$target_file" "$MASTER_FASTA"
        else
            first_cons=$(echo "$CONSENSUS_FILES" | head -n 1)
            cp "$first_cons" "$MASTER_FASTA"
        fi
    done

    echo "✅ [FIN] Traitement terminé pour $sample_id"
}

export -f process_single_sample

# =====================================================================
# 5. EXECUTION EN PARALLÈLE AVEC XARGS (-P 4)
# =====================================================================
SAMPLE_LIST=$(mktemp)

tail -n +2 "$CSV_FILE" | while IFS="$SEP" read -r plate_pos sample component forward reverse flowcell kit_id primers || [ -n "$sample" ]; do
    sample=$(echo "$sample" | tr -d '\r\n ')
    component=$(echo "$component" | tr -d '\r\n ')
    plate_pos=$(echo "$plate_pos" | tr -d '\r\n ')

    [ -z "$sample" ] && continue
    [[ "$plate_pos" == *"Plate"* ]] && continue
    [[ "$sample" == *"Sample"* ]] && continue

    num_barcode=$(echo "$component" | grep -o '[0-9]\+' | head -n 1 || true)
    [ -z "$num_barcode" ] && continue
    
    formatted_barcode=$(printf "%02d" "$num_barcode")
    echo "barcode_${formatted_barcode}_${sample}" >> "$SAMPLE_LIST"
done

# Parallélisation native xargs
cat "$SAMPLE_LIST" | xargs -I {} -P "$PARALLEL_JOBS" bash -c 'process_single_sample "$@"' _ {} "$min_length" "$max_length" "$DEHOST_DIR" "$TRIMMED_DIR" "$AMPLICON_DIR" "$GENOTYPING_DIR" "$DB_PREFIX" "$MI_CUTOFF" "$MIN_READS" "$R_SCRIPT" "$THREADS_PER_JOB"

rm -f "$SAMPLE_LIST"

# =====================================================================
# 6. CONSOLIDATION DU FICHIER GLOBAL DE MULTI-INFECTION
# =====================================================================
echo "--> Aggregation du fichier global SUMMARY_Multi_Infection.tsv..."
find "$GENOTYPING_DIR" -name "*_validated_genotypes.tsv" | while read -r summary_file; do
    awk -F'\t' 'NR>1 {
        ratio = $5; gsub(/%/, "", ratio);
        print $1"\t"$2"\t"$4"\t"ratio
    }' "$summary_file" >> "$SUMMARY_MULTIINF"
done

echo "====================================================================="
echo " 🎉 Étape 02 terminée avec succès sur l'ensemble du run !"
echo "====================================================================="