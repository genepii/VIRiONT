#!/bin/bash

# =====================================================================
# CONFIGURATION DES CHEMINS & AUTONOMIE DES DROITS
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

chmod +x "$SCRIPT_DIR"/*.sh "$SCRIPT_DIR"/*.R 2>/dev/null || true

RESULTS_DIR="results_test"
GENOTYPING_DIR="${RESULTS_DIR}/05_GENOTYPING"
CONSENSUS_DIR="${RESULTS_DIR}/07_CONSENSUS"
PHYLO_DIR="${RESULTS_DIR}/09_PHYLOGENY"

REF_DATABASE=$(find "${RESULTS_DIR}/00_SUPDATA" -type f -name "*.fasta" 2>/dev/null | head -n 1)

ALL_CONS_FASTA="${CONSENSUS_DIR}/all_samples_consensus.fasta"
FILTERED_CONS_FASTA="${PHYLO_DIR}/filtered_consensus.fasta"
SELECTED_REFS_FASTA="${PHYLO_DIR}/selected_references.fasta"
ALL_SEQ_FASTA="${PHYLO_DIR}/all_sequences_for_tree.fasta"
ALIGNED_FASTA="${PHYLO_DIR}/aligned_sequences.fasta"
TREE_FILE="${PHYLO_DIR}/IQtree_analysis.treefile"
TREE_PDF="${PHYLO_DIR}/PHYLOGRAM_tree.pdf"

R_PLOT_TREE="${SCRIPT_DIR}/plot_tree.R"

mkdir -p "$PHYLO_DIR"

echo "====================================================================="
echo "  ÉTAPE 05 : PHYLOGÉNIE & DÉTECTION DE CONTAMINATION (IQ-TREE 3.1.3)"
echo "  Dossier de sortie : $PHYLO_DIR"
echo "====================================================================="

if [ ! -f "$ALL_CONS_FASTA" ]; then
    echo "❌ Erreur : Fichier consensus général introuvable ($ALL_CONS_FASTA)."
    exit 1
fi

# =====================================================================
# 1. FILTRAGE DES CONSENSUS RETENUS (Couverture >= 90%)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 1. Filtrage des consensus retenus (Exclusion si > 10% de N)..."

awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1;
    seq="";
    for(i=2;i<=NF;i++) seq=seq $i;
    gsub(/[ \t\r\n]/, "", seq);
    len=length(seq);
    n_count=gsub(/[Nn]/, "", seq);
    
    if (len > 0 && ((len - n_count) / len * 100) >= 90) {
        printf ">%s\n%s\n", header, seq;
    }
}' "$ALL_CONS_FASTA" > "$FILTERED_CONS_FASTA"

num_kept=$(grep -c "^>" "$FILTERED_CONS_FASTA" || echo 0)
echo "   ✅ $num_kept séquences consensus qualifiées conservées."

if [ "$num_kept" -eq 0 ]; then
    echo "⚠️ Aucun consensus n'a atteint le seuil de 90% de couverture."
    exit 0
fi

# =====================================================================
# 2. INGESTION DE TOUTES LES RÉFÉRENCES
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 2. Chargement de l'intégralité du panel de références..."

> "$SELECTED_REFS_FASTA"

if [ -n "$REF_DATABASE" ] && [ -f "$REF_DATABASE" ]; then
    cat "$REF_DATABASE" > "$SELECTED_REFS_FASTA"
fi

num_refs=$(grep -c "^>" "$SELECTED_REFS_FASTA" || echo 0)
echo "   ✅ $num_refs séquences de référence intégrées à l'analyse."


cat "$SELECTED_REFS_FASTA" "$FILTERED_CONS_FASTA" | awk '/^>/ {print $0; next} {gsub(/[^ATGCNatgcn-]/, "N"); print $0}' > "$ALL_SEQ_FASTA"

# =====================================================================
# 3. ALIGNEMENT MULTIPLE (MAFFT)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 3. Alignement Multiple via MAFFT..."

mafft --auto "$ALL_SEQ_FASTA" > "$ALIGNED_FASTA" 2>"${PHYLO_DIR}/mafft.log" || true

if [ -s "$ALIGNED_FASTA" ]; then
    echo "   ✅ Alignement MAFFT terminé."
else
    echo "❌ Erreur : L'alignement MAFFT est vide."
    cat "${PHYLO_DIR}/mafft.log"
    exit 1
fi

# =====================================================================
# 4. INFÉRENCE PHYLOGÉNÉTIQUE (IQ-TREE)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 4. Inférence de l'arbre via IQ-TREE..."

IQCMD="iqtree"
command -v iqtree3 &>/dev/null && IQCMD="iqtree3"

$IQCMD -s "$ALIGNED_FASTA" \
       -m MFP \
       -B 1000 \
       -alrt 1000 \
       -T AUTO \
       --prefix "${PHYLO_DIR}/IQtree_analysis" \
       -redo > "${PHYLO_DIR}/iqtree.log" 2>&1 || true

if [ -f "$TREE_FILE" ] && [ -s "$TREE_FILE" ]; then
    echo "   ✅ Arbre phylogénétique généré avec succès ($TREE_FILE)."
else
    echo "❌ Erreur IQ-TREE."
    cat "${PHYLO_DIR}/iqtree.log"
    exit 1
fi

# =====================================================================
# 5. RENDU GRAPHIQUE RECTANGULAIRE (R)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 5. Génération de l'arbre PDF rectangulaire..."

Rscript "$R_PLOT_TREE" "$TREE_FILE" "$TREE_PDF"

echo "====================================================================="
echo " 🎉 Étape 05 (Phylogénie & Arbre) terminée avec succès !"
echo " 📄 Arbre Rectangulaire PDF : $TREE_PDF"
echo "====================================================================="