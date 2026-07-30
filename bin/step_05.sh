#!/bin/bash

# =====================================================================
# CONFIGURATION DES CHEMINS & AUTONOMIE DES DROITS (VIRiONT V2)
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Auto-attribution des droits d'exécution sur les scripts bin/
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
RAW_ALIGNED_FASTA="${PHYLO_DIR}/raw_aligned.fasta"
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
# 1. FILTRAGE ET PARSING ROBUSTE DES CONSENSUS RETENUS (Seuil >= 50%)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 1. Filtrage et nettoyage des consensus retenus..."

awk '
/^>/ {
    if (NR > 1) {
        gsub(/[ \t\r\n]/, "", seq);
        len = length(seq);
        n_count = gsub(/[Nn]/, "", seq);
        if (len > 0 && ((len - n_count) / len * 100) >= 50) {
            print header "\n" seq;
        }
    }
    header = $0;
    seq = "";
    next;
}
{
    seq = seq $0;
}
END {
    if (header != "") {
        gsub(/[ \t\r\n]/, "", seq);
        len = length(seq);
        n_count = gsub(/[Nn]/, "", seq);
        if (len > 0 && ((len - n_count) / len * 100) >= 50) {
            print header "\n" seq;
        }
    }
}' "$ALL_CONS_FASTA" > "$FILTERED_CONS_FASTA"

num_kept=$(grep -c "^>" "$FILTERED_CONS_FASTA" || echo 0)
echo "   ✅ $num_kept séquences consensus qualifiées conservées."

if [ "$num_kept" -eq 0 ]; then
    echo "⚠️ Aucun consensus n'a atteint le seuil de couverture."
    exit 0
fi

# =====================================================================
# 2. FUSION SÉCURISÉE AVEC SAUT DE LIGNE GARANTI
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 2. Chargement du panel de références et fusion..."

> "$SELECTED_REFS_FASTA"

if [ -n "$REF_DATABASE" ] && [ -f "$REF_DATABASE" ]; then
    cat "$REF_DATABASE" > "$SELECTED_REFS_FASTA"
fi

num_refs=$(grep -c "^>" "$SELECTED_REFS_FASTA" || echo 0)
echo "   ✅ $num_refs séquences de référence intégrées à l'analyse."

# Le 'echo ""' empêche la fusion de la dernière ligne de référence avec le premier header consensus
(cat "$SELECTED_REFS_FASTA"; echo ""; cat "$FILTERED_CONS_FASTA") | \
awk '/^>/ {print $0; next} {gsub(/[^ATGCNatgcn-]/, "N"); print $0}' | \
grep -v '^$' > "$ALL_SEQ_FASTA"

# =====================================================================
# 3. ALIGNEMENT MULTIPLE & RÉORIENTATION GLOBALE DE TOUTES LES SÉQUENCES
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 3. Réorientation globale et Alignement Multiple via MAFFT..."

# 1. Alignement avec détection et réorientation automatique du sens (+ / -)
mafft --auto --adjustdirection "$ALL_SEQ_FASTA" > "$RAW_ALIGNED_FASTA" 2>"${PHYLO_DIR}/mafft.log" || true

# 2. Nettoyage strict des headers FASTA (suppression de _R_ sans altérer la séquence inversée)
awk '/^>/ {gsub(/^>_R_/, ">"); gsub(/^>_R/, ">"); print; next} {print}' "$RAW_ALIGNED_FASTA" > "$ALIGNED_FASTA"

rm -f "$RAW_ALIGNED_FASTA"

if [ -s "$ALIGNED_FASTA" ]; then
    echo "   ✅ Toutes les séquences ont été réorientées dans le bon sens et alignées."
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