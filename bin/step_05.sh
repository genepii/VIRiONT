#!/bin/bash

# Interruption en cas d'erreur globale bloquante
set -eo pipefail

# =====================================================================
# CONFIGURATION DES CHEMINS & RESSOURCES (VIRiONT V2)
# =====================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

chmod +x "$SCRIPT_DIR"/*.sh "$SCRIPT_DIR"/*.R 2>/dev/null || true

RESULTS_DIR="results_test"
GENOTYPING_DIR="${RESULTS_DIR}/05_GENOTYPING"
CONSENSUS_DIR="${RESULTS_DIR}/07_CONSENSUS"
PHYLO_DIR="${RESULTS_DIR}/09_PHYLOGENY"

# --- ALLOCATION DYNAMIQUE DES RESSOURCES MULTI-THREADS ---
TOTAL_THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 16)
THREADS=$TOTAL_THREADS

REF_DATABASE=$(find "${RESULTS_DIR}/00_SUPDATA" -type f -name "*.fasta" 2>/dev/null | head -n 1)

ALL_CONS_FASTA="${CONSENSUS_DIR}/all_samples_consensus.fasta"
FILTERED_CONS_FASTA="${PHYLO_DIR}/filtered_consensus.fasta"
SELECTED_REFS_FASTA="${PHYLO_DIR}/selected_references.fasta"
ALL_SEQ_FASTA="${PHYLO_DIR}/all_sequences_for_tree.fasta"
RAW_ALIGNED_FASTA="${PHYLO_DIR}/raw_aligned.fasta"
ALIGNED_FASTA="${PHYLO_DIR}/aligned_sequences.fasta"
TREE_FILE="${PHYLO_DIR}/IQtree_analysis.treefile"
TREE_PDF="${PHYLO_DIR}/PHYLOGRAM_tree.pdf"

R_PLOT_TREE="${SCRIPT_DIR}/plot_tree.R"

mkdir -p "$PHYLO_DIR"

echo "====================================================================="
echo "   ÉTAPE 05 : PHYLOGÉNIE RIGOUREUSE HAUTE PRÉCISION (IQ-TREE & MAFFT)"
echo "   Dossier de sortie : $PHYLO_DIR"
echo "   Ressources CPU exploitées : $THREADS Threads"
echo "====================================================================="

if [ ! -f "$ALL_CONS_FASTA" ]; then
    echo "❌ Erreur : Fichier consensus général introuvable ($ALL_CONS_FASTA)."
    exit 1
fi

# =====================================================================
# 1. FILTRAGE ET PARSING DES CONSENSUS DU RUN (Seuil >= 50%)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 1. Filtrage des consensus retenus du run..."

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
{ seq = seq $0; }
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
echo "    ✅ $num_kept séquences consensus qualifiées conservées."

if [ "$num_kept" -eq 0 ]; then
    echo "⚠️ Aucun consensus n'a atteint le seuil de couverture."
    exit 0
fi

# =====================================================================
# 2. SELECTION DES RÉFÉRENCES (PANEL COMPLET & EQUILIBRE)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 2. Sélection des références..."

DETECTED_LETTRES=$(grep "^>" "$FILTERED_CONS_FASTA" | grep -oP '[A-H](?=[0-9]*)' | sort -u | tr '\n' ' ' || true)
echo "    🎯 Génotypes majeurs détectés dans le run : $DETECTED_LETTRES"

> "$SELECTED_REFS_FASTA"

if [ -n "$REF_DATABASE" ] && [ -f "$REF_DATABASE" ]; then
    awk -v det_g="$DETECTED_LETTRES" '
    BEGIN {
        split(det_g, dg, " ");
        for (i in dg) focus[dg[i]] = 1;
    }
    /^>/ {
        header = $0;
        match($0, /[A-H]/);
        geno_lettre = (RLENGTH > 0) ? substr($0, RSTART, 1) : "OTHER";
        
        # Sélection des références pour un arbre équilibré
        if (focus[geno_lettre] == 1 || count[geno_lettre] < 3) {
            keep = 1;
            count[geno_lettre]++;
            print header;
        } else {
            keep = 0;
        }
        next;
    }
    {
        if (keep == 1) print $0;
    }
    ' "$REF_DATABASE" > "$SELECTED_REFS_FASTA"
fi

num_refs=$(grep -c "^>" "$SELECTED_REFS_FASTA" || echo 0)
echo "    ✅ $num_refs séquences de référence sélectionnées."

# Fusion des séquences
(cat "$SELECTED_REFS_FASTA"; echo ""; cat "$FILTERED_CONS_FASTA") | \
awk '/^>/ {print $0; next} {gsub(/[^ATGCNatgcn-]/, "N"); print $0}' | \
grep -v '^$' > "$ALL_SEQ_FASTA"

# =====================================================================
# 3. ALIGNEMENT MULTIPLE AVEC MAFFT (MAXIMALE RIGUEUR)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 3. Alignement & Réorientation globale des brins (MAFFT multi-threads)..."

mafft --auto --thread "$THREADS" --adjustdirection "$ALL_SEQ_FASTA" > "$RAW_ALIGNED_FASTA" 2>"${PHYLO_DIR}/mafft.log" || true

sed -E 's/>_R_/>/g; s/>_R/>/g' "$RAW_ALIGNED_FASTA" > "$ALIGNED_FASTA"
rm -f "$RAW_ALIGNED_FASTA"

if [ -s "$ALIGNED_FASTA" ]; then
    echo "    ✅ Alignement terminé avec succès."
else
    echo "❌ Erreur alignement MAFFT."
    cat "${PHYLO_DIR}/mafft.log"
    exit 1
fi

# =====================================================================
# 4. INFÉRENCE PHYLOGÉNÉTIQUE (MODELFINDER RIGOUREUX -m MFP)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 4. Inférence de l'arbre via IQ-TREE (ModelFinder complet)..."

IQCMD="iqtree"
command -v iqtree3 &>/dev/null && IQCMD="iqtree3"

# -m MFP : Teste TOUS les modèles possibles pour trouver la vraie matrice de substitution évolutive.
# -B 1000 -alrt 1000 : Double validation statistique (Ultrafast Bootstrap + SH-aLRT)
$IQCMD -s "$ALIGNED_FASTA" \
       -m MFP \
       -B 1000 \
       -alrt 1000 \
       -T AUTO \
       -ntmax "$THREADS" \
       --prefix "${PHYLO_DIR}/IQtree_analysis" \
       -redo > "${PHYLO_DIR}/iqtree.log" 2>&1 || true

if [ -f "$TREE_FILE" ] && [ -s "$TREE_FILE" ]; then
    echo "    ✅ Arbre phylogénétique haute précision généré ($TREE_FILE)."
else
    echo "❌ Erreur IQ-TREE."
    cat "${PHYLO_DIR}/iqtree.log"
    exit 1
fi

# =====================================================================
# 5. RENDU GRAPHIQUE (R)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "--> 5. Génération du PDF..."

Rscript "$R_PLOT_TREE" "$TREE_FILE" "$TREE_PDF" &>/dev/null || true

echo "====================================================================="
echo " 🎉 Étape 05 terminée avec succès !"
echo " 📄 Arbre Rectangulaire PDF : $TREE_PDF"
echo "====================================================================="