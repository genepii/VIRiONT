#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), "configs")

import numpy as np
from itertools import combinations
import pandas as pd
from Bio import SeqIO, Align
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#########################################
##### Prise en charge des arguments #####
#########################################
parser = argparse.ArgumentParser(description='Génération de la matrice de distance pairwise')
parser.add_argument('--out_dir', default="./")
parser.add_argument('--file', required=True)

args = parser.parse_args()

out_folder = os.path.abspath(args.out_dir)
os.makedirs(out_folder, exist_ok=True)

clufile = args.file
clustername = os.path.splitext(os.path.basename(clufile))[0]

# Chargement et déduplication des noms de séquences
records = {}
suivlist = ["", "_bis", "_ter", "_qua", "_qui", "_sex", "_sep", "_oct", "_non"]

def iter_num(seqname):
    for sfx in suivlist:
        candidate = f"{seqname}{sfx}"
        if candidate not in records:
            return candidate
    return f"{seqname}_{len(records)}"

for seq in SeqIO.parse(clufile, "fasta"):
    records[iter_num(seq.id)] = seq.seq

n = len(records)
if n < 2:
    print(f"Moins de 2 séquences dans {clufile}, calcul impossible.")
    sys.exit(0)

# Moteur d'alignement global Biopython
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.wildcard = "N"
aligner.open_internal_gap_score = -15

matrix = np.full((n, n), np.nan)
noms = list(records.keys())

def matching_length(alignmt):
    start_seq1 = alignmt.aligned[0][0][0]
    start_seq2 = alignmt.aligned[1][0][0]
    end_seq1 = alignmt.aligned[0][-1][1]
    end_seq2 = alignmt.aligned[1][-1][1]
    return max(end_seq1 - start_seq1, end_seq2 - start_seq2)

csv_out_path = os.path.join(out_folder, f"{clustername}_pw-brief.csv")
with open(csv_out_path, "w") as f:
    f.write("strain1;strain2;diffperc;cluster;match_length\n")

for i, j in combinations(range(n), 2):
    s1 = records[noms[i]]
    s2 = records[noms[j]]
    s2_rc = s2.reverse_complement()

    score_fwd = aligner.score(s1, s2)
    score_rev = aligner.score(s1, s2_rc)

    if score_fwd >= score_rev:
        best_aln = aligner.align(s1, s2)[0]
    else:
        best_aln = aligner.align(s1, s2_rc)[0]

    try:
        mismatches = best_aln.counts().mismatches
        lensite = matching_length(best_aln)
        if lensite > 0:
            diff_pct = (mismatches / lensite) * 100.0
            matrix[j][i] = diff_pct

            with open(csv_out_path, "a") as f:
                f.write(f"{noms[i]};{noms[j]};{diff_pct:.2f};{clustername};{lensite}\n")
    except Exception as e:
        print(f"Erreur comparaison {noms[i]} vs {noms[j]}: {e}")

######################################################
##### Heatmap de la matrice de distance #####
######################################################
matrice = pd.DataFrame(matrix, index=noms, columns=noms)

# Dimensionnement dynamique
figsize = max(11, 7 + int(n * 0.7))
font_size = max(9, min(15, int(130 / n)))

fig, ax = plt.subplots(figsize=(figsize, figsize))

# 1. Tracé de la heatmap SANS annot seaborn (pour éviter le bug des NaN)
sns.heatmap(
    matrice,
    annot=False,
    linewidths=1.0,
    linecolor="white",
    square=True,
    cmap="YlOrRd",
    cbar=True,
    cbar_kws={'label': '% Divergence', 'shrink': 0.8},
    ax=ax
)

# Fond gris clair pour les NaN
ax.collections[0].cmap.set_bad('#e8e8e8')

# 2. Écriture manuelle et forcée du texte dans chaque case calculée
for row_idx in range(n):
    for col_idx in range(n):
        val = matrix[row_idx, col_idx]
        if not np.isnan(val):
            # Couleur du texte : blanc sur fond très rouge (> 5%), noir sinon
            text_color = "white" if val > 5.0 else "black"
            ax.text(
                col_idx + 0.5,
                row_idx + 0.5,
                f"{val:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=font_size,
                fontweight="bold"
            )

plt.xticks(rotation=45, ha='right', fontsize=10, fontweight='bold')
plt.yticks(rotation=0, fontsize=10, fontweight='bold')

plt.tight_layout()

png_out_path = os.path.join(out_folder, f"distmat_{clustername}.png")
plt.savefig(png_out_path, dpi=300, bbox_inches='tight')
plt.close()