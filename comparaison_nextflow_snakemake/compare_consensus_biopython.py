#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# Répertoires
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SNAKE_DIR = os.path.join(SCRIPT_DIR, "snakemake")
NF_MULTI_FASTA = os.path.join(SCRIPT_DIR, "nextflow/validated_consensus_all.fasta")
OUT_DIR = SCRIPT_DIR

IUPAC_BASES = set("RYSWKMBDHVN")

def extract_meta_snakemake(filename):
    base = os.path.basename(filename)
    if "FJ349241" in base: return "FJ349241", "C2"
    if "377549" in base:   return "QC377549", "B2"
    if "BC11" in base:     return "BC11", "C2"
    if "DQ089802" in base: return "DQ089802", "C2"

    m_id = re.search(r"(\d{9,12})", base)
    sid = m_id.group(1).lstrip('0') if m_id else base.split('_')[0]

    m_geno = re.search(r"_([A-I]\d?)(?:_cons|\.fasta|$)", base)
    geno = m_geno.group(1) if m_geno else "UNK"
    return sid, geno

def extract_meta_nextflow(header):
    tokens = header.split()[0].split('_')
    
    if "FJ349241" in header: return "FJ349241", "C2", header.split()[0]
    if "377549" in header:   return "QC377549", "B2", header.split()[0]
    if "26000000097" in header or "BC11" in header: return "BC11", "C2", header.split()[0]
    if "DQ089802" in header or "REF_UA2" in header: return "DQ089802", "C2", header.split()[0]

    m_id = re.search(r"(\d{9,12})", header)
    sid = m_id.group(1).lstrip('0') if m_id else tokens[2]

    m_geno = re.search(r"Geno_([A-I]\d?)", header)
    geno = m_geno.group(1) if m_geno else "UNK"
    return sid, geno, tokens[0] + "_" + tokens[1]

# 1. Charger les FASTA Snakemake
snake_records = {}
for f in glob.glob(os.path.join(SNAKE_DIR, "*.fa*")):
    sid, geno = extract_meta_snakemake(f)
    for rec in SeqIO.parse(f, "fasta"):
        snake_records[sid] = {
            "sid": sid,
            "geno": geno,
            "seq": str(rec.seq).upper(),
            "len": len(rec.seq),
            "file": os.path.basename(f)
        }
        break

# 2. Charger le multi-FASTA Nextflow
nf_records = {}
for rec in SeqIO.parse(NF_MULTI_FASTA, "fasta"):
    sid, geno, tag = extract_meta_nextflow(rec.id)
    nf_records[sid] = {
        "sid": sid,
        "geno": geno,
        "seq": str(rec.seq).upper(),
        "len": len(rec.seq),
        "tag": tag,
        "header": rec.id
    }

print(f"[+] Snakemake : {len(snake_records)} séquences indexées.")
print(f"[+] Nextflow  : {len(nf_records)} séquences indexées.")

# 3. Association des paires
pairs = []
for sid, sn_data in snake_records.items():
    if sid in nf_records:
        pairs.append((sid, sn_data, nf_records[sid]))
    else:
        found = False
        for nf_sid, nf_data in nf_records.items():
            if sid in nf_sid or nf_sid in sid:
                pairs.append((sid, sn_data, nf_data))
                found = True
                break
        if not found:
            print(f"[-] Pas de match Nextflow pour {sid} ({sn_data['file']})")

print(f"[+] Paires appariées : {len(pairs)}\n")

# 4. Alignement Pairwise FASTA vs FASTA (filtrage exclusif SNVs)
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -1

summary_list = []
snv_divergences_list = []

for sid, sn_data, nf_data in sorted(pairs, key=lambda x: x[0]):
    seq_sn = sn_data["seq"]
    seq_nf = nf_data["seq"]
    geno = sn_data["geno"]

    aln = aligner.align(seq_sn, seq_nf)[0]
    aln_sn, aln_nf = str(aln[0]), str(aln[1])

    pos_sn = 0
    pos_nf = 0
    matches = 0
    snvs = 0

    for b_sn, b_nf in zip(aln_sn, aln_nf):
        if b_sn != '-': pos_sn += 1
        if b_nf != '-': pos_nf += 1

        if b_sn == b_nf and b_sn != '-':
            matches += 1
        elif b_sn != '-' and b_nf != '-':
            snvs += 1
            snv_divergences_list.append({
                "Sample": sid,
                "Genotype": geno,
                "Pos_Snakemake": pos_sn,
                "Base_Snakemake": b_sn,
                "Pos_Nextflow": pos_nf,
                "Base_Nextflow": b_nf,
                "Type": "SNV",
                "Quasispecies_IUPAC": "OUI" if b_sn in IUPAC_BASES else "NON"
            })

    total_sites = matches + snvs
    ident = round((matches / total_sites * 100), 2) if total_sites > 0 else 0.0
    div = round(100.0 - ident, 2)

    summary_list.append({
        "Sample": f"{sid}_{geno}",
        "Length_Snakemake": sn_data["len"],
        "Length_Nextflow": nf_data["len"],
        "Identity_Pct": ident,
        "Divergence_Pct": div,
        "SNVs": snvs
    })

df_summary = pd.DataFrame(summary_list)
df_snvs = pd.DataFrame(snv_divergences_list)

summary_tsv = os.path.join(OUT_DIR, "summary_concordance_biopython.tsv")
detail_tsv = os.path.join(OUT_DIR, "divergent_positions_detail.tsv")

df_summary.to_csv(summary_tsv, sep="\t", index=False)
df_snvs.to_csv(detail_tsv, sep="\t", index=False)

print(df_summary[["Sample", "Length_Snakemake", "Length_Nextflow", "Identity_Pct", "Divergence_Pct", "SNVs"]].to_string(index=False))
print(f"\n[✔] Résumé écrit dans : {summary_tsv}")
print(f"[✔] Détail des SNVs écrit dans : {detail_tsv} ({len(df_snvs)} SNVs relevés)")

# 5. Graphique de concordance sans indels ni chevauchement
fig, axes = plt.subplots(1, 2, figsize=(22, max(8, len(df_summary) * 0.55)), gridspec_kw={'width_ratios': [1.1, 1.4]})

df_melt = df_summary.melt(id_vars=['Sample'], value_vars=['Length_Snakemake', 'Length_Nextflow'], var_name='Pipeline', value_name='Taille_bp')
df_melt['Pipeline'] = df_melt['Pipeline'].replace({'Length_Snakemake': 'Snakemake (v1)', 'Length_Nextflow': 'Nextflow (VIRiONT_V2)'})

sns.barplot(data=df_melt, y='Sample', x='Taille_bp', hue='Pipeline', palette=['#5DADE2', '#1B4F72'], ax=axes[0])
axes[0].set_title("Tailles des consensus FASTA (bp)", fontsize=13, fontweight="bold")
axes[0].set_xlabel("Taille (bp)", fontweight="bold")
axes[0].set_ylabel("Échantillon / Génotype", fontweight="bold")
axes[0].grid(axis="x", linestyle="--", alpha=0.6)
axes[0].legend(loc="lower right")

df_annot = df_summary.set_index("Sample")
sns.heatmap(
    df_annot[["Identity_Pct"]], 
    annot=True, 
    fmt=".2f", 
    cmap="Blues", 
    vmin=95, 
    vmax=100, 
    cbar=True, 
    cbar_kws={'label': '% Identité', 'shrink': 0.8, 'pad': 0.35},
    ax=axes[1]
)

for i, sample in enumerate(df_summary["Sample"]):
    div_val = df_summary.loc[df_summary["Sample"] == sample, "Divergence_Pct"].values[0]
    snv_val = df_summary.loc[df_summary["Sample"] == sample, "SNVs"].values[0]
    axes[1].text(
        1.05, 
        i + 0.5, 
        f"Div: {div_val:.2f}%  [SNVs: {snv_val}]", 
        va='center', 
        ha='left', 
        fontsize=10, 
        fontweight="bold"
    )

axes[1].set_xlim(0, 3.4)
axes[1].set_title("Identité nucléotidique FASTA (SNVs uniquement)", fontsize=13, fontweight="bold")
axes[1].set_ylabel("")
axes[1].set_xticks([0.5])
axes[1].set_xticklabels(["% Identité"], fontweight="bold")

plt.suptitle(f"Validation VIRiONT_V2 : Snakemake vs Nextflow ({len(df_summary)} séquences)", fontsize=15, fontweight="bold", y=0.99)
plt.tight_layout()

png_out = os.path.join(OUT_DIR, "comparaison_biopython_fasta.png")
plt.savefig(png_out, dpi=300, bbox_inches='tight')
plt.close()
print(f"[✔] Graphique généré (SNVs uniquement) : {png_out}")
