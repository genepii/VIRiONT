#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script de comparaison de consensus génomiques (Snakemake vs Nextflow)
100 % Biopython natif via Bio.Align.PairwiseAligner.
Gère les génomes circulaires (duplication virtuelle) et le brin inverse.
"""

import os
import sys
import re
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# ==============================================================================
# 1. CLI ARGUMENTS
# ==============================================================================
def parse_args():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description="Comparaison de consensus avec Bio.Align (Biopython).")
    parser.add_argument(
        "--snake_dir", "-s",
        default=os.path.join(script_dir, "snakemake"),
        help="Dossier des consensus Snakemake"
    )
    parser.add_argument(
        "--nf_dir", "-n",
        default="/srv/scratch/chu-lyon.fr/alamiso/VIRiONT_PRONAME_REC/results/07_GENOTYPING/renamed_consensus",
        help="Dossier des consensus Nextflow"
    )
    parser.add_argument(
        "--out_dir", "-o",
        default=script_dir,
        help="Dossier de sortie (TSV et figure)"
    )
    parser.add_argument(
        "--circular", "-c",
        action="store_true",
        default=True,
        help="Duplication de séquence pour génomes circulaires (défaut: True)"
    )
    return parser.parse_args()

# ==============================================================================
# 2. MÉTADONNÉES
# ==============================================================================
def parse_metadata(filename):
    base = os.path.basename(filename)

    if "FJ349241" in base.upper(): return "26018776295", "C2", "BC95"
    if "377549" in base.upper():   return "26018776296", "B2", "BC96"
    if "BC11" in base.upper() or "RECOMB" in base.upper(): return "26018776297", "C2", "BC97"
    if "DQ089802" in base.upper(): return "26018776298", "C2", "BC98"

    m_bc = re.search(r"(?:barcode_?|BC)(\d{1,2})(?:_|\.|\b)", base, re.IGNORECASE)
    bc = f"BC{int(m_bc.group(1)):02d}" if m_bc else None

    m_num = re.search(r"(\d{10,12})", base)
    if m_num:
        sample_id = m_num.group(1).lstrip('0')
    else:
        clean_base = re.sub(r"^(?:barcode_?\d{1,2}_|BC\d{1,2}_)", "", base, flags=re.IGNORECASE)
        sample_id = clean_base.split('_')[0].split('.')[0]

    if bc == "BC97" or "26018776297" in base: sample_id = "26018776297"
    elif bc == "BC98" or "26018776298" in base: sample_id = "26018776298"

    m_geno_vhb = re.search(r"(?:_Geno_|_|-)([A-I]\d*)(?:_|\.|\s|$)", base)
    m_geno_vhd = re.search(r"(?:_Geno_|_|HDV)(\d[a-z]?)(?:_|\.|\s|$)", base, re.IGNORECASE)

    if m_geno_vhb:
        geno = m_geno_vhb.group(1).upper()
    elif m_geno_vhd:
        geno = m_geno_vhd.group(1).lower()
    else:
        m_iso = re.search(r"_([A-I])(?:_|\.|\b)", base)
        geno = m_iso.group(1).upper() if m_iso else "UNK"

    return sample_id, geno, bc

# ==============================================================================
# 3. CHARGEMENT
# ==============================================================================
def load_sequences(directory, label):
    records = {}
    files = glob.glob(os.path.join(directory, "*.fasta")) + \
            glob.glob(os.path.join(directory, "*.fa")) + \
            glob.glob(os.path.join(directory, "*.fna"))

    for f in sorted(files):
        sid, geno, bc = parse_metadata(f)
        for rec in SeqIO.parse(f, "fasta"):
            key = f"{sid}_{geno}"
            if key in records and bc:
                key = f"{sid}_{geno}_{bc}"
            records[key] = {
                "id": rec.id,
                "seq": str(rec.seq).upper(),
                "len": len(rec.seq),
                "file": os.path.basename(f),
                "sid": sid,
                "geno": geno,
                "bc": bc,
                "path": f
            }
            break
    print(f"[+] {label:10s} : {len(records)} séquences indexées.")
    return records

# ==============================================================================
# 4. APPARIEMENT
# ==============================================================================
def match_pairs(snake_dict, nf_dict):
    pairs = []
    nf_by_id = {}
    for k, val in nf_dict.items():
        nf_by_id.setdefault(val["sid"], []).append((k, val))

    for key, sn_val in sorted(snake_dict.items()):
        sid = sn_val["sid"]
        geno = sn_val["geno"]

        if key in nf_dict:
            pairs.append((key, sn_val, nf_dict[key]))
        elif sid in nf_by_id:
            candidates = nf_by_id[sid]
            matched = None
            for _, c_val in candidates:
                if geno != "UNK" and (c_val["geno"].startswith(geno[0]) or geno.startswith(c_val["geno"][0])):
                    matched = c_val
                    break
            pairs.append((f"{sid}_{geno}", sn_val, matched if matched else candidates[0][1]))
        else:
            print(f"[-] Pas de correspondance pour : {sn_val['file']}")

    print(f"\n[+] Paires concordantes : {len(pairs)} / {len(snake_dict)}")
    return pairs

# ==============================================================================
# 5. ALIGNEMENT VIA BIO.ALIGN.PAIRWISEALIGNER
# ==============================================================================
def run_alignment_biopython(pairs, circular_mode=True):
    aligner = PairwiseAligner()
    # Mode local (Smith-Waterman) pour permettre l'alignement sur cible doublée
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -2.0

    results = []

    for label, sn_val, nf_val in pairs:
        seq_sn = sn_val["seq"]
        seq_nf = nf_val["seq"]

        # Cible Nextflow doublée pour s'affranchir de la coupure d'origine circulaire
        target_seq = (seq_nf + seq_nf) if circular_mode else seq_nf

        best_alignment = None
        best_strand = "Forward"
        best_score = -float("inf")

        # Évaluation des deux brins
        from Bio.Seq import Seq
        seq_sn_obj = Seq(seq_sn)
        orientations = [
            ("Forward", str(seq_sn_obj)),
            ("Reverse-Complement", str(seq_sn_obj.reverse_complement()))
        ]

        for strand_name, query_seq in orientations:
            alignments = aligner.align(target_seq, query_seq)
            aln = next(iter(alignments), None)
            if aln and aln.score > best_score:
                best_score = aln.score
                best_alignment = aln
                best_strand = strand_name

        if best_alignment:
            counts = best_alignment.counts()
            matches = counts.identities
            mismatches = counts.mismatches
            gaps = counts.gaps
            aln_len = matches + mismatches + gaps

            # Identité calculée sur la longueur alignée (similaire à BLAST)
            pident = round((matches / (matches + mismatches + gaps)) * 100.0, 3) if aln_len > 0 else 0.0
            divergence = round(100.0 - pident, 3)
        else:
            pident, divergence, aln_len, mismatches, gaps, best_strand = 0.0, 100.0, 0, 0, 0, "Unknown"

        results.append({
            "Sample": label,
            "File_Snakemake": sn_val["file"],
            "File_Nextflow": nf_val["file"],
            "Length_Snakemake": sn_val["len"],
            "Length_Nextflow": nf_val["len"],
            "Delta_Length": nf_val["len"] - sn_val["len"],
            "Orientation": best_strand,
            "Identity_Pct": pident,
            "Divergence_Pct": divergence,
            "SNVs_Count": mismatches,
            "Internal_Gaps": gaps,
            "Overlap_Length": aln_len
        })

    return pd.DataFrame(results)

# ==============================================================================
# 6. VISUALISATION BI-PANNEAU
# ==============================================================================
def generate_plot(df, out_png):
    n_samples = len(df)
    fig_height = max(8, n_samples * 0.55)
    fig, axes = plt.subplots(1, 2, figsize=(22, fig_height), gridspec_kw={'width_ratios': [1.1, 1.4]})

    df_melt = df.melt(
        id_vars=['Sample'],
        value_vars=['Length_Snakemake', 'Length_Nextflow'],
        var_name='Pipeline',
        value_name='Taille_bp'
    )
    df_melt['Pipeline'] = df_melt['Pipeline'].replace({
        'Length_Snakemake': 'Snakemake (v1)',
        'Length_Nextflow': 'Nextflow (VIRiONT_V2)'
    })

    sns.barplot(
        data=df_melt,
        y='Sample',
        x='Taille_bp',
        hue='Pipeline',
        palette=['#5DADE2', '#1B4F72'],
        ax=axes[0]
    )
    axes[0].set_title("Comparaison des Longueurs (bp)", fontsize=13, fontweight="bold")
    axes[0].set_xlabel("Taille du consensus (bp)", fontweight="bold")
    axes[0].set_ylabel("Échantillon / Génotype", fontweight="bold")
    axes[0].grid(axis='x', linestyle='--', alpha=0.6)
    axes[0].legend(title="Pipeline", loc="lower right")

    df_annot = df.set_index("Sample")
    sns.heatmap(
        df_annot[["Identity_Pct"]],
        annot=True,
        fmt=".2f",
        cmap="Blues",
        vmin=90,
        vmax=100,
        cbar=True,
        cbar_kws={'label': '% Identité', 'shrink': 0.8, 'pad': 0.38},
        ax=axes[1]
    )

    for i, sample in enumerate(df["Sample"]):
        div = df.loc[df["Sample"] == sample, "Divergence_Pct"].values[0]
        if div < 1.0:
            statut = "(Quasi-identique <1%)"
        elif div < 4.0:
            statut = "(Intra-sous-génotype <4%)"
        elif div < 7.5:
            statut = "(Divergence sous-type 4-7.5%)"
        else:
            statut = "(Divergence génotype >7.5%)"
        axes[1].text(1.06, i + 0.5, f"Div: {div:.2f}%  {statut}", va='center', ha='left', fontsize=9.5, fontweight="bold")

    axes[1].set_xlim(0, 3.4)
    axes[1].set_title("Concordance de Séquence (PairwiseAligner)", fontsize=13, fontweight="bold")
    axes[1].set_ylabel("")
    axes[1].set_xticks([0.5])
    axes[1].set_xticklabels(["% Identité"], fontweight="bold")

    plt.suptitle(f"Validation VIRiONT_V2 : Biopython PairwiseAligner ({n_samples} séquences)", fontsize=15, fontweight="bold", y=0.99)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

# ==============================================================================
# MAIN
# ==============================================================================
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if not os.path.isdir(args.snake_dir) or not os.path.isdir(args.nf_dir):
        sys.exit("[-] Erreur : répertoire d'entrée introuvable.")

    snake_dict = load_sequences(args.snake_dir, "Snakemake")
    nf_dict = load_sequences(args.nf_dir, "Nextflow")

    pairs = match_pairs(snake_dict, nf_dict)
    if not pairs:
        sys.exit("[-] Aucune paire correspondante trouvée.")

    df = run_alignment_biopython(pairs, circular_mode=args.circular)

    tsv_out = os.path.join(args.out_dir, "comparison_metrics_biopython.tsv")
    df.to_csv(tsv_out, sep="\t", index=False)
    print(f"\n[+] Tableau TSV écrit : {tsv_out}\n")
    print(df[["Sample", "Length_Snakemake", "Length_Nextflow", "Identity_Pct", "Divergence_Pct", "SNVs_Count", "Internal_Gaps", "Overlap_Length"]].to_string(index=False))

    png_out = os.path.join(args.out_dir, "comparaison_biopython_nextflow_vs_snakemake.png")
    generate_plot(df, png_out)
    print(f"\n[✔] Graphique généré avec succès : {png_out}")

if __name__ == "__main__":
    main()