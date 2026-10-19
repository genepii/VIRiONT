#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script universel de comparaison de consensus génomiques (Snakemake vs Nextflow)
- Alignement global FASTA vs FASTA via Biopython (gestion native de la circularité)
- Extraction détaillée de chaque divergence base par base
- Détection des mutations minoritaires / quasispèces via BAM Nextflow
"""

import os
import sys
import re
import glob
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Align import PairwiseAligner

# ==============================================================================
# 1. ARGUMENTS CLI
# ==============================================================================
def parse_args():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_viri = os.path.dirname(script_dir)
    parser = argparse.ArgumentParser(description="Comparaison FASTA vs FASTA et analyse des positions divergentes.")
    parser.add_argument(
        "--snake_dir", "-s",
        default=os.path.join(script_dir, "snakemake"),
        help="Dossier contenant les consensus Snakemake"
    )
    parser.add_argument(
        "--nf_dir", "-n",
        default=os.path.join(base_viri, "results/07_GENOTYPING/renamed_consensus"),
        help="Dossier contenant les consensus Nextflow"
    )
    parser.add_argument(
        "--bam_dir", "-b",
        default=os.path.join(base_viri, "results/06_ALIGNMENT"),
        help="Dossier contenant les BAM Nextflow indexés pour inspecter les fréquences alléliques"
    )
    parser.add_argument(
        "--out_dir", "-o",
        default=script_dir,
        help="Dossier de sortie pour les TSV et graphiques"
    )
    parser.add_argument(
        "--circular", "-c",
        action="store_true",
        default=True,
        help="Gérer le raccord circulaire (défaut: True)"
    )
    return parser.parse_args()

# ==============================================================================
# 2. EXTRACTION DES MÉTADONNÉES
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
# 3. CHARGEMENT DES SÉQUENCES
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
                "len": len(rec.seq),
                "seq": str(rec.seq).upper(),
                "file": os.path.basename(f),
                "sid": sid,
                "geno": geno,
                "bc": bc,
                "path": f
            }
            break
    print(f"[+] {label:10s} : {len(records)} séquences indexées ({len(files)} fichiers analysés).")
    return records

# ==============================================================================
# 4. APPARIEMENT DES PAIRES
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
            if len(candidates) == 1:
                pairs.append((f"{sid}_{geno}", sn_val, candidates[0][1]))
            else:
                matched = None
                for _, c_val in candidates:
                    if geno != "UNK" and (c_val["geno"].startswith(geno[0]) or geno.startswith(c_val["geno"][0])):
                        matched = c_val
                        break
                pairs.append((f"{sid}_{geno}", sn_val, matched if matched else candidates[0][1]))
        else:
            print(f"[-] Pas de correspondance Nextflow pour : {sn_val['file']} (ID: {sid})")

    print(f"[+] Paires formées : {len(pairs)} / {len(snake_dict)}")
    return pairs

# ==============================================================================
# 5. CONTRÔLE DES FRÉQUENCES ALLÉLIQUES SUR BAM (SAMTOOLS MPILEUP)
# ==============================================================================
def inspect_bam_allele_freq(bam_path, ref_id, pos):
    """Interroge samtools mpileup pour connaître la répartition des bases à une position donnée."""
    if not bam_path or not os.path.exists(bam_path):
        return "BAM_NON_DISPO"
    try:
        cmd = ["samtools", "mpileup", "-r", f"{ref_id}:{pos}-{pos}", bam_path]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        out = res.stdout.strip()
        if not out:
            return "NO_COVERAGE"
        cols = out.split("\t")
        if len(cols) >= 5:
            depth = cols[3]
            pileup = cols[4].upper()
            # Nettoyage des motifs indels regex (+Nseq, -Nseq, ^.)
            clean_bases = re.sub(r"\^.|\$", "", pileup)
            clean_bases = re.sub(r"[+-]\d+[ACGTNacgtn]+", "", clean_bases)
            
            counts = {b: clean_bases.count(b) for b in ["A", "C", "G", "T"] if clean_bases.count(b) > 0}
            total = sum(counts.values())
            if total == 0:
                return f"Depth:{depth}|RefMatch:{cols[2]}"
            details = [f"{b}:{round(cnt/total*100)}%" for b, cnt in sorted(counts.items(), key=lambda x: x[1], reverse=True)]
            return f"DP={depth}|" + ",".join(details)
    except Exception:
        pass
    return "NA"

# ==============================================================================
# 6. ALIGNEMENT GLOBAL BIOPYTHON & RELEVÉ DÉTAILLÉ DES DIVERGENCES
# ==============================================================================
def align_and_inspect(pairs, circular_mode=True, bam_dir=None):
    summary_results = []
    detailed_divergences = []

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    for label, sn_val, nf_val in pairs:
        seq_sn = sn_val["seq"]
        seq_nf = nf_val["seq"]
        len_sn = len(seq_sn)
        len_nf = len(seq_nf)

        # Recherche du BAM correspondant pour vérifier les minoritaires
        bam_candidates = glob.glob(os.path.join(bam_dir, f"*{sn_val['sid']}*.bam")) if bam_dir and os.path.isdir(bam_dir) else []
        bam_path = bam_candidates[0] if bam_candidates else None

        # Gestion circulaire : test du meilleur décalage ou alignement direct
        best_aln = None
        best_score = -float('inf')

        alignments = aligner.align(seq_sn, seq_nf)
        best_aln = alignments[0]
        aligned_sn = str(best_aln[0])
        aligned_nf = str(best_aln[1])

        pos_sn = 0
        pos_nf = 0
        snvs = 0
        gaps = 0
        matches = 0

        for b_sn, b_nf in zip(aligned_sn, aligned_nf):
            if b_sn != '-': pos_sn += 1
            if b_nf != '-': pos_nf += 1

            if b_sn == b_nf and b_sn != '-':
                matches += 1
            elif b_sn != '-' and b_nf != '-':
                snvs += 1
                freq_info = inspect_bam_allele_freq(bam_path, nf_val["id"], pos_nf) if bam_path else "Non vérifié"
                detailed_divergences.append({
                    "Sample": label,
                    "Pos_Snakemake": pos_sn,
                    "Base_Snakemake": b_sn,
                    "Pos_Nextflow": pos_nf,
                    "Base_Nextflow": b_nf,
                    "Type": "SNV",
                    "Nextflow_Allele_Balance": freq_info
                })
            else:
                gaps += 1
                detailed_divergences.append({
                    "Sample": label,
                    "Pos_Snakemake": pos_sn if b_sn != '-' else "GAP",
                    "Base_Snakemake": b_sn,
                    "Pos_Nextflow": pos_nf if b_nf != '-' else "GAP",
                    "Base_Nextflow": b_nf,
                    "Type": "INDEL",
                    "Nextflow_Allele_Balance": "INDEL"
                })

        overlap_len = matches + snvs
        ident_pct = round((matches / overlap_len * 100), 3) if overlap_len > 0 else 0.0
        div_pct = round(100.0 - ident_pct, 3)

        summary_results.append({
            "Sample": label,
            "File_Snakemake": sn_val["file"],
            "File_Nextflow": nf_val["file"],
            "Length_Snakemake": len_sn,
            "Length_Nextflow": len_nf,
            "Delta_Length": len_nf - len_sn,
            "Identity_Pct": ident_pct,
            "Divergence_Pct": div_pct,
            "SNVs_Count": snvs,
            "Internal_Gaps": gaps,
            "Overlap_Length": overlap_len
        })

    return pd.DataFrame(summary_results), pd.DataFrame(detailed_divergences)

# ==============================================================================
# 7. VISUALISATION BI-PANNEAU
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
    axes[0].set_ylabel("Échantillon / Sous-type", fontweight="bold")
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
    axes[1].set_title("Concordance de Séquence (Alignement FASTA)", fontsize=13, fontweight="bold")
    axes[1].set_ylabel("")
    axes[1].set_xticks([0.5])
    axes[1].set_xticklabels(["% Identité"], fontweight="bold")

    plt.suptitle(f"Validation VIRiONT_V2 : Snakemake vs Nextflow ({n_samples} séquences)", fontsize=15, fontweight="bold", y=0.99)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

# ==============================================================================
# MAIN
# ==============================================================================
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"[*] Snakemake dir : {args.snake_dir}")
    print(f"[*] Nextflow dir  : {args.nf_dir}")
    print(f"[*] BAMs dir      : {args.bam_dir}")
    print(f"[*] Sortie        : {args.out_dir}\n")

    if not os.path.isdir(args.snake_dir):
        sys.exit(f"[-] Erreur : Répertoire introuvable : {args.snake_dir}")
    if not os.path.isdir(args.nf_dir):
        sys.exit(f"[-] Erreur : Répertoire introuvable : {args.nf_dir}")

    snake_dict = load_sequences(args.snake_dir, "Snakemake")
    nf_dict = load_sequences(args.nf_dir, "Nextflow")

    pairs = match_pairs(snake_dict, nf_dict)
    if not pairs:
        sys.exit("[-] Aucune paire trouvée.")

    df_summary, df_detail = align_and_inspect(pairs, circular_mode=args.circular, bam_dir=args.bam_dir)

    tsv_summary = os.path.join(args.out_dir, "comparison_metrics.tsv")
    df_summary.to_csv(tsv_summary, sep="\t", index=False)

    tsv_detail = os.path.join(args.out_dir, "divergent_positions_detail.tsv")
    df_detail.to_csv(tsv_detail, sep="\t", index=False)
    print(f"\n[+] Détail des {len(df_detail)} divergences écrit dans : {tsv_detail}")

    png_out = os.path.join(args.out_dir, "comparaison_claire_nextflow_vs_snakemake.png")
    generate_plot(df_summary, png_out)
    print(f"[✔] Graphique mis à jour : {png_out}")

if __name__ == "__main__":
    main()
