#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script universel et robuste de comparaison de consensus génomiques (ex: Snakemake vs Nextflow)
Conçu pour gérer nativement les génomes circulaires (VHB, VHD, plasmides) sans perte de bases.
"""

import os
import sys
import re
import glob
import argparse
import subprocess
import tempfile
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# ==============================================================================
# 1. PARSING DES ARGUMENTS CLI
# ==============================================================================
def parse_args():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description="Comparaison robuste et universelle de deux pipelines génomiques.")
    parser.add_argument(
        "--snake_dir", "-s",
        default=os.path.join(script_dir, "snakemake"),
        help="Dossier contenant les consensus du pipeline 1 (Snakemake)"
    )
    parser.add_argument(
        "--nf_dir", "-n",
        default="/srv/scratch/chu-lyon.fr/alamiso/VIRiONT_PRONAME_REC/results/07_GENOTYPING/renamed_consensus",
        help="Dossier contenant les consensus du pipeline 2 (Nextflow)"
    )
    parser.add_argument(
        "--out_dir", "-o",
        default=script_dir,
        help="Dossier de sortie pour le TSV et le graphique"
    )
    parser.add_argument(
        "--circular", "-c",
        action="store_true",
        default=True,
        help="Activer la duplication virtuelle systématique pour génomes circulaires (défaut: True)"
    )
    return parser.parse_args()

# ==============================================================================
# 2. EXTRACTION ROBUSTE DES MÉTADONNÉES
# ==============================================================================
def parse_metadata(filename):
    base = os.path.basename(filename)

    # 1. Contrôles spécifiques reconnus par ID de référence
    if "FJ349241" in base.upper(): return "26018776295", "C2", "BC95"
    if "377549" in base.upper():   return "26018776296", "B2", "BC96"
    if "BC11" in base.upper() or "RECOMB" in base.upper(): return "26018776297", "C2", "BC97"
    if "DQ089802" in base.upper(): return "26018776298", "C2", "BC98"

    # 2. Barcode (ex: barcode_01, barcode83, BC02)
    m_bc = re.search(r"(?:barcode_?|BC)(\d{1,2})(?:_|\.|\b)", base, re.IGNORECASE)
    bc = f"BC{int(m_bc.group(1)):02d}" if m_bc else None

    # 3. Patient ID (numéro long prioritaire, sinon premier token propre)
    m_num = re.search(r"(\d{10,12})", base)
    if m_num:
        sample_id = m_num.group(1).lstrip('0')
    else:
        clean_base = re.sub(r"^(?:barcode_?\d{1,2}_|BC\d{1,2}_)", "", base, flags=re.IGNORECASE)
        sample_id = clean_base.split('_')[0].split('.')[0]

    # Barcodes contrôles Nextflow
    if bc == "BC97" or "26018776297" in base: sample_id = "26018776297"
    elif bc == "BC98" or "26018776298" in base: sample_id = "26018776298"

    # 4. Sous-type / Génotype
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
# 3. CHARGEMENT ET INDEXATION DES SÉQUENCES
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
                "file": os.path.basename(f),
                "sid": sid,
                "geno": geno,
                "bc": bc,
                "path": f
            }
            break
    print(f"[+] {label:10s} : {len(records)} séquences indexées depuis {len(files)} fichiers.")
    return records

# ==============================================================================
# 4. APPARIEMENT ROBUSTE DES ÉCHANTILLONS
# ==============================================================================
def match_pairs(snake_dict, nf_dict):
    pairs = []
    nf_by_id = {}
    for k, val in nf_dict.items():
        nf_by_id.setdefault(val["sid"], []).append((k, val))

    for key, sn_val in sorted(snake_dict.items()):
        sid = sn_val["sid"]
        geno = sn_val["geno"]

        # Match direct clé complète
        if key in nf_dict:
            pairs.append((key, sn_val, nf_dict[key]))
        # Match par identifiant patient
        elif sid in nf_by_id:
            candidates = nf_by_id[sid]
            if len(candidates) == 1:
                pairs.append((f"{sid}_{geno}", sn_val, candidates[0][1]))
            else:
                # Filtrage par compatibilité de sous-type si ambiguïté
                matched = None
                for _, c_val in candidates:
                    if geno != "UNK" and (c_val["geno"].startswith(geno[0]) or geno.startswith(c_val["geno"][0])):
                        matched = c_val
                        break
                if matched:
                    pairs.append((f"{sid}_{geno}", sn_val, matched))
                else:
                    # Repli sur le premier candidat si génotype non discriminant
                    pairs.append((f"{sid}_{geno}", sn_val, candidates[0][1]))
        else:
            print(f"[-] Pas de correspondance Nextflow trouvée pour : {sn_val['file']} (ID: {sid})")

    print(f"\n[+] Paires concordantes formées : {len(pairs)} / {len(snake_dict)}")
    return pairs

# ==============================================================================
# 5. COMPARAISON BLASTN UNIVERSELLE AVEC DUPLICATION CIRCULAIRE
# ==============================================================================
def run_alignment(pairs, circular_mode=True):
    results = []
    for label, sn_val, nf_val in pairs:
        path_sn = sn_val["path"]
        path_nf = nf_val["path"]

        target_path = path_nf
        tmp_target = None

        try:
            # Gestion systématique de la circularité :
            # En doublant la séquence cible, BLAST traverse le point de raccord sans aucune interruption
            if circular_mode:
                rec_nf = next(SeqIO.parse(path_nf, "fasta"))
                doubled_seq = rec_nf.seq + rec_nf.seq
                rec_doubled = rec_nf[:]
                rec_doubled.seq = doubled_seq
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tf:
                    tmp_target = tf.name
                    SeqIO.write(rec_doubled, tmp_target, "fasta")
                target_path = tmp_target

            cmd = [
                "blastn",
                "-query", path_sn,
                "-subject", target_path,
                "-outfmt", "6 pident length mismatch gaps sstrand qlen",
                "-dust", "no"
            ]

            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            lines = [l.strip() for l in proc.stdout.strip().split("\n") if l.strip()]

            if lines:
                hsps = []
                for l in lines:
                    cols = l.split("\t")
                    aln_len = int(cols[1])
                    qlen = int(cols[5])
                    # On s'assure qu'un alignement sur cible doublée ne dépasse pas la taille physique réelle
                    bounded_len = min(aln_len, qlen)
                    hsps.append({
                        "pident": float(cols[0]),
                        "length": bounded_len,
                        "mismatches": int(cols[2]),
                        "gaps": int(cols[3]),
                        "strand": "Forward" if cols[4] == "plus" else "Reverse-Complement"
                    })
                # Meilleur HSP couvrant le génome
                best = max(hsps, key=lambda x: x["length"])
                pident = best["pident"]
                divergence = round(100.0 - pident, 3)
                aln_len = best["length"]
                mismatches = best["mismatches"]
                gaps = best["gaps"]
                strand = best["strand"]
            else:
                pident, divergence, aln_len, mismatches, gaps, strand = 0.0, 100.0, 0, 0, 0, "Unknown"

        except Exception as e:
            print(f"[-] Erreur BLAST sur {label}: {e}")
            pident, divergence, aln_len, mismatches, gaps, strand = 0.0, 100.0, 0, 0, 0, "Error"

        finally:
            if tmp_target and os.path.exists(tmp_target):
                os.remove(tmp_target)

        results.append({
            "Sample": label,
            "File_Snakemake": sn_val["file"],
            "File_Nextflow": nf_val["file"],
            "Length_Snakemake": sn_val["len"],
            "Length_Nextflow": nf_val["len"],
            "Delta_Length": nf_val["len"] - sn_val["len"],
            "Orientation": strand,
            "Identity_Pct": round(pident, 3),
            "Divergence_Pct": divergence,
            "SNVs_Count": mismatches,
            "Internal_Gaps": gaps,
            "Overlap_Length": aln_len
        })

    return pd.DataFrame(results)

# ==============================================================================
# 6. GÉNÉRATION DE LA VISUALISATION BI-PANNEAU
# ==============================================================================
def generate_plot(df, out_png):
    n_samples = len(df)
    fig_height = max(8, n_samples * 0.55)
    fig, axes = plt.subplots(1, 2, figsize=(22, fig_height), gridspec_kw={'width_ratios': [1.1, 1.4]})

    # Panneau 1 : Longueurs
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

    # Panneau 2 : Heatmap de concordance
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
    axes[1].set_title("Concordance de Séquence", fontsize=13, fontweight="bold")
    axes[1].set_ylabel("")
    axes[1].set_xticks([0.5])
    axes[1].set_xticklabels(["% Identité"], fontweight="bold")

    plt.suptitle(f"Validation VIRiONT_V2 : Snakemake vs Nextflow ({n_samples} séquences)", fontsize=15, fontweight="bold", y=0.99)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

# ==============================================================================
# MAIN ENTRYPOINT
# ==============================================================================
def main():
    args = parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"[*] Dossier Snakemake : {args.snake_dir}")
    print(f"[*] Dossier Nextflow  : {args.nf_dir}")
    print(f"[*] Dossier de sortie : {args.out_dir}")
    print(f"[*] Mode circulaire   : {args.circular}\n")

    if not os.path.isdir(args.snake_dir):
        sys.exit(f"[-] Erreur : Répertoire introuvable : {args.snake_dir}")
    if not os.path.isdir(args.nf_dir):
        sys.exit(f"[-] Erreur : Répertoire introuvable : {args.nf_dir}")

    snake_dict = load_sequences(args.snake_dir, "Snakemake")
    nf_dict = load_sequences(args.nf_dir, "Nextflow")

    pairs = match_pairs(snake_dict, nf_dict)
    if not pairs:
        sys.exit("[-] Aucune paire correspondante trouvée.")

    df = run_alignment(pairs, circular_mode=args.circular)

    tsv_out = os.path.join(args.out_dir, "comparison_metrics.tsv")
    df.to_csv(tsv_out, sep="\t", index=False)
    print(f"\n[+] Tableau TSV écrit : {tsv_out}\n")

    print(df[["Sample", "Length_Snakemake", "Length_Nextflow", "Identity_Pct", "Divergence_Pct", "SNVs_Count", "Internal_Gaps", "Overlap_Length"]].to_string(index=False))

    png_out = os.path.join(args.out_dir, "comparaison_claire_nextflow_vs_snakemake.png")
    generate_plot(df, png_out)
    print(f"\n[✔] Graphique généré avec succès : {png_out}")

if __name__ == "__main__":
    main()