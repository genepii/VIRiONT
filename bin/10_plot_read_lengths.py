#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    parser = argparse.ArgumentParser(description="Rapport PDF : distribution de taille des reads avec repère biologique dynamique.")
    parser.add_argument("--work-dir", default=".", help="Dossier contenant les FASTQ (*_merged, *_dehosted, *_trimmed)")
    parser.add_argument("--min-len", type=int, default=1000, help="Seuil de coupure basse (pb)")
    parser.add_argument("--virus", default=None, help="Virus (VHB ou VHD)")
    parser.add_argument("--tech", default=None, help="Protocole (WG, POL, R0)")
    parser.add_argument("--out-pdf", default="RUN_READ_LENGTHS_DISTRIBUTION.pdf", help="Fichier PDF de sortie")
    return parser.parse_args()

def extract_lengths(filepath):
    lengths = []
    if not filepath or not os.path.exists(filepath):
        return lengths
    opener = gzip.open if filepath.endswith(".gz") else open
    try:
        with opener(filepath, "rt", errors="ignore") as handle:
            for idx, line in enumerate(handle):
                if idx % 4 == 1:
                    lengths.append(len(line.strip()))
    except Exception as e:
        print(f"Erreur de lecture sur {filepath} : {e}")
    return lengths

def main():
    args = parse_args()

    raw_files = sorted(glob.glob(os.path.join(args.work_dir, "*_merged.fastq.gz")))
    if not raw_files:
        raw_files = sorted(glob.glob(os.path.join(args.work_dir, "*_merged.fastq")))

    if not raw_files:
        print("Aucun fichier FASTQ brut/merged trouvé.")
        return

    path_context = (args.work_dir + " " + " ".join(raw_files)).upper()
    virus = args.virus.upper() if args.virus else ("VHD" if ("VHD" in path_context or "HDV" in path_context) else "VHB")
    
    if args.tech:
        tech = args.tech.upper()
    elif "POL" in path_context:
        tech = "POL"
    elif "R0" in path_context:
        tech = "R0"
    else:
        tech = "WG"

    if virus == "VHD":
        if tech == "R0" or args.min_len < 500:
            target_len = 400
            target_label = "VHD R0 (~400 pb)"
            max_x = 1000
        else:
            target_len = 1700
            target_label = "VHD WG (~1,7 kb)"
            max_x = 2500
    else:
        if tech == "POL":
            target_len = 1200
            target_label = "VHB POL (~1,2 kb)"
            max_x = 2500
        else:
            target_len = 3200
            target_label = "VHB (~3,2 kb)"
            max_x = 4000

    bins = np.linspace(0, max_x, 80)
    stages = [
        ("01_RAW / MERGED", "_merged", "#1f77b4"),
        ("02_DEHOSTING", "_dehosted", "#ff7f0e"),
        ("03_FILTERED / TRIMMED", "_trimmed", "#2ca02c")
    ]

    print(f"Génération PDF [{virus} - {tech}] (Repère : {target_label}) -> {args.out_pdf}")

    retention_rows = []

    with PdfPages(args.out_pdf) as pdf:
        for raw_path in raw_files:
            filename = os.path.basename(raw_path)
            sample_id = re.sub(r"_merged\.fastq(\.gz)?$", "", filename)

            dehost_path = os.path.join(args.work_dir, f"{sample_id}_dehosted.fastq.gz")
            if not os.path.exists(dehost_path):
                dehost_path = os.path.join(args.work_dir, f"{sample_id}_dehosted.fastq")

            trimmed_path = os.path.join(args.work_dir, f"{sample_id}_trimmed.fastq.gz")
            if not os.path.exists(trimmed_path):
                trimmed_path = os.path.join(args.work_dir, f"{sample_id}_trimmed.fastq")

            paths = [raw_path, dehost_path, trimmed_path]
            data = [extract_lengths(p) for p in paths]

            n_raw = len(data[0])
            n_dehost = len(data[1])
            n_trim = len(data[2])

            ret_host = round((n_dehost / n_raw * 100), 2) if n_raw > 0 else 0.0
            ret_trim = round((n_trim / n_raw * 100), 2) if n_raw > 0 else 0.0

            retention_rows.append({
                "Sample": sample_id,
                "Rétention hôte %": ret_host,
                "Reads conservés post-filtre %": ret_trim
            })

            fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=False)

            for ax, (title, _, color), lens in zip(axes, stages, data):
                if lens:
                    ax.hist(lens, bins=bins, color=color, alpha=0.75, edgecolor="black", linewidth=0.5)

                ax.axvline(args.min_len, color="#d62728", linestyle="--", linewidth=1.2, label=f"Seuil ({args.min_len} pb)")
                ax.axvline(target_len, color="#6a3d9a", linestyle=":", linewidth=1.5, label=target_label)

                ax.set_title(f"{title}\n({len(lens):,} reads)", fontsize=11, fontweight="bold")
                ax.set_xlabel("Longueur (pb)", fontsize=9.5)
                ax.set_ylabel("Nombre de reads", fontsize=9.5)
                ax.set_xlim(0, max_x)
                ax.grid(axis="y", linestyle="--", alpha=0.3)
                ax.legend(loc="upper right", fontsize=8)

            titre_principal = f"Échantillon : {sample_id}\nRétention hôte : {ret_host:.1f} %  |  Reads conservés post-filtre : {ret_trim:.1f} %"
            fig.suptitle(titre_principal, fontsize=12, fontweight="bold", y=1.05)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    # Sauvegarde automatique du TSV de rétention pour injection dans la table finale
    tsv_out = os.path.join(args.work_dir, "read_retention_metrics.tsv")
    pd.DataFrame(retention_rows).to_csv(tsv_out, sep="\t", index=False)
    print(f"Table intermédiaire des métriques de rétention enregistrée : {tsv_out}")
    print("Rapport PDF généré avec succès.")

if __name__ == "__main__":
    main()