#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import csv
from collections import defaultdict


def parse_fasta(fasta_file):
    sequences = {}
    current_sample = None
    seq_lines = []

    if not fasta_file or not os.path.exists(fasta_file):
        return sequences

    with open(fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_sample and seq_lines:
                    if current_sample not in sequences:
                        sequences[current_sample] = "".join(seq_lines).replace("\r", "").replace("\n", "")
                header = line[1:].split()[0]
                if "_Geno_" in header:
                    current_sample = header.split('_Geno_')[0]
                else:
                    current_sample = header
                seq_lines = []
            else:
                seq_lines.append(line)

        if current_sample and seq_lines:
            if current_sample not in sequences:
                sequences[current_sample] = "".join(seq_lines).replace("\r", "").replace("\n", "")

    return sequences


def clean_sample_id(raw_sample):
    parts = raw_sample.split('_')
    if len(parts) >= 3 and parts[0] == 'barcode':
        return "_".join(parts[2:])
    return raw_sample


def safe_reads_int(val):
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return 0


def load_alerts(alert_file, alert_prefix="AvisBio_Mutation"):
    """
    Lit un fichier d'alertes GLIMS (variants ou indels)
    et renvoie un dictionnaire associant l'ID épuré et le nom brut
    au libellé clinique formaté :
    { "25063706601": "AvisBio_Indel (Domaine HBx, AUTRE_REGION)" }
    """
    alerts_by_sample = defaultdict(set)
    if not alert_file or not os.path.exists(alert_file) or os.path.getsize(alert_file) == 0:
        return {}

    with open(alert_file, 'r', encoding='utf-8') as f:
        sample_chunk = f.read(2048)
        f.seek(0)
        delim = ';' if ';' in sample_chunk else ','
        reader = csv.DictReader(f, delimiter=delim)

        for row in reader:
            raw_sample = row.get("SAMPLE", "").strip() or row.get("Sample", "").strip() or row.get("SAMPLE_ID", "").strip()
            gene = row.get("GENE", "").strip()
            if not raw_sample:
                continue

            cleaned = clean_sample_id(raw_sample)
            if gene:
                alerts_by_sample[cleaned].add(gene)
                alerts_by_sample[raw_sample].add(gene)
            else:
                alerts_by_sample[cleaned].add("Inconnu")
                alerts_by_sample[raw_sample].add("Inconnu")

    formatted = {}
    for sid, genes in alerts_by_sample.items():
        genes_list = [g for g in sorted(genes) if g != "Inconnu"]
        if genes_list:
            formatted[sid] = f"{alert_prefix} ({', '.join(genes_list)})"
        else:
            formatted[sid] = alert_prefix
    return formatted


def main():
    parser = argparse.ArgumentParser(description="Génère le CSV FastFinder/GLIMS pour VIRiONT V2/V3")
    parser.add_argument("--summary", required=True, help="Chemin vers SUMMARY_Multi_Infection.tsv")
    parser.add_argument("--fasta", required=True, help="Chemin vers validated_consensus_all.fasta")
    parser.add_argument("--mutations", default=None, help="Chemin vers ALL_SAMPLES_glims_variants_alerts.csv")
    parser.add_argument("--indels", default=None, help="Chemin vers ALL_SAMPLES_glims_indels_alerts.csv")
    parser.add_argument("--run-id", required=True, help="Identifiant du run (ex: 260931_VHB_WG)")
    parser.add_argument("--output", required=True, help="Fichier CSV de sortie")
    args = parser.parse_args()

    fasta_seqs = parse_fasta(args.fasta)
    mut_alerts = load_alerts(args.mutations, alert_prefix="AvisBio_Mutation")
    indel_alerts = load_alerts(args.indels, alert_prefix="AvisBio_Indel")

    samples_data = defaultdict(list)
    if os.path.exists(args.summary):
        with open(args.summary, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample_col = row.get('sample') or row.get('Sample')
                if sample_col:
                    samples_data[sample_col].append(row)

    headers = [
        "Sample ID",
        "AssayResultTargetCode",
        "Instrument_Id",
        "Target_1_result",
        "Target_2_result",
        "Target_3_result",
        "Target_4_result",
        "Target_5_result"
    ]

    run_upper = args.run_id.upper()

    with open(args.output, 'w', newline='', encoding='utf-8') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(headers)

        for raw_sample, rows in samples_data.items():
            sample_id = clean_sample_id(raw_sample)
            protocol = rows[0].get('protocol', 'WG').upper()

            if "VHD" in run_upper or "HDV" in run_upper or "VHD" in protocol or "HDV" in protocol:
                if "R0" in run_upper or "R0" in protocol:
                    assay_code = "GENOVHDR0"
                else:
                    assay_code = "GENOVHDWG"
            elif "POL" in run_upper or "POL" in protocol:
                assay_code = "GENOVHBPOL"
            else:
                assay_code = "GENOVHBWG"

            total_reads = sum(safe_reads_int(r.get('total_reads', 0)) for r in rows)

            coinf_candidates = [
                r for r in rows
                if total_reads > 0 and (safe_reads_int(r.get('total_reads', 0)) / total_reads * 100.0) > 30.0
            ]

            validated_geno = "ND"
            for r in rows:
                if r.get('status', '').upper() == 'VALIDATED':
                    validated_geno = r.get('genotype', 'ND')
                    break

            if total_reads < 50 or validated_geno == "ND":
                target_2 = "ININT"
            elif 50 <= total_reads <= 500:
                target_2 = "AVIS_BIO_DEPTH"
            elif len(coinf_candidates) >= 2:
                target_2 = "AVIS_BIO_COINF"
            else:
                target_2 = validated_geno

            # Target_3 : Substitutions / Variants
            target_3 = mut_alerts.get(sample_id, mut_alerts.get(raw_sample, "RAS"))
           
            # Target_4 : INDELs
            target_4 = indel_alerts.get(sample_id, indel_alerts.get(raw_sample, "RAS"))
           
            # Target_5 : Séquence consensus
            target_5 = fasta_seqs.get(raw_sample, "")
            if not target_5:
                for seq_key, seq_val in fasta_seqs.items():
                    if sample_id in seq_key or raw_sample in seq_key:
                        target_5 = seq_val
                        break

            writer.writerow([
                sample_id,
                assay_code,
                "ONT",
                args.run_id,
                target_2,
                target_3,
                target_4,
                target_5
            ])

    print(f"[FASTFINDER] Export GLIMS (Target_1 à 5) généré avec succès : {args.output}")


if __name__ == "__main__":
    main()