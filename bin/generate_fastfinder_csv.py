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

    if not os.path.exists(fasta_file):
        return sequences

    with open(fasta_file, 'r') as f:
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
    # Transforme barcode_01_25063706601 -> 25063706601
    parts = raw_sample.split('_')
    if len(parts) >= 3 and parts[0] == 'barcode':
        return "_".join(parts[2:])
    return raw_sample


def safe_reads_int(val):
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return 0


def main():
    parser = argparse.ArgumentParser(description="Génère le CSV FastFinder/GLIMS pour VIRiONT V2")
    parser.add_argument("--summary", required=True, help="Chemin vers SUMMARY_Multi_Infection.tsv")
    parser.add_argument("--fasta", required=True, help="Chemin vers validated_consensus_all.fasta")
    parser.add_argument("--run-id", required=True, help="Identifiant du run (ex: 260931_VHD_R0)")
    parser.add_argument("--output", required=True, help="Fichier CSV de sortie")
    args = parser.parse_args()

    # 1. Lecture des consensus validés
    fasta_seqs = parse_fasta(args.fasta)

    # 2. Lecture du tableau de génotypage
    samples_data = defaultdict(list)
    if os.path.exists(args.summary):
        with open(args.summary, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('sample'):
                    samples_data[row['sample']].append(row)

    # 3. Export GLIMS / FastFinder
    headers = [
        "Sample ID",
        "AssayResultTargetCode",
        "Instrument_Id",
        "Target_1_result",
        "Target_2_result",
        "Target_3_result"
    ]

    run_upper = args.run_id.upper()

    with open(args.output, 'w', newline='', encoding='utf-8') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(headers)

        for raw_sample, rows in samples_data.items():
            sample_id = clean_sample_id(raw_sample)
            protocol = rows[0].get('protocol', 'WG').upper()

            # Attribution dynamique du code d'analyse GLIMS
            if "VHD" in run_upper or "HDV" in run_upper or "VHD" in protocol or "HDV" in protocol:
                if "R0" in run_upper or "R0" in protocol:
                    assay_code = "GENOVHDR0"
                else:
                    assay_code = "GENOVHDWG"
            elif "POL" in run_upper or "POL" in protocol:
                assay_code = "GENOVHBPOL"
            else:
                assay_code = "GENOVHBWG"

            # Somme des reads attribués
            total_reads = sum(safe_reads_int(r.get('total_reads', 0)) for r in rows)

            # Candidats co-infection (> 30 % des reads totaux)
            coinf_candidates = [
                r for r in rows
                if total_reads > 0 and (safe_reads_int(r.get('total_reads', 0)) / total_reads * 100.0) > 30.0
            ]

            # Recherche exclusive du génotype validé cliniquement
            validated_geno = "ND"
            for r in rows:
                if r.get('status', '').upper() == 'VALIDATED':
                    validated_geno = r.get('genotype', 'ND')
                    break

            # Arbre de décision biologique GLIMS
            if total_reads < 50 or validated_geno == "ND":
                target_2 = "ININT"
            elif 50 <= total_reads <= 500:
                target_2 = "AVIS_BIO_DEPTH"
            elif len(coinf_candidates) >= 2:
                target_2 = "AVIS_BIO_COINF"
            else:
                target_2 = validated_geno

            # Consensus associé (vide si non validé ou rejeté)
            target_3 = fasta_seqs.get(raw_sample, "")

            writer.writerow([
                sample_id,
                assay_code,
                "ONT",
                args.run_id,
                target_2,
                target_3
            ])

    print(f"Fichier FastFinder généré avec succès : {args.output}")


if __name__ == "__main__":
    main()