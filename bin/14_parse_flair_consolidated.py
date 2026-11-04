#!/usr/bin/env python3
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Consolidation et filtrage du splicing FLAIR.")
    parser.add_argument("sample_id", help="Identifiant de l'échantillon")
    parser.add_argument("counts_file", help="Fichier counts.txt de FLAIR")
    parser.add_argument("bed_file", help="Fichier isoforms.bed de FLAIR")
    parser.add_argument("catalog_file", help="Catalogue d'épissage TSV")
    parser.add_argument("output_file", help="Fichier CSV de sortie")
    parser.add_argument("--min-pct", type=float, default=10.0, help="Pourcentage minimum de reads (défaut: 10.0)")
    parser.add_argument("--min-size", type=int, default=150, help="Taille minimale de l'intron en pb (défaut: 150)")
    parser.add_argument("--min-reads", type=int, default=1000, help="Nombre minimum absolu de reads (défaut: 1000)")
    parser.add_argument("--tolerance", type=int, default=15, help="Tolérance de position en pb (défaut: 15)")
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Charger le catalogue d'épissage VHB
    catalog = []
    if os.path.exists(args.catalog_file):
        with open(args.catalog_file, "r") as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 7:
                    name = parts[0]
                    p_start = int(parts[3])  # PIPE_START
                    p_end = int(parts[4])    # PIPE_END
                    tr_type = parts[6]       # TYPE
                    catalog.append({'name': name, 'start': p_start, 'end': p_end, 'type': tr_type})

    # 2. Charger les comptes globaux
    counts = {}
    total_reads = 0
    if os.path.exists(args.counts_file):
        with open(args.counts_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    counts[parts[0]] = int(parts[1])
                    total_reads += int(parts[1])

    # 3. Extraire les introns du BED FLAIR
    raw_introns = []
    if os.path.exists(args.bed_file):
        with open(args.bed_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 12:
                    isoform_id = parts[3]
                    block_count = int(parts[9])
                    block_sizes = [int(x) for x in parts[10].rstrip(',').split(',')]
                    block_starts = [int(x) for x in parts[11].rstrip(',').split(',')]
                    chrom_start = int(parts[1])
                    
                    if block_count > 1:
                        for i in range(block_count - 1):
                            e_end = chrom_start + block_starts[i] + block_sizes[i]
                            e_start = chrom_start + block_starts[i+1]
                            reads = counts.get(isoform_id, 0)
                            raw_introns.append((e_end, e_start, reads))

    # 4. Clustering initial des introns avec tolérance de position
    clusters = []
    for start, end, reads in sorted(raw_introns, key=lambda x: x[2], reverse=True):
        merged = False
        for cluster in clusters:
            if abs(cluster['start'] - start) <= args.tolerance and abs(cluster['end'] - end) <= args.tolerance:
                cluster['reads'] += reads
                merged = True
                break
        if not merged:
            clusters.append({'start': start, 'end': end, 'reads': reads})

    clusters.sort(key=lambda x: x['reads'], reverse=True)

    # 5. Filtrage (Reads min, Taille min, Pourcentage min) et annotation
    with open(args.output_file, "w") as out:
        out.write("Sample,Intron_Start,Intron_End,Size_bp,Reads,Percentage,Transcript,Type\n")
        for c in clusters:
            size = c['end'] - c['start']
            pct = (c['reads'] / total_reads) * 100 if total_reads > 0 else 0.0

            # FILTRES ANTI-BRUIT DE FOND
            if c['reads'] < args.min_reads:
                continue
            if size < args.min_size:
                continue
            if pct < args.min_pct:
                continue

            matched_name = "UNKNOWN"
            matched_type = "NOVEL"
            
            for cat in catalog:
                if (abs(c['start'] - cat['start']) <= args.tolerance) and (abs(c['end'] - cat['end']) <= args.tolerance):
                    matched_name = cat['name']
                    matched_type = cat['type']
                    break
            
            out.write(f"{args.sample_id},{c['start']},{c['end']},{size},{c['reads']},{pct:.2f},{matched_name},{matched_type}\n")

    print(f"[FLAIR] {args.sample_id} filtré : min_size={args.min_size}bp, min_pct={args.min_pct}%")

if __name__ == "__main__":
    main()