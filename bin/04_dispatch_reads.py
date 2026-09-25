#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import subprocess
from collections import defaultdict
import pysam
from Bio import SeqIO

def main():
    if len(sys.argv) < 5:
        print("Usage: 04_dispatch_reads.py <bam_file> <ref_panel_fasta> <sample_id> <output_dir> [cutoff_ratio=0.05] [min_reads=20]")
        sys.exit(1)

    bam_path = sys.argv[1]
    ref_panel_path = sys.argv[2]
    sample_id = sys.argv[3]
    out_dir = sys.argv[4]
    cutoff_ratio = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
    min_reads_abs = int(sys.argv[6]) if len(sys.argv) > 6 else 20

    os.makedirs(out_dir, exist_ok=True)

    # 1. Indexation locale des références du panel
    ref_sequences = {rec.id: rec for rec in SeqIO.parse(ref_panel_path, "fasta")}

    # 2. Comptage et assignation des reads par référence
    bam = pysam.AlignmentFile(bam_path, "rb")
    reads_per_ref = defaultdict(list)

    for read in bam.fetch(until_eof=True):
        if not read.is_unmapped and read.mapping_quality >= 10:
            reads_per_ref[read.reference_name].append(read.query_name)
    bam.close()

    total_assigned = sum(len(r) for r in reads_per_ref.values())
    print(f"Total reads assignes (MAPQ >= 10) : {total_assigned}")

    # 3. Export du tableau de comptage effectif
    counts_file = f"{sample_id}_real_counts.tsv"
    with open(counts_file, "w") as cf:
        cf.write("cluster\treads\tpercentage\n")
        for ref_name, qnames in sorted(reads_per_ref.items(), key=lambda x: len(x[1]), reverse=True):
            cnt = len(qnames)
            pct = (cnt / total_assigned * 100.0) if total_assigned > 0 else 0.0
            cf.write(f"{ref_name}\t{cnt}\t{pct:.2f}\n")

    # 4. Partitionnement des reads pour chaque variant retenu
    for ref_name, qnames in reads_per_ref.items():
        cnt = len(qnames)
        ratio = cnt / total_assigned if total_assigned > 0 else 0.0
        if ratio >= cutoff_ratio and cnt >= min_reads_abs:
            m_geno = re.search(r"(?:HBV|HDV)[_-]?([A-Za-z0-9]+)", ref_name, re.IGNORECASE)
            geno_tag = m_geno.group(1) if m_geno else re.sub(r"[^A-Za-z0-9]", "", ref_name)

            # Écriture de la référence dédiée au sous-génotype
            sub_ref_path = os.path.join(out_dir, f"{sample_id}_{geno_tag}.ref.fasta")
            if ref_name in ref_sequences:
                SeqIO.write([ref_sequences[ref_name]], sub_ref_path, "fasta")

            # Extraction des reads correspondants
            qnames_txt = os.path.join(out_dir, f"{geno_tag}_names.txt")
            with open(qnames_txt, "w") as qf:
                for name in qnames:
                    qf.write(f"{name}\n")

            fq_out = os.path.join(out_dir, f"{sample_id}_{geno_tag}.fastq.gz")
            cmd = f"samtools view -N {qnames_txt} -b {bam_path} | samtools fastq - | gzip -c > {fq_out}"
            subprocess.run(cmd, shell=True, check=True)

            if os.path.exists(qnames_txt):
                os.remove(qnames_txt)
            print(f"--> Sous-population retenue : {geno_tag} ({cnt} reads, {ratio * 100.0:.1f}%)")

if __name__ == "__main__":
    main()