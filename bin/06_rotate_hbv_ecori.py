#!/usr/bin/env python3
import sys
import subprocess

def rotate_to_canonical(input_fasta, ref_fasta, output_fasta):
    # 1. Lecture du consensus requête
    header = ""
    seq_chunks = []
    with open(input_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
            else:
                seq_chunks.append(line.strip())
    query_seq = "".join(seq_chunks).upper()

    if not query_seq:
        sys.exit(f"Erreur : séquence vide dans {input_fasta}")

    # 2. Extraire les 80 premières bases de la ref EcoRI jpHMM
    ref_chunks = []
    with open(ref_fasta, "r") as f:
        for line in f:
            if not line.startswith(">"):
                ref_chunks.append(line.strip())
                if sum(len(c) for c in ref_chunks) >= 80:
                    break
    seed = "".join(ref_chunks)[:80].upper()

    # 3. Aligner le fragment initial sur le consensus
    # On double virtuellement la séquence requête pour gérer la circularité
    doubled_seq = query_seq + query_seq
    
    cmd = [
        "minimap2", "-x", "sr", "-k", "12", "-w", "4", "--secondary=no",
        input_fasta, "-"
    ]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    paf_out, _ = p.communicate(input=f">seed\n{seed}\n")

    cut_pos = 0
    for line in paf_out.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 9:
            target_strand = parts[4]
            target_start = int(parts[7])
            query_start = int(parts[2])
            
            if target_strand == "+":
                cut_pos = target_start - query_start
            else:
                cut_pos = target_start + query_start
            break

    # 4. Rotation circulaire
    if 0 < cut_pos < len(query_seq):
        rotated_seq = query_seq[cut_pos:] + query_seq[:cut_pos]
    else:
        rotated_seq = query_seq

    with open(output_fasta, "w") as out:
        out.write(f"{header}\n")
        for i in range(0, len(rotated_seq), 60):
            out.write(rotated_seq[i:i+60] + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit("Usage: 06_rotate_hbv_ecori.py <in.fasta> <ref.fasta> <out.fasta>")
    rotate_to_canonical(sys.argv[1], sys.argv[2], sys.argv[3])