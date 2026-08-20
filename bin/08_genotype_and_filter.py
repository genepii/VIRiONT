#!/usr/bin/env python3
import sys
import os
import glob
import re
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(seq_str):
    return str(Seq(seq_str).reverse_complement())


def load_ref_lengths(ref_fasta_path):
    ref_lengths = {}
    if ref_fasta_path and os.path.exists(ref_fasta_path):
        for rec in SeqIO.parse(ref_fasta_path, "fasta"):
            ref_lengths[rec.id] = len(rec.seq)
    return ref_lengths


def main():
    if len(sys.argv) < 5:
        print("Usage: python3 08_genotype_and_filter.py <consensus_dir> <csvs_dir> <blast_db> <cutoff_percent> [ref_fasta_for_lengths]")
        sys.exit(1)

    consensus_dir = sys.argv[1]
    csvs_dir      = sys.argv[2]
    blast_db      = sys.argv[3]
    cutoff        = float(sys.argv[4])
    ref_fasta_for_lengths = sys.argv[5] if len(sys.argv) > 5 else None

    ref_lengths = load_ref_lengths(ref_fasta_for_lengths)

    summary_rows = []
    validated_records = []
    manifest_rows = []

    fasta_files = sorted(glob.glob(os.path.join(consensus_dir, "*.fasta")))

    for f_path in fasta_files:
        sample_id = os.path.basename(f_path).replace(".fasta", "")
        csv_file = os.path.join(csvs_dir, f"{sample_id}_results.csv")

        # 1. Lecture des effectifs réels
        reads_dict = {}
        if os.path.exists(csv_file):
            try:
                df_csv = pd.read_csv(csv_file, sep="\t", header=0)
                for _, row in df_csv.iterrows():
                    c_name = str(row["cluster_id"]).strip()
                    try:
                        c_reads = int(row["real_read_count"])
                        if c_name and c_reads > 0:
                            reads_dict[c_name] = c_reads
                    except Exception:
                        pass
            except Exception:
                try:
                    df_csv = pd.read_csv(csv_file, header=None)
                    for _, row in df_csv.iterrows():
                        c_name = str(row[0]).strip()
                        c_reads = int(row[1])
                        if c_name and c_name != "Total" and c_reads > 0:
                            reads_dict[c_name] = c_reads
                except Exception:
                    pass

        # 2. BLASTN
        blast_out = f"{sample_id}_blast.tsv"
        cmd = [
            "blastn", "-query", f_path, "-db", blast_db,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-max_target_seqs", "1"
        ]
        with open(blast_out, "w") as out_f:
            subprocess.run(cmd, stdout=out_f, check=True)

        blast_hits = {}
        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
            with open(blast_out) as bf:
                for line in bf:
                    parts = line.strip().split("\t")
                    if parts[0] not in blast_hits:
                        blast_hits[parts[0]] = {
                            "genotype": parts[1],
                            "pident": float(parts[2]),
                            "length": int(parts[3]),
                            "sstart": int(parts[8]),
                            "send": int(parts[9])
                        }

        # 3. Attribution des clusters bruts (On ne retient QUE les vrais hits BLAST)
        records = list(SeqIO.parse(f_path, "fasta"))
        raw_clusters = []
        for idx, rec in enumerate(records):
            # Si aucun hit BLAST n'existe pour cette séquence, on l'exclut directement
            if rec.id not in blast_hits:
                continue

            hit = blast_hits[rec.id]
            geno = hit["genotype"]

            # Exclusion des Unknown ou No_Match
            if geno.upper() in ["UNKNOWN", "NO_MATCH", "NONE", ""]:
                continue

            c_name = rec.id
            reads = 0
            if c_name in reads_dict:
                reads = reads_dict[c_name]
            else:
                for k, v in reads_dict.items():
                    if k in c_name or c_name in k:
                        reads = v
                        break
            if reads == 0:
                m = re.search(r'\((\d+)\)', c_name)
                if m:
                    reads = int(m.group(1))
            if reads == 0 and len(reads_dict) > idx:
                reads = list(reads_dict.values())[idx]
            if reads == 0:
                reads = 1

            strand = "minus" if hit["sstart"] > hit["send"] else "plus"

            raw_clusters.append({
                "record": rec,
                "reads": reads,
                "genotype": geno,
                "pident": hit["pident"],
                "length_seq": len(rec.seq),
                "strand": strand,
                "cluster_label": f"c{idx + 1}"
            })

        if not raw_clusters:
            continue

        # 4. AGRÉGATION PAR GÉNOTYPE
        genotype_groups = {}
        for c in raw_clusters:
            geno = c["genotype"]
            genotype_groups.setdefault(geno, []).append(c)

        aggregated_genotypes = []
        for geno, members in genotype_groups.items():
            total_reads = sum(m["reads"] for m in members)
            target_len = ref_lengths.get(geno, None)
            if target_len is not None:
                best_member = min(members, key=lambda x: abs(x["length_seq"] - target_len))
            else:
                best_member = max(members, key=lambda x: x["length_seq"])
            weighted_pident = sum(m["pident"] * m["reads"] for m in members) / total_reads
            merged_clusters = ",".join(m["cluster_label"] for m in members)
            aggregated_genotypes.append({
                "genotype": geno,
                "record": best_member["record"],
                "reads": total_reads,
                "pident": round(weighted_pident, 3),
                "length_seq": best_member["length_seq"],
                "strand": best_member["strand"],
                "cluster_label": "c1" if len(members) == 1 else f"merged({merged_clusters})",
                "is_merged": len(members) > 1,
                "merged_from": merged_clusters
            })

        # 5. Calcul des ratios cliniques
        max_reads = max([g["reads"] for g in aggregated_genotypes]) or 1
        for g in aggregated_genotypes:
            ratio = (g["reads"] / max_reads) * 100.0
            status = "VALIDATED" if (g["pident"] >= 70.0 and ratio >= cutoff) else "REJECTED"

            summary_rows.append({
                "sample": sample_id,
                "genotype": g["genotype"],
                "best_cluster": f"{sample_id}_{g['genotype']}",
                "total_reads": g["reads"],
                "ratio_percent": f"{ratio:.2f}",
                "ratio_num": ratio,
                "pident": g["pident"],
                "length": g["length_seq"],
                "strand": g["strand"],
                "status": status,
                "merge_note": f"Regroupement de {g['merged_from']}" if g["is_merged"] else ""
            })

            if status == "VALIDATED":
                final_seq = reverse_complement(str(g["record"].seq)) if g["strand"] == "minus" else str(g["record"].seq)
                final_id = f"{sample_id}_Geno_{g['genotype']}"
                new_rec = g["record"]
                new_rec.id = final_id
                new_rec.description = ""
                new_rec.seq = Seq(final_seq)
                validated_records.append(new_rec)
                manifest_rows.append({
                    "sample_id": sample_id,
                    "genotype": g["genotype"],
                    "filename": f"{final_id}.fasta"
                })

    # 6. Écriture des fichiers de sortie
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv("SUMMARY_Multi_Infection.tsv", sep="\t", index=False)

    SeqIO.write(validated_records, "validated_consensus_all.fasta", "fasta")

    os.makedirs("renamed_consensus", exist_ok=True)
    for rec in validated_records:
        SeqIO.write([rec], os.path.join("renamed_consensus", f"{rec.id}.fasta"), "fasta")

    df_manifest = pd.DataFrame(manifest_rows, columns=["sample_id", "genotype", "filename"])
    df_manifest.to_csv(os.path.join("renamed_consensus", "manifest.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()