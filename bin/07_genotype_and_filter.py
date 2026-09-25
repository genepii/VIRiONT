#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import glob
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def load_references(ref_fasta_path):
    refs = {}
    if ref_fasta_path and os.path.exists(ref_fasta_path):
        for rec in SeqIO.parse(ref_fasta_path, "fasta"):
            clean_name = re.sub(r"[^A-Za-z0-9]", "", rec.id).upper()
            refs[rec.id] = rec
            refs[clean_name] = rec
    return refs


def find_ref_record(genotype, refs_dict):
    clean_g = re.sub(r"[^A-Za-z0-9]", "", genotype).upper()
    if clean_g in refs_dict:
        return refs_dict[clean_g]
    for k, v in refs_dict.items():
        if clean_g in k or k in clean_g:
            return v
    return None


def ensure_plus_strand(consensus_seq, ref_seq):
    """
    Vérifie l'orientation par k-mer matching rapide contre la référence.
    Si le brin inverse donne un score supérieur, on reverse-complémente.
    """
    seq_str = str(consensus_seq).upper()
    ref_str = str(ref_seq).upper()
    k = 21

    def count_kmers(q_seq, target):
        score = 0
        step = max(1, len(q_seq) // 200)
        for i in range(0, len(q_seq) - k, step):
            if q_seq[i:i + k] in target:
                score += 1
        return score

    plus_score = count_kmers(seq_str, ref_str)
    rev_seq = str(Seq(seq_str).reverse_complement())
    minus_score = count_kmers(rev_seq, ref_str)

    if minus_score > plus_score:
        return rev_seq, True
    return seq_str, False


def trim_concatemer(seq_str, expected_len):
    """
    Raccourcit le consensus s'il dépasse 115% de la taille attendue (rolling-circle/chimeric).
    """
    if not expected_len or len(seq_str) <= int(expected_len * 1.15):
        return seq_str
    print(f"✂️ Curation concatémère : {len(seq_str)} bp -> ramené à {expected_len} bp")
    return seq_str[:expected_len]


def main():
    if len(sys.argv) < 7:
        print("Usage: 07_genotype_and_filter.py <consensus_dir> <fastqs_dir> <counts_dir> <ref_fasta> <cutoff_percent> [virus_name]")
        sys.exit(1)

    consensus_dir = sys.argv[1]
    fastqs_dir    = sys.argv[2]
    counts_dir    = sys.argv[3]
    ref_path      = sys.argv[4]
    cutoff        = float(sys.argv[5])
    virus_name    = sys.argv[6] if len(sys.argv) > 6 else "VHB"

    active_refs = load_references(ref_path)

    summary_rows = []
    validated_records = []
    manifest_rows = []
    os.makedirs("renamed_consensus", exist_ok=True)

    # Récupérer tous les fichiers de comptage réel issus de PRONAME (04_competitive_align)
    count_files = sorted(glob.glob(os.path.join(counts_dir, "*_real_counts.tsv")))

    for c_file in count_files:
        sample_id = os.path.basename(c_file).replace("_real_counts.tsv", "")
        protocol = "WG"

        try:
            df_counts = pd.read_csv(c_file, sep="\t")
        except Exception as e:
            print(f"⚠️ Erreur lecture {c_file}: {e}")
            continue

        if df_counts.empty:
            continue

        # Calcul du total des reads assignés pour ce barcode
        total_sample_reads = df_counts["reads"].sum()
        max_reads = df_counts["reads"].max() or 1

        for _, row in df_counts.iterrows():
            ref_name = str(row["cluster"]).strip()
            reads = int(row["reads"])

            # Nettoyage du génotype (ex: HBV_A2 -> A2)
            m_geno = re.search(r"(?:HBV|HDV)[_-]?([A-Za-z0-9]+)", ref_name, re.IGNORECASE)
            genotype = m_geno.group(1) if m_geno else re.sub(r"[^A-Za-z0-9]", "", ref_name)

            # Ratio par rapport au variant dominant (%)
            ratio_relative = (reads / max_reads) * 100.0
            # Pourcentage global sur le barcode (%)
            fraction_global = (reads / total_sample_reads) * 100.0 if total_sample_reads > 0 else 0.0

            # Bruit négligeable (< 5% des reads ou < 20 reads)
            if fraction_global < 5.0 and reads < 20:
                continue

            # Règle clinique
            status = "VALIDATED" if ratio_relative >= cutoff else "REJECTED"

            # Recherche du fichier consensus poli par Medaka
            possible_fasta = [
                os.path.join(consensus_dir, sample_id, f"{sample_id}_Geno_{genotype}.fasta"),
                os.path.join(consensus_dir, f"{sample_id}_Geno_{genotype}.fasta"),
                os.path.join(consensus_dir, f"{sample_id}.fasta")
            ]

            cons_file = next((p for p in possible_fasta if os.path.exists(p) and os.path.getsize(p) > 0), None)

            seq_str = ""
            final_len = 0
            if cons_file:
                recs = list(SeqIO.parse(cons_file, "fasta"))
                if recs:
                    seq_str = str(recs[0].seq)

            ref_record = find_ref_record(genotype, active_refs)
            expected_len = len(ref_record.seq) if ref_record else None

            if seq_str:
                if ref_record:
                    seq_str, flipped = ensure_plus_strand(seq_str, str(ref_record.seq))
                    if flipped:
                        print(f"🔄 {sample_id}_{genotype} : réorienté sur le brin (+)")
                seq_str = trim_concatemer(seq_str, expected_len)
                final_len = len(seq_str)

            summary_rows.append({
                "sample": sample_id,
                "genotype": genotype,
                "protocol": protocol,
                "best_cluster": f"{sample_id}_{genotype}",
                "total_reads": reads,
                "ratio_percent": f"{ratio_relative:.2f}",
                "ratio_num": ratio_relative,
                "pident": 100.0,
                "length": final_len,
                "strand": "plus",
                "status": status,
                "merge_note": f"Alignement compétitif ({fraction_global:.1f}% des reads totaux)"
            })

            # Exportation du consensus individuel
            if seq_str:
                final_id = f"{sample_id}_Geno_{genotype}"
                out_rec = SeqIO.SeqRecord(
                    Seq(seq_str),
                    id=final_id,
                    description=f"reads={reads} ratio={ratio_relative:.2f}% status={status}"
                )

                SeqIO.write([out_rec], os.path.join("renamed_consensus", f"{final_id}.fasta"), "fasta")

                if status == "VALIDATED":
                    validated_records.append(out_rec)
                    manifest_rows.append({
                        "sample_id": sample_id,
                        "genotype": genotype,
                        "protocol": protocol,
                        "filename": f"{final_id}.fasta"
                    })

    # Écriture des rapports cliniques finaux
    pd.DataFrame(summary_rows).to_csv("SUMMARY_Multi_Infection.tsv", sep="\t", index=False)
    SeqIO.write(validated_records, "validated_consensus_all.fasta", "fasta")

    pd.DataFrame(
        manifest_rows,
        columns=["sample_id", "genotype", "protocol", "filename"]
    ).to_csv(os.path.join("renamed_consensus", "manifest.tsv"), sep="\t", index=False)

    print(f"Terminé : {len(summary_rows)} variants documentés ({len(validated_records)} validés).")


if __name__ == "__main__":
    main()