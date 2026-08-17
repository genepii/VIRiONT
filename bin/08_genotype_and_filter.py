#!/usr/bin/env python3
import sys
import os
import glob
import re
import subprocess
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(seq_str):
    return str(Seq(seq_str).reverse_complement())


def blast_pairwise_identity_coverage(seq1, seq2):
    """Compare deux séquences via blastn local alignment, en tenant compte de
    TOUS les HSPs (pas seulement le meilleur), car un gros gap interne
    (consensus partiel/tronqué) génère plusieurs blocs d'alignement séparés
    plutôt qu'un seul alignement continu. On agrège la couverture totale
    (non chevauchante) sur la query et l'identité moyenne pondérée par
    la longueur de chaque bloc."""
    if not seq1 or not seq2:
        return 0.0, 0.0

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qf, \
         tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as sf:
        qf.write(f">q\n{seq1}\n")
        sf.write(f">s\n{seq2}\n")
        q_path, s_path = qf.name, sf.name

    try:
        cmd = [
            "blastn", "-query", q_path, "-subject", s_path,
            "-outfmt", "6 pident length qlen slen qstart qend",
            "-evalue", "1e-10"
            # Pas de -max_hsps : on veut TOUS les blocs d'alignement
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        lines = [l for l in result.stdout.strip().split("\n") if l.strip()]
        if not lines:
            return 0.0, 0.0

        qlen, slen = None, None
        intervals = []       # (qstart, qend) de chaque HSP, pour couverture non-chevauchante
        weighted_identity_sum = 0.0
        total_aln_len = 0

        for line in lines:
            parts = line.split("\t")
            pident = float(parts[0])
            aln_len = int(parts[1])
            qlen, slen = int(parts[2]), int(parts[3])
            qstart, qend = int(parts[4]), int(parts[5])
            if qstart > qend:
                qstart, qend = qend, qstart
            intervals.append((qstart, qend))
            weighted_identity_sum += pident * aln_len
            total_aln_len += aln_len

        if total_aln_len == 0 or qlen is None:
            return 0.0, 0.0

        # Fusion des intervalles qui se chevauchent, pour ne pas compter deux fois
        # la même portion de query si plusieurs HSPs se recoupent partiellement
        intervals.sort()
        merged = [intervals[0]]
        for start, end in intervals[1:]:
            last_start, last_end = merged[-1]
            if start <= last_end + 1:
                merged[-1] = (last_start, max(last_end, end))
            else:
                merged.append((start, end))
        covered_len = sum(end - start + 1 for start, end in merged)

        min_len = min(qlen, slen)
        coverage = (covered_len / min_len) * 100.0 if min_len > 0 else 0.0
        avg_identity = weighted_identity_sum / total_aln_len

        return avg_identity, coverage
    except Exception as e:
        print(f"⚠️ Erreur blastn pairwise: {e}")
        return 0.0, 0.0
    finally:
        os.remove(q_path)
        os.remove(s_path)


def load_ref_lengths(ref_fasta_path):
    ref_lengths = {}
    for rec in SeqIO.parse(ref_fasta_path, "fasta"):
        ref_lengths[rec.id] = len(rec.seq)
    return ref_lengths


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[ry] = rx


def merge_split_clusters(cluster_info, ref_lengths, identity_threshold=97.0,
                          min_coverage=85.0, min_overlap=500):
    """Fusionne les clusters VALIDATED du même génotype qui représentent en réalité
    le même variant, scindé artificiellement par amplicon_sorter (amplicons
    chevauchants, ou split par bruit de séquençage). Utilise blastn local pour
    trouver la vraie région de recouvrement (au lieu d'un alignement global forcé),
    et exige à la fois une identité élevée ET une couverture suffisante de la
    portion commune, pour ne pas fusionner deux amplicons qui se chevauchent
    seulement partiellement mais représentent des positions génomiques différentes."""
    validated_idx = [i for i, c in enumerate(cluster_info) if c["status_tmp"] == "VALIDATED"]
    uf = UnionFind(len(cluster_info))

    for a in range(len(validated_idx)):
        for b in range(a + 1, len(validated_idx)):
            i, j = validated_idx[a], validated_idx[b]
            c1, c2 = cluster_info[i], cluster_info[j]
            if c1["genotype"] == "Unknown" or c1["genotype"] != c2["genotype"]:
                continue
            if min(c1["length_seq"], c2["length_seq"]) < min_overlap:
                continue

            identity, coverage = blast_pairwise_identity_coverage(
                str(c1["record"].seq), str(c2["record"].seq)
            )
            c1.setdefault("compare_log", []).append(
                f"vs {c2['cluster_label']}: id={identity:.2f}% cov={coverage:.1f}%"
            )

            if identity >= identity_threshold and coverage >= min_coverage:
                uf.union(i, j)

    groups = {}
    for i in range(len(cluster_info)):
        root = uf.find(i)
        groups.setdefault(root, []).append(i)

    merged_cluster_info = []
    for root, members in groups.items():
        if len(members) == 1:
            merged_cluster_info.append(cluster_info[members[0]])
            continue

        group = [cluster_info[m] for m in members]
        total_reads = sum(c["reads"] for c in group)
        genotype = group[0]["genotype"]
        target_len = ref_lengths.get(genotype, None)

        if target_len is not None:
            representative = min(group, key=lambda c: abs(c["length_seq"] - target_len))
        else:
            representative = max(group, key=lambda c: c["length_seq"])

        merged_ids = ",".join(c["cluster_label"] for c in group if c is not representative)

        new_cluster = dict(representative)
        new_cluster["reads"] = total_reads
        new_cluster["merged_from"] = merged_ids
        new_cluster["is_merged"] = True
        merged_cluster_info.append(new_cluster)

    return merged_cluster_info


def main():
    if len(sys.argv) < 5:
        print("Usage: python3 08_genotype_and_filter.py <consensus_dir> <csvs_dir> <blast_db> <cutoff_percent> [ref_fasta_for_lengths]")
        sys.exit(1)

    consensus_dir = sys.argv[1]
    csvs_dir      = sys.argv[2]
    blast_db      = sys.argv[3]
    cutoff        = float(sys.argv[4])
    ref_fasta_for_lengths = sys.argv[5] if len(sys.argv) > 5 else None

    ref_lengths = {}
    if ref_fasta_for_lengths and os.path.exists(ref_fasta_for_lengths):
        ref_lengths = load_ref_lengths(ref_fasta_for_lengths)

    IDENTITY_MERGE_THRESHOLD = 97.0   # % identité sur la région alignée
    MIN_COVERAGE_THRESHOLD   = 85.0   # % de la séquence la plus courte devant être couverte
    MIN_OVERLAP_LENGTH       = 500    # longueur minimale pour tenter une comparaison

    summary_rows = []
    validated_records = []
    fasta_files = sorted(glob.glob(os.path.join(consensus_dir, "*.fasta")))

    for f_path in fasta_files:
        sample_id = os.path.basename(f_path).replace(".fasta", "")
        csv_file = os.path.join(csvs_dir, f"{sample_id}_results.csv")

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
            except Exception as e:
                print(f"⚠️ Erreur lecture comptage réel {csv_file}: {e}")

        blast_out = f"{sample_id}_blast.tsv"
        cmd = (
            f"blastn -query '{f_path}' -db '{blast_db}' "
            f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' "
            f"-max_target_seqs 1 > '{blast_out}'"
        )
        os.system(cmd)

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

        records = list(SeqIO.parse(f_path, "fasta"))
        cluster_info = []
        for idx, rec in enumerate(records):
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

            hit = blast_hits.get(rec.id, {
                "genotype": "Unknown", "pident": 0.0, "length": len(rec.seq),
                "sstart": 1, "send": 2
            })
            strand = "minus" if hit["sstart"] > hit["send"] else "plus"

            cluster_info.append({
                "record": rec,
                "reads": reads,
                "genotype": hit["genotype"],
                "pident": hit["pident"],
                "length": hit["length"],
                "length_seq": len(rec.seq),
                "strand": strand,
                "cluster_label": f"c{idx + 1}",
            })

        if not cluster_info:
            continue

        max_reads_tmp = max([c["reads"] for c in cluster_info]) or 1
        for c in cluster_info:
            ratio_tmp = (c["reads"] / max_reads_tmp) * 100.0
            c["status_tmp"] = "VALIDATED" if (c["pident"] >= 70.0 and ratio_tmp >= cutoff) else "REJECTED"

        cluster_info = merge_split_clusters(
            cluster_info, ref_lengths,
            IDENTITY_MERGE_THRESHOLD, MIN_COVERAGE_THRESHOLD, MIN_OVERLAP_LENGTH
        )

        max_reads = max([c["reads"] for c in cluster_info]) or 1
        for c_idx, c in enumerate(cluster_info, start=1):
            ratio = (c["reads"] / max_reads) * 100.0
            status = "VALIDATED" if (c["pident"] >= 70.0 and ratio >= cutoff) else "REJECTED"
            canonical_cluster_id = f"{sample_id}_c{c_idx}"

            note = ""
            if c.get("is_merged"):
                note = f"Fusion split-artefact (clusters {c.get('merged_from')})"

            summary_rows.append({
                "sample": sample_id,
                "genotype": c["genotype"],
                "best_cluster": canonical_cluster_id,
                "total_reads": c["reads"],
                "ratio_percent": f"{ratio:.2f}",
                "ratio_num": ratio,
                "pident": c["pident"],
                "length": c["length_seq"],
                "strand": c["strand"],
                "status": status,
                "merge_note": note,
            })

            if status == "VALIDATED":
                final_seq = reverse_complement(str(c["record"].seq)) if c["strand"] == "minus" else str(c["record"].seq)
                final_id = f"{sample_id}_c{c_idx}_Geno_{c['genotype']}"
                new_rec = c["record"]
                new_rec.id = final_id
                new_rec.description = ""
                new_rec.seq = Seq(final_seq)
                validated_records.append(new_rec)

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv("SUMMARY_Multi_Infection.tsv", sep="\t", index=False)

    SeqIO.write(validated_records, "validated_consensus_all.fasta", "fasta")

    os.makedirs("renamed_consensus", exist_ok=True)
    for rec in validated_records:
        SeqIO.write([rec], os.path.join("renamed_consensus", f"{rec.id}.fasta"), "fasta")


if __name__ == "__main__":
    main()