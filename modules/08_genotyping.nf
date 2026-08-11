nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 08: GENOTYPING BLAST & GLOBAL REPORT (08_GENOTYPING)
========================================================================================
*/
process GENOTYPING {
    publishDir "${params.outdir}/08_GENOTYPING", mode: 'copy'

    input:
    path consensus_fastas  // Fichiers consensus collectés
    path ref_fasta         // Fichier FASTA de référence ciblé

    output:
    path "GLOBAL_Genotyping_Report.pdf"      , emit: global_pdf
    path "SUMMARY_Multi_Infection.tsv"        , emit: summary_tsv
    path "*_blastnR.tsv"                      , emit: sample_blast_tsv
    path "*_fmt0.txt"                         , emit: sample_blast_fmt0

    script:
    def base_name = ref_fasta.baseName
    """
    echo "=== 1. Indexation de la base de référence ciblée (${base_name}) ==="

    DB_DIR="${projectDir}/00_SUPDATA/DB"
    mkdir -p \$DB_DIR

    if [ ! -f "\$DB_DIR/${base_name}.nhr" ]; then
        echo "⚙️ Génération de la base BLAST pour ${base_name}..."
        makeblastdb -in "${ref_fasta}" -dbtype nucl -out "\$DB_DIR/${base_name}"
    fi

    # En-tête exact attendu par 06_generate_report.R
    echo -e "sample\tgenotype\tbest_cluster\ttotal_reads\tratio_percent\tratio_num\tpident\tlength\tstrand\tstatus" > SUMMARY_Multi_Infection.tsv

    echo "=== 2. Génotypage & Alignements Détaillés ==="

    for fq in ${consensus_fastas}; do
        sample_id=\$(basename "\$fq" _consensus.fasta)
        sample_id=\$(echo "\$sample_id" | sed 's/.fasta//g')

        # En-tête fichier BLAST individuel
        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > \${sample_id}_blastnR.tsv

        # A. BLAST Tabulaire (fmt 6)
        blastn \\
            -query "\$fq" \\
            -db "\$DB_DIR/${base_name}" \\
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
            -max_target_seqs 5 \\
            -num_threads ${task.cpus} > tmp_hits.txt

        cat tmp_hits.txt >> \${sample_id}_blastnR.tsv

        # B. BLAST Pairwise (fmt 0)
        echo "==========================================================================" > \${sample_id}_fmt0.txt
        echo "BLAST against database: ${base_name}" >> \${sample_id}_fmt0.txt
        echo "==========================================================================" >> \${sample_id}_fmt0.txt
        blastn \\
            -query "\$fq" \\
            -db "\$DB_DIR/${base_name}" \\
            -outfmt 0 \\
            -max_target_seqs 3 \\
            -num_threads ${task.cpus} >> \${sample_id}_fmt0.txt

        # C. Identification du Majeur (100%) et calcul du ratio vs Seuil Clinique 30%
        if [ -s tmp_hits.txt ]; then
            awk -v sample="\$sample_id" -v cutoff="${params.mi_cutoff}" '
            BEGIN {
                OFS = "\t";
                n_hits = 0;
            }
            {
                contig = \$1;
                target_genotype = \$2;
                pident = \$3;
                len = \$4;

                # Garde uniquement le meilleur hit BLAST pour chaque contig
                if (!seen_contig[contig]++) {
                    split(contig, parts, "[()]");
                    reads = (length(parts[2]) > 0) ? parts[2] + 0 : len + 0;

                    n_hits++;
                    contigs[n_hits] = contig;
                    genotypes[n_hits] = target_genotype;
                    read_counts[n_hits] = reads;
                    pidents[n_hits] = pident;
                    lengths[n_hits] = len;
                }
            }
            END {
                # 1. Nombre de reads du Génotype Majeur
                max_reads = 0;
                for (i = 1; i <= n_hits; i++) {
                    if (read_counts[i] > max_reads) {
                        max_reads = read_counts[i];
                    }
                }
                if (max_reads == 0) max_reads = 1;

                # 2. Ratio relatif par rapport au Majeur (%) & Statut (Validation si ratio >= 30%)
                for (i = 1; i <= n_hits; i++) {
                    ratio_val = (read_counts[i] / max_reads) * 100;
                    ratio = sprintf("%.2f", ratio_val);

                    status = (pidents[i] >= 70.0 && ratio_val >= cutoff) ? "VALIDATED" : "REJECTED";

                    print sample, genotypes[i], contigs[i], read_counts[i], ratio, ratio, pidents[i], lengths[i], "plus", status;
                }
            }' tmp_hits.txt >> SUMMARY_Multi_Infection.tsv
        else
            echo -e "\${sample_id}\tNo_Match\tNone\t0\t0\t0\t0\t0\tplus\tREJECTED" >> SUMMARY_Multi_Infection.tsv
        fi

        rm -f tmp_hits.txt
    done

    echo "=== 3. Génération du Rapport PDF Global ==="
    Rscript ${projectDir}/bin/06_generate_report.R SUMMARY_Multi_Infection.tsv ${params.mi_cutoff} GLOBAL_Genotyping_Report.pdf
    """
}