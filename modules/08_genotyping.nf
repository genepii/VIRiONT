nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 08: GENOTYPING BLAST & GLOBAL REPORT (08_GENOTYPING)
========================================================================================
*/
process GENOTYPING {
    publishDir "${params.outdir}/08_GENOTYPING", mode: 'copy', pattern: '{GLOBAL_Genotyping_Report.pdf,SUMMARY_Multi_Infection.tsv,*_blastnR.tsv,*_fmt0.txt,renamed_consensus/*.fasta}'
    publishDir "${params.outdir}/00_SUPDATA", mode: 'copy', pattern: 'DB/*'

    input:
    path consensus_fastas  // Fichiers consensus émis par Medaka
    path ref_fasta         // FASTA de référence
    path 'csvs/*'          // CSV d'AmpliconSorter re-nommés par main.nf

    output:
    path "GLOBAL_Genotyping_Report.pdf"      , emit: global_pdf
    path "SUMMARY_Multi_Infection.tsv"        , emit: summary_tsv
    path "*_blastnR.tsv"                      , emit: sample_blast_tsv
    path "*_fmt0.txt"                         , emit: sample_blast_fmt0
    path "renamed_consensus/*.fasta"          , emit: genotyped_fastas
    path "DB/*"                               , emit: db_files

    script:
    def base_name = ref_fasta.baseName
    """
    echo "=== 1. Indexation de la base de référence (${base_name}) ==="

    mkdir -p DB
    mkdir -p renamed_consensus

    if [ ! -f "DB/${base_name}.nhr" ]; then
        makeblastdb -in "${ref_fasta}" -dbtype nucl -out "DB/${base_name}"
    fi

    echo -e "sample\tgenotype\tbest_cluster\ttotal_reads\tratio_percent\tratio_num\tpident\tlength\tstrand\tstatus" > SUMMARY_Multi_Infection.tsv

    echo "=== 2. Génotypage & Extraction des effectifs réels ==="

    for fq in ${consensus_fastas}; do
        if [ -s "\$fq" ]; then
            sample_id=\$(basename "\$fq" .fasta | sed 's/_consensus//g')

            # Ciblage direct du CSV renommé transmis par le workflow
            csv_file="csvs/\${sample_id}_results.csv"

            echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > \${sample_id}_blastnR.tsv

            # BLAST Tabulaire
            blastn \\
                -query "\$fq" \\
                -db "DB/${base_name}" \\
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
                -max_target_seqs 5 \\
                -num_threads ${task.cpus} > tmp_hits.txt

            cat tmp_hits.txt >> \${sample_id}_blastnR.tsv

            # BLAST Pairwise
            echo "==========================================================================" > \${sample_id}_fmt0.txt
            echo "BLAST against database: ${base_name}" >> \${sample_id}_fmt0.txt
            echo "==========================================================================" >> \${sample_id}_fmt0.txt
            blastn \\
                -query "\$fq" \\
                -db "DB/${base_name}" \\
                -outfmt 0 \\
                -max_target_seqs 3 \\
                -num_threads ${task.cpus} >> \${sample_id}_fmt0.txt

            if [ -s tmp_hits.txt ]; then
                awk -v sample="\$sample_id" -v cutoff="${params.mi_cutoff}" -v csv="\$csv_file" '
                BEGIN {
                    OFS = "\t";
                    n_hits = 0;
                    c_idx = 0;
                    
                    # Lecture directe et robuste du CSV
                    if (csv != "") {
                        while ((getline line < csv) > 0) {
                            split(line, fields, ",");
                            c_name = fields[1];
                            gsub(/^[ \t]+|[ \t]+\$/, "", c_name);
                            
                            c_val = fields[2];
                            gsub(/^[ \t]+|[ \t]+\$/, "", c_val);
                            val_num = c_val + 0;
                            
                            if (c_name != "" && c_name != "Total" && val_num > 0) {
                                c_idx++;
                                csv_reads[c_idx] = val_num;
                            }
                        }
                        close(csv);
                    }
                }
                {
                    contig = \$1;
                    target_genotype = \$2;
                    pident = \$3;
                    len = \$4;
                    sstart = \$9 + 0;
                    send = \$10 + 0;

                    if (!seen_contig[contig]++) {
                        n_hits++;
                        
                        # Récupération ordonnée directe par index de cluster
                        reads_count = (n_hits in csv_reads) ? csv_reads[n_hits] : 1;

                        contigs[n_hits] = contig;
                        genotypes[n_hits] = target_genotype;
                        read_counts[n_hits] = reads_count;
                        pidents[n_hits] = pident;
                        lengths[n_hits] = len;
                        strands[n_hits] = (sstart > send) ? "minus" : "plus";
                    }
                }
                END {
                    max_reads = 0;
                    for (i = 1; i <= n_hits; i++) {
                        if (read_counts[i] > max_reads) max_reads = read_counts[i];
                    }
                    if (max_reads == 0) max_reads = 1;

                    for (i = 1; i <= n_hits; i++) {
                        ratio_val = (read_counts[i] / max_reads) * 100;
                        ratio = sprintf("%.2f", ratio_val);
                        status = (pidents[i] >= 70.0 && ratio_val >= cutoff) ? "VALIDATED" : "REJECTED";

                        print sample, genotypes[i], contigs[i], read_counts[i], ratio, ratio, pidents[i], lengths[i], strands[i], status;
                    }
                }' tmp_hits.txt >> SUMMARY_Multi_Infection.tsv

                # Éclatement des contigs
                rm -f tmp_seq_*
                csplit -s -z -f tmp_seq_ "\$fq" '/^>/' '{*}'

                for sf in tmp_seq_*; do
                    if [ -s "\$sf" ]; then
                        seq_hdr=\$(head -n 1 "\$sf" | cut -d' ' -f1 | sed 's/>//g')
                        
                        hit_line=\$(grep -w "\$seq_hdr" tmp_hits.txt | head -n 1)
                        if [ -z "\$hit_line" ]; then
                            hit_line=\$(head -n 1 tmp_hits.txt)
                        fi

                        sub_geno=\$(echo "\$hit_line" | awk '{print \$2}')
                        sstart=\$(echo "\$hit_line" | awk '{print \$9}')
                        send=\$(echo "\$hit_line" | awk '{print \$10}')

                        clean_hdr=">\${seq_hdr}_Geno_\${sub_geno}"

                        if [ "\$sstart" -gt "\$send" ]; then
                            echo "🔄 Re-orientation 5'->3' pour \${seq_hdr}"
                            seqkit seq --reverse --complement "\$sf" | awk -v h="\$clean_hdr" '/^>/ {print h; next} {print}' > "renamed_consensus/\${seq_hdr}_Geno_\${sub_geno}.fasta"
                        else
                            awk -v h="\$clean_hdr" '/^>/ {print h; next} {print}' "\$sf" > "renamed_consensus/\${seq_hdr}_Geno_\${sub_geno}.fasta"
                        fi

                        rm -f "\$sf"
                    fi
                done

            else
                echo -e "\${sample_id}\tNo_Match\tNone\t0\t0\t0\t0\t0\tplus\tREJECTED" >> SUMMARY_Multi_Infection.tsv
                clean_hdr=">\${sample_id}_Geno_Unknown"
                awk -v h="\$clean_hdr" '/^>/ {print h; next} {print}' "\$fq" > "renamed_consensus/\${sample_id}_Geno_Unknown.fasta"
            fi

            rm -f tmp_hits.txt
        fi
    done

    echo "=== 3. Génération du Rapport PDF Global ==="
    Rscript ${projectDir}/bin/08_generate_report.R SUMMARY_Multi_Infection.tsv ${params.mi_cutoff} GLOBAL_Genotyping_Report.pdf
    """
}