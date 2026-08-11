nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 09: QC ANALYSIS (MOSDEPTH & METRICS SUMMARY)
========================================================================================
*/
process QC_ANALYSIS {
    publishDir "${params.outdir}/09_QC_ANALYSIS", mode: 'copy'

    input:
    path raw_fastqs        // FASTQ bruts fusionnés (Step 01)
    path dehosted_fastqs   // FASTQ déshôtés (Step 02)
    path trimmed_fastqs    // FASTQ filtrés (Step 03)
    path bams              // BAM alignés (Step 06)
    path bais              // BAI index (Step 06)
    path vcfs              // VCF Clair3 (Step 07)
    path genotyping_tsv    // SUMMARY_Multi_Infection.tsv (Step 08)

    output:
    path "METRIC_summary_table.csv" , emit: metrics_csv
    path "ALIGNMENT_QC/*"            , emit: alignment_qc

    script:
    """
    echo "=== 1. Execution de Mosdepth (ALIGNMENT_QC) ==="
    
    for bam in ${bams}; do
        if [ -f "\$bam" ]; then
            sample_id=\$(basename "\$bam" .bam)
            mkdir -p "ALIGNMENT_QC/\${sample_id}"

            mosdepth \\
                --threads ${task.cpus} \\
                --no-per-base \\
                --by 50 \\
                "ALIGNMENT_QC/\${sample_id}/\${sample_id}" \\
                "\$bam" &>/dev/null || true
        fi
    done

    echo "=== 2. Génération de la table METRIC_summary_table.csv ==="

    HEADER="sample;step;assignedref;read_count;minlengthread;maxlengthread;meanread_length;medianread_length;pident_blast;mean_depth_coverage;clair3_variants;assigned_reads"
    echo "\$HEADER" > METRIC_summary_table.csv

    # A. Stats 01_RAW
    for fq in ${raw_fastqs}; do
        s_id=\$(basename "\$fq" | sed -E 's/(_merged|_dehosted|_trimmed)?\\.fastq\\.gz//g')
        if [ -s "\$fq" ]; then
            st=\$(seqkit stats -T "\$fq" | tail -n 1)
            rc=\$(echo "\$st" | awk '{print \$4}')
            min_l=\$(echo "\$st" | awk '{print \$5}')
            mean_l=\$(echo "\$st" | awk '{print \$6}')
            max_l=\$(echo "\$st" | awk '{print \$7}')
            med_l=\$(echo "\$st" | awk '{print \$13}')
            echo "\${s_id};01_RAW;NONE;\${rc};\${min_l};\${max_l};\${mean_l};\${med_l};NA;NA;NA;NA" >> METRIC_summary_table.csv
        fi
    done

    # B. Stats 02_DEHOSTING
    for fq in ${dehosted_fastqs}; do
        s_id=\$(basename "\$fq" | sed -E 's/(_merged|_dehosted|_trimmed)?\\.fastq\\.gz//g')
        if [ -s "\$fq" ]; then
            st=\$(seqkit stats -T "\$fq" | tail -n 1)
            rc=\$(echo "\$st" | awk '{print \$4}')
            min_l=\$(echo "\$st" | awk '{print \$5}')
            mean_l=\$(echo "\$st" | awk '{print \$6}')
            max_l=\$(echo "\$st" | awk '{print \$7}')
            med_l=\$(echo "\$st" | awk '{print \$13}')
            echo "\${s_id};02_DEHOSTING;NONE;\${rc};\${min_l};\${max_l};\${mean_l};\${med_l};NA;NA;NA;NA" >> METRIC_summary_table.csv
        fi
    done

    # C. Stats 03_FILTERED_TRIMMED
    for fq in ${trimmed_fastqs}; do
        s_id=\$(basename "\$fq" | sed -E 's/(_merged|_dehosted|_trimmed)?\\.fastq\\.gz//g')
        if [ -s "\$fq" ]; then
            st=\$(seqkit stats -T "\$fq" | tail -n 1)
            rc=\$(echo "\$st" | awk '{print \$4}')
            min_l=\$(echo "\$st" | awk '{print \$5}')
            mean_l=\$(echo "\$st" | awk '{print \$6}')
            max_l=\$(echo "\$st" | awk '{print \$7}')
            med_l=\$(echo "\$st" | awk '{print \$13}')
            echo "\${s_id};03_FILTERED_TRIMMED;NONE;\${rc};\${min_l};\${max_l};\${mean_l};\${med_l};NA;NA;NA;NA" >> METRIC_summary_table.csv
        fi
    done

    # D. Stats 05_GENOTYPING
    if [ -s "${genotyping_tsv}" ]; then
        awk -F'\t' 'NR>1 {
            sample = \$1
            geno = \$2
            reads = \$4
            pident = \$7
            
            depth = "NA"
            mos_summary = "ALIGNMENT_QC/" sample "/" sample ".mosdepth.summary.txt"
            if (system("[ -f " mos_summary " ]") == 0) {
                cmd_d = "tail -n 1 " mos_summary " | awk \'{print \$4}\'"
                cmd_d | getline depth
                close(cmd_d)
            }

            n_var = "0"
            vcf_f = sample ".vcf.gz"
            if (system("[ -f " vcf_f " ]") == 0) {
                cmd_v = "bcftools view -H " vcf_f " 2>/dev/null | wc -l"
                cmd_v | getline n_var
                close(cmd_v)
            }

            printf "%s;05_GENOTYPING;%s;%s;1000;5000;NA;NA;%.2f;%s;%s;NA\n", sample, geno, reads, pident, depth, n_var
        }' "${genotyping_tsv}" >> METRIC_summary_table.csv
    fi
    """
}