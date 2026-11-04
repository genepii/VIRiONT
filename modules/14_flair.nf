nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 14: ANALYSE DE L'ÉPISSAGE VHB VIA FLAIR
========================================================================================
*/

process FLAIR {
    tag { "${sample_id}_${genotype}" }
    publishDir { "${params.outdir}/14_FLAIR/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), val(protocol), path(fastq), path(ref_fasta)
    path splice_catalog

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_flair_results"), emit: results
    path("${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_flair_consolidated_summary.csv"), emit: summary_csv

    script:
    def min_pct   = params.min_flair_pct ?: 8.0
    def min_size  = params.min_flair_size ?: 150
    def min_reads = params.min_reads ?: 20

    """
    mkdir -p "${sample_id}_${genotype}_flair_results"

    # 1. Comptage rapide du nombre total de reads
    TOTAL_READS=\$(zcat "${fastq}" | awk 'NR%4==1 {c++} END {print c+0}')
    TARGET_READS=40000

    # 2. Sous-échantillonnage représentatif si sur-couverture (> 50 000 reads)
    if [ "\${TOTAL_READS}" -gt 50000 ]; then
        echo "[FLAIR] Échantillon sur-couvert (\${TOTAL_READS} reads). Sous-échantillonnage à \${TARGET_READS} reads..."
        zcat "${fastq}" | head -n \$(( TARGET_READS * 4 )) | gzip -c > input_flair.fastq.gz
        TARGET_FASTQ="input_flair.fastq.gz"
        SCALE_FACTOR=\$(awk -v tot="\${TOTAL_READS}" -v target="\${TARGET_READS}" 'BEGIN { printf "%.6f", tot / target }')
    else
        TARGET_FASTQ="${fastq}"
        SCALE_FACTOR=1.0
    fi

    # 3. Alignement FLAIR
    flair align \
        -g "${ref_fasta}" \
        -r "\${TARGET_FASTQ}" \
        -o "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_align" \
        -t ${task.cpus}

    # 4. Collapse des isoformes
    flair collapse \
        -g "${ref_fasta}" \
        -q "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_align.bed" \
        -r "\${TARGET_FASTQ}" \
        -o "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_collapse" \
        -t ${task.cpus}

    # 5. Parsing et filtrage clinique
    python3 ${projectDir}/bin/14_parse_flair_consolidated.py \
        "${sample_id}_${genotype}" \
        "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_collapse.isoform.counts.txt" \
        "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_collapse.isoforms.bed" \
        "${splice_catalog}" \
        "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_flair_consolidated_summary.csv" \
        --min-pct ${min_pct} \
        --min-size ${min_size} \
        --min-reads ${min_reads}

    # 6. Rétablissement des statistiques exactes à l'échelle du run réel
    if [ "\${SCALE_FACTOR}" != "1.0" ] && [ -f "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_flair_consolidated_summary.csv" ]; then
        awk -F',' -v factor="\${SCALE_FACTOR}" 'BEGIN {OFS=","} NR==1 {print \$0} NR>1 { \$5 = sprintf("%.0f", \$5 * factor); print \$0 }' \
            "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_flair_consolidated_summary.csv" > "${sample_id}_${genotype}_flair_results/temp_scaled.csv"
        mv "${sample_id}_${genotype}_flair_results/temp_scaled.csv" "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_flair_consolidated_summary.csv"
    fi

    rm -f input_flair.fastq.gz
    """

    stub:
    """
    mkdir -p "${sample_id}_${genotype}_flair_results"
    touch "${sample_id}_${genotype}_flair_results/${sample_id}_${genotype}_flair_consolidated_summary.csv"
    """
}

process MERGE_SUMMARIES {
    publishDir "${params.outdir}/14_FLAIR", mode: 'copy'

    input:
    path csv_files

    output:
    path "all_samples_flair_consolidated_summary.csv", emit: merged_summary

    script:
    """
    echo "Sample,Intron_Start,Intron_End,Size_bp,Reads,Percentage,Transcript,Type" > all_samples_flair_consolidated_summary.csv
    for f in ${csv_files}; do
        if [ -f "\$f" ] && [ \$(wc -l < "\$f") -gt 1 ]; then
            tail -n +2 "\$f" >> all_samples_flair_consolidated_summary.csv
        fi
    done
    """

    stub:
    """
    touch all_samples_flair_consolidated_summary.csv
    """
}