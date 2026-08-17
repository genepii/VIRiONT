nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 06: ALIGNMENT & BAM GENERATION (06_BAM)
========================================================================================
*/
process ALIGN_BAM {
    tag "$sample_id"
    publishDir { "${params.outdir}/06_BAM/${sample_id}" }, mode: 'copy'
    input:
    tuple val(sample_id), path(trimmed_fastq), path(final_consensus)
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam_bai
    tuple val(sample_id), path(final_consensus), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: for_vcf
    script:
    """
    echo "=== Alignement minimap2 pour ${sample_id} ==="
    # Indexation du consensus
    samtools faidx ${final_consensus}
    # Alignement Nanopore map-ont + tri + indexation
    minimap2 -ax map-ont -t ${task.cpus} ${final_consensus} ${trimmed_fastq} | \
        samtools view -bS -F 4 - | \
        samtools sort -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """
}

process COUNT_REAL_READS {
    tag "$sample_id"
    publishDir { "${params.outdir}/06_BAM/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_real_counts.tsv"), emit: real_counts

    script:
    """
    echo "=== Comptage exhaustif des reads réels par cluster pour ${sample_id} (tous les reads du barcode, sans sous-échantillonnage) ==="
    echo -e "cluster_id\\treal_read_count" > ${sample_id}_real_counts.tsv
    samtools idxstats ${bam} | awk 'BEGIN{OFS="\\t"} \$1 != "*" && \$3 > 0 {print \$1, \$3}' >> ${sample_id}_real_counts.tsv
    """
}