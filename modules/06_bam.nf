nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 06: ALIGNMENT & BAM GENERATION (06_BAM)
========================================================================================
*/
process ALIGN_BAM {
    tag "$sample_id"
    // Remplacement des guillemets par une closure { }
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