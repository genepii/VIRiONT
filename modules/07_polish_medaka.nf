/*
========================================================================================
    MODULE 07: POLISSAGE MEDAKA
========================================================================================
*/

process POLISH_MEDAKA {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/07_POLISH_MEDAKA/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), val(cluster_id), path(draft_fasta), path(reads_fastq)
    path ref_seq

    output:
    tuple val(sample_id), val(genotype), path("*_polished_consensus.fasta"), optional: true, emit: final_consensus
    tuple val(sample_id), val(genotype), path(draft_fasta), path("calls_to_ref.bam"), path("calls_to_ref.bam.bai"), optional: true, emit: bam_tuple

    script:
    """
    # 1. Alignement des reads sur la référence/draft avec minimap2
    minimap2 -ax map-ont -t ${task.cpus} ${draft_fasta} ${reads_fastq} | samtools sort -o calls_to_ref.bam -
    samtools index calls_to_ref.bam

    # 2. Polissage Medaka
    medaka consensus calls_to_ref.bam consensus.hdf5 --model ${params.medaka_model} --threads ${task.cpus}
    medaka stitch consensus.hdf5 ${draft_fasta} ${sample_id}_${genotype}_polished_consensus.fasta
    """
}