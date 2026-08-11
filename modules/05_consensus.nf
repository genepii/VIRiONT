nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 05: MEDAKA POLISHING (05_CONSENSUS)
========================================================================================
*/
process MEDAKA_CONSENSUS {
    tag "$sample_id"
    publishDir { "${params.outdir}/05_CONSENSUS/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq), path(preconsensus_fasta)
    val medaka_model

    output:
    tuple val(sample_id), path("${sample_id}_consensus.fasta"), emit: final_consensus
    tuple val(sample_id), path(trimmed_fastq), path("${sample_id}_consensus.fasta"), emit: for_bam

    script:
    """
    echo "=== Polissage Medaka pour ${sample_id} ==="

    medaka_consensus \
        -i ${trimmed_fastq} \
        -d ${preconsensus_fasta} \
        -o medaka_out \
        -m ${medaka_model} \
        -t ${task.cpus}

    mv medaka_out/consensus.fasta ${sample_id}_consensus.fasta
    """
}