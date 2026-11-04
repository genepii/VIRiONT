nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 03b: ROGNAGE DES AMORCES AVEC PORECHOP_ABI (03b_PRIMERS_TRIMMED)
========================================================================================
*/

process TRIM_PRIMERS {
    tag "$sample_id"
    publishDir "${params.outdir}/03b_PRIMERS_TRIMMED", mode: 'copy'

    input:
    tuple val(sample_id), path(chopper_fastq)
    path primers_fasta

    output:
    tuple val(sample_id), path("${sample_id}_primertrimmed.fastq.gz"), emit: fastq
    tuple val(sample_id), path("${sample_id}_porechop.log")          , emit: log

    script:
    """
    echo "=== Rognage des amorces avec Porechop_ABI pour ${sample_id} ==="

    porechop_abi \\
        -i "${chopper_fastq}" \\
        -o "${sample_id}_primertrimmed.fastq.gz" \\
        -cap "${primers_fasta}" \\
        -ddb \\
        --no_split \\
        --extra_end_trim 0 \\
        -t ${task.cpus} \\
        -v 2 > "${sample_id}_porechop.log"
    """

    stub:
    """
    touch "${sample_id}_primertrimmed.fastq.gz"
    echo "=== Porechop_ABI stub ===" > "${sample_id}_porechop.log"
    """
}