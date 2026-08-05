/*
========================================================================================
    MODULE 03: CHOPPER (03_FILTERED_TRIMMED)
========================================================================================
*/

process TRIM_CHOPPER {
    tag "$sample_id"
    publishDir "${params.outdir}/03_FILTERED_TRIMMED", mode: 'copy'

    input:
    tuple val(sample_id), path(dehosted_fastq)
    val min_length
    val max_length

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_fastq

    script:
    """
    gunzip -c ${dehosted_fastq} | \
    chopper -l ${min_length} --maxlength ${max_length} | \
    gzip -c > ${sample_id}_trimmed.fastq.gz
    """
}