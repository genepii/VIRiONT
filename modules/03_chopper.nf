nextflow.enable.dsl=2

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
    echo "=== Filtrage Chopper (Q12) pour ${sample_id} [${min_length} bp - ${max_length} bp] ==="

    gunzip -c ${dehosted_fastq} | \
    chopper -q 12 -l ${min_length} --maxlength ${max_length} | \
    gzip -c > ${sample_id}_trimmed.fastq.gz

    # Si le fichier est vide ou quasi-vide, on émet un warning sans faire crasher Nextflow (exit 0)
    if [ ! -s "${sample_id}_trimmed.fastq.gz" ] || [ \$(stat -c%s "${sample_id}_trimmed.fastq.gz") -lt 50 ]; then
        echo "⚠️ WARNING : Aucun read conservé par Chopper pour ${sample_id}."
    fi
    """
}