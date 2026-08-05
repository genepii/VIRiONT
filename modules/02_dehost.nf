/*
========================================================================================
    MODULE 02: DEHOSTING HOSTILE (02_DEHOSTING)
========================================================================================
*/

process DEHOST_HOSTILE {
    tag "$sample_id"
    publishDir "${params.outdir}/02_DEHOSTING", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_dehosted.fastq.gz"), emit: dehosted_fastq

    script:
    """
    hostile clean \
        --fastq1 ${merged_fastq} \
        --index ${params.hostile_index} \
        --aligner minimap2 \
        --threads ${task.cpus} \
        -o .

    if [ -f "${sample_id}_merged.clean.fastq.gz" ]; then
        mv "${sample_id}_merged.clean.fastq.gz" "${sample_id}_dehosted.fastq.gz"
    else
        alt_file=\$(find . -maxdepth 1 -name "*.clean.fastq.gz" | head -n 1)
        if [ -n "\$alt_file" ]; then
            mv "\$alt_file" "${sample_id}_dehosted.fastq.gz"
        else
            echo "❌ Fichier Hostile non généré pour ${sample_id}"
            exit 1
        fi
    fi
    """
}