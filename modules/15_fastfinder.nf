nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 15: EXPORT FASTFINDER / GLIMS (15_FASTFINDER)
========================================================================================
*/

process FASTFINDER {
    tag "$run_id"
    publishDir "${params.outdir}/15_FASTFINDER", mode: 'copy'

    input:
    path summary_tsv
    path consensus_fasta
    val run_id

    output:
    path "fastfinder_export_${run_id}.csv", emit: csv

    script:
    """
    python3 ${projectDir}/bin/generate_fastfinder_csv.py \\
        --summary ${summary_tsv} \\
        --fasta ${consensus_fasta} \\
        --run-id "${run_id}" \\
        --output "fastfinder_export_${run_id}.csv"
    """
}