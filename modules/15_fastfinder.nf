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
    path glims_mutations
    path glims_indels
    val run_id

    output:
    path "fastfinder_export_${run_id}.csv", emit: csv

    script:
    def mut_opt = (glims_mutations && glims_mutations.name != 'NO_FILE' && glims_mutations.exists()) ? "--mutations ${glims_mutations}" : ""
    def indel_opt = (glims_indels && glims_indels.name != 'NO_FILE' && glims_indels.exists()) ? "--indels ${glims_indels}" : ""
    """
    python3 ${projectDir}/bin/15_generate_fastfinder_csv.py \\
        --summary "${summary_tsv}" \\
        --fasta "${consensus_fasta}" \\
        ${mut_opt} \\
        ${indel_opt} \\
        --run-id "${run_id}" \\
        --output "fastfinder_export_${run_id}.csv"
    """

    stub:
    """
    touch "fastfinder_export_${run_id}.csv"
    """
}