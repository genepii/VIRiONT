/*
========================================================================================
    MODULE 15: TRACÉ DES PROFILS DE COUVERTURE GLOBALE (GGPLOT2)
========================================================================================
*/

process PLOT_GLOBAL_COVERAGE {
    publishDir "${params.outdir}/08_QC_ANALYSIS", mode: 'copy'

    input:
    path cov_files

    output:
    path "VIRiONT_coverage_profiles.pdf", optional: true, emit: coverage_pdf

    script:
    """
    cat ${cov_files} > all_samples_coverage.cov
    15_plot_cov_MI.R all_samples_coverage.cov VIRiONT_coverage_profiles.pdf
    """
}