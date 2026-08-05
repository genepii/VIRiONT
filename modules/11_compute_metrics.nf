/*
========================================================================================
    MODULE 11: CALCUL DU TABLEAU RÉCAPITULATIF METROLOGIQUE (08_QC_ANALYSIS)
========================================================================================
*/

process COMPUTE_QC_METRICS {
    publishDir "${params.outdir}/08_QC_ANALYSIS", mode: 'copy'

    input:
    path raw_csvs
    path dehost_csvs
    path trimm_csvs
    path geno_csvs

    output:
    path "METRIC_summary_table.csv", emit: summary_table

    script:
    """
    ALL_RAW="all_rawcount.csv"
    ALL_DEHOST="all_dehostcount.csv"
    ALL_TRIMM="all_trimcount.csv"
    ALL_GENO="all_genocount.csv"

    cat ${raw_csvs} > "\$ALL_RAW" 2>/dev/null || touch "\$ALL_RAW"
    cat ${dehost_csvs} > "\$ALL_DEHOST" 2>/dev/null || touch "\$ALL_DEHOST"
    cat ${trimm_csvs} > "\$ALL_TRIMM" 2>/dev/null || touch "\$ALL_TRIMM"
    cat ${geno_csvs} > "\$ALL_GENO" 2>/dev/null || touch "\$ALL_GENO"

    Rscript ${projectDir}/bin/11_compute_metrics.R "\$ALL_RAW" "\$ALL_DEHOST" "\$ALL_TRIMM" "\$ALL_GENO" "METRIC_summary_table.csv"
    """
}