/*
========================================================================================
    MODULE 06: RAPPORT PDF GLOBAL (05_GENOTYPING/GLOBAL_Genotyping_Report.pdf)
========================================================================================
*/

process GLOBAL_REPORT {
    publishDir "${params.outdir}/05_GENOTYPING", mode: 'copy'

    input:
    path validated_tsvs

    output:
    path "GLOBAL_Genotyping_Report.pdf", emit: global_pdf
    path "SUMMARY_Multi_Infection.tsv", emit: summary_tsv

    script:
    """
    # 1. Agrégation du TSV global multi-infection
    echo -e "sample\tgenotype\ttotal_reads\tratio_percent" > SUMMARY_Multi_Infection.tsv
    for file in ${validated_tsvs}; do
        awk -F'\\t' 'NR>1 {
            ratio = \$5; gsub(/%/, "", ratio);
            print \$1"\\t"\$2"\\t"\$4"\\t"ratio
        }' \$file >> SUMMARY_Multi_Infection.tsv
    done

    # 2. Génération du rapport PDF global multi-pages
    Rscript ${projectDir}/bin/generate_global_report.R . ${params.mi_cutoff} "GLOBAL_Genotyping_Report.pdf" || true
    """
}