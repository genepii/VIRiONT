nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 07: FILTRAGE CLINIQUE & GÉNOTYPAGE (07_CLINICAL_FILTER)
========================================================================================
*/

process CLINICAL_FILTER {
    publishDir "${params.outdir}/07_GENOTYPING", mode: 'copy'

    input:
    path consensus_fastas
    path trimmed_fastqs
    path ref_db
    path count_tsvs
    val virus_name
    val tech_name

    output:
    path "SUMMARY_Multi_Infection.tsv"        , emit: summary_tsv
    path "validated_consensus_all.fasta"      , emit: validated_all_fasta
    path "renamed_consensus/*.fasta"          , emit: genotyped_fastas
    path "renamed_consensus/manifest.tsv"      , emit: genotype_manifest
    path "Rapport_Repartition_References.pdf" , emit: report_pdf

    script:
    def cutoff = params.mi_cutoff ?: 30.0
    """
    # 1. Analyse et génération du SUMMARY_Multi_Infection.tsv
    07_genotype_and_filter.py \\
        "." \\
        "." \\
        "." \\
        "${ref_db}" \\
        ${cutoff} \\
        "${virus_name}" \\
        "${tech_name}"

    # 2. Rapport graphique combinant profils complets et tableau simple (>= 5 reads)
    07_generate_report.R \\
        "." \\
        ${cutoff} \\
        "Rapport_Repartition_References.pdf"
    """
}