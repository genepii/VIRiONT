nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 09: CONTROL QUALITÉ CENTRALISÉ (09_QC_ANALYSIS)
========================================================================================
*/
process QC_ANALYSIS {
    publishDir "${params.outdir}/09_QC_ANALYSIS", mode: 'copy'

    input:
    path raw_fastqs
    path dehosted_fastqs
    path trimmed_fastqs
    path bams
    path bais
    path vcfs
    path summary_tsv
    path renamed_consensus
    val min_length
    val max_length

    output:
    path "QC_Metrics_Summary.tsv", emit: qc_table
    path "matrix_table.tsv"      , emit: matrix_table, optional: true
    path "matrix_comp.pdf"       , emit: matrix_pdf  , optional: true
    path "*.mosdepth.*"          , optional: true

    script:
    """
    echo "=== 1. Calcul des couvertures via Mosdepth ==="
    for bam in ${bams}; do
        if [ -f "\$bam" ]; then
            prefix=\$(basename "\$bam" .bam)
            mosdepth -n --fast-mode "\${prefix}_cov" "\$bam"
        fi
    done

    echo "=== 2. Génération de la table de synthèse QC ==="
    Rscript ${projectDir}/bin/09_qc_analysis.R "${summary_tsv}" "${min_length}" "${max_length}"

    echo "=== 3. Génération de la matrice pairwise & Heatmap ==="
    Rscript ${projectDir}/bin/09_pairwise_matrix.R
    """
}