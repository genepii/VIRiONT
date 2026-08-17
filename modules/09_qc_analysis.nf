nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 09: CONTRÔLE QUALITÉ CENTRALISÉ (09_QC_ANALYSIS)
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
    path validated_consensus_fasta
    val min_length
    val max_length

    output:
    path "RUN_METRICS_SUMMARY_TABLE.tsv", emit: metrics_table
    path "matrix_table.tsv"              , emit: pairwise_tsv, optional: true
    path "matrix_comp.pdf"               , emit: pairwise_pdf, optional: true

    script:
    """
    echo "=== 1. Calcul des couvertures via Mosdepth ==="
    for bam in *.sorted.bam; do
        if [ -f "\$bam" ]; then
            prefix=\$(basename "\$bam" .bam)
            mosdepth -n --fast-mode "\${prefix}_cov" "\$bam"
        fi
    done
    
    echo "=== 2. Génération de la table de synthèse QC ==="
    Rscript ${projectDir}/bin/09_qc_analysis.R "${summary_tsv}" "${min_length}" "${max_length}"
    
    echo "=== 3. Génération de la matrice pairwise & Heatmap ==="
    Rscript ${projectDir}/bin/09_pairwise_matrix.R "${validated_consensus_fasta}"
    """
}