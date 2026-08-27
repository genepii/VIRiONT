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
    path "matrix_table.csv"              , emit: pairwise_csv, optional: true
    path "matrix_comp.png"               , emit: pairwise_png, optional: true
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

    echo "=== 3. Génération de la matrice pairwise & Heatmap (Python / Biopython) ==="
    N_SEQ=\$(grep -c "^>" "${validated_consensus_fasta}" || echo 0)
    if [ "\$N_SEQ" -ge 2 ]; then
        mkdir -p distmatrix_out
        python3 ${projectDir}/bin/09_gen_distmatrix.py \\
            --file "${validated_consensus_fasta}" \\
            --out_dir ./distmatrix_out/

        # Le script nomme ses sorties d'après le nom du fichier fasta d'entrée
        # (clustername) : on les renomme vers des noms stables pour publishDir.
        CLUSTERNAME=\$(basename "${validated_consensus_fasta}" .fasta)
        CSV_OUT="distmatrix_out/\${CLUSTERNAME}_pw-brief.csv"
        PNG_OUT="distmatrix_out/distmat_\${CLUSTERNAME}.png"

        if [ -f "\$CSV_OUT" ]; then
            mv "\$CSV_OUT" matrix_table.csv
        else
            echo "⚠️ WARNING : CSV de matrice pairwise non trouvé (\$CSV_OUT)"
        fi
        if [ -f "\$PNG_OUT" ]; then
            mv "\$PNG_OUT" matrix_comp.png
        else
            echo "⚠️ WARNING : PNG de heatmap non trouvé (\$PNG_OUT)"
        fi
    else
        echo "⚠️ Moins de 2 séquences validées : matrice pairwise non générée."
    fi
    """
}