nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 10: CONTRÔLE QUALITÉ CENTRALISÉ (MOSDEPTH + ANALYSE R/PYTHON)
========================================================================================
*/

process QC_MOSDEPTH {
    tag "mosdepth_all_bams"

    input:
    path bams
    path bais

    output:
    path "*_cov*" , emit: cov_files

    script:
    """
    echo "=== 1. Calcul des couvertures via Mosdepth ==="
    for bam in *.sorted.bam; do
        if [ -f "\$bam" ]; then
            prefix=\$(basename "\$bam" .bam)
            mosdepth -n --fast-mode "\${prefix}_cov" "\$bam"
        fi
    done
    """

    stub:
    """
    echo "=== [STUB] Calcul des couvertures via Mosdepth ==="
    for bam in *.sorted.bam; do
        if [ -f "\$bam" ]; then
            prefix=\$(basename "\$bam" .bam)
            touch "\${prefix}_cov.mosdepth.global.dist.txt"
            touch "\${prefix}_cov.mosdepth.summary.txt"
        fi
    done
    """
}

process QC_ANALYSIS {
    publishDir "${params.outdir}/10_QC_ANALYSIS", mode: 'copy'

    input:
    path raw_fastqs
    path dehosted_fastqs
    path trimmed_fastqs
    path bams
    path bais
    path vcfs
    path summary_tsv
    path validated_consensus_fasta
    path cov_files
    val min_length
    val max_length

    output:
    path "RUN_METRICS_SUMMARY_TABLE.tsv"     , emit: metrics_table
    path "RUN_READ_LENGTHS_DISTRIBUTION.pdf" , emit: lengths_pdf, optional: true
    path "matrix_table.csv"                  , emit: pairwise_csv, optional: true
    path "matrix_comp.png"                   , emit: pairwise_png, optional: true

    script:
    """
    # Résolution des droits de cache pour Matplotlib et Fontconfig dans Singularity
    export MPLCONFIGDIR=\$(pwd)/.matplotlib_cache
    export XDG_CACHE_HOME=\$(pwd)/.fontconfig_cache
    mkdir -p \$MPLCONFIGDIR \$XDG_CACHE_HOME

    echo "=== 2. Génération de la table de synthèse QC ==="
    Rscript ${projectDir}/bin/10_qc_analysis.R "${summary_tsv}" "${min_length}" "${max_length}"

    echo "=== 2.bis Génération du rapport PDF multipage (3 graphes par barcode) ==="
    python3 ${projectDir}/bin/10_plot_read_lengths.py \\
        --work-dir . \\
        --min-len "${min_length}" \\
        --out-pdf "RUN_READ_LENGTHS_DISTRIBUTION.pdf"

    echo "=== 3. Génération de la matrice pairwise & Heatmap ==="
    N_SEQ=\$(grep -c "^>" "${validated_consensus_fasta}" || echo 0)
    if [ "\$N_SEQ" -ge 2 ]; then
        mkdir -p distmatrix_out
        python3 ${projectDir}/bin/10_gen_distmatrix.py \\
            --file "${validated_consensus_fasta}" \\
            --out_dir ./distmatrix_out/

        CLUSTERNAME=\$(basename "${validated_consensus_fasta}" .fasta)
        CSV_OUT="distmatrix_out/\${CLUSTERNAME}_pw-brief.csv"
        PNG_OUT="distmatrix_out/distmat_\${CLUSTERNAME}.png"

        [ -f "\$CSV_OUT" ] && mv "\$CSV_OUT" matrix_table.csv
        [ -f "\$PNG_OUT" ] && mv "\$PNG_OUT" matrix_comp.png
    fi
    """

    stub:
    """
    echo "=== [STUB] Génération des livrables QC ==="
    touch "RUN_METRICS_SUMMARY_TABLE.tsv"
    touch "RUN_READ_LENGTHS_DISTRIBUTION.pdf"
    touch "matrix_table.csv"
    touch "matrix_comp.png"
    """
}