nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 11: CALCUL ET TRACÉ DES PROFILS DE COUVERTURE (11_COVERAGE, PAR GÉNOTYPE)
========================================================================================
*/
process COMPUTE_BAM_COVERAGE {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/11_COVERAGE/${sample_id}/${genotype}" }, mode: 'copy'
    input:
    tuple val(sample_id), val(genotype), path(bam), path(bai)
    output:
    path "${sample_id}_${genotype}.cov", emit: sample_cov
    script:
    """
    echo "=== Calcul de la couverture bedtools pour ${sample_id} - Génotype ${genotype} ==="
    if [ -f "${bam}" ] && [ -s "${bam}" ]; then
        bedtools genomecov -ibam "${bam}" -d -split | sed "s/\$/\t${sample_id}\t${genotype}\tMINION/" > "${sample_id}_${genotype}.cov"
    else
        echo "⚠️ WARNING: BAM vide ou inexistant pour ${sample_id} - ${genotype}"
        touch "${sample_id}_${genotype}.cov"
    fi
    """
}

process PLOT_GLOBAL_COVERAGE {
    publishDir "${params.outdir}/11_COVERAGE", mode: 'copy'
    input:
    path cov_files
    output:
    path "VIRiONT_coverage_profiles.pdf", optional: true, emit: coverage_pdf
    path "cov_sum.cov"                 , optional: true, emit: global_cov
    script:
    """
    echo "=== Concaténation des couvertures et tracé des profils R ==="
    cat ${cov_files} > cov_sum.cov
    if [ -s "cov_sum.cov" ]; then
        Rscript ${projectDir}/bin/11_plot_cov_MI.R cov_sum.cov VIRiONT_coverage_profiles.pdf || true
    else
        echo "⚠️ WARNING: Aucun profil de couverture à tracer."
    fi
    """
}