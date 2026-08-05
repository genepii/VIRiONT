/*
========================================================================================
    MODULE 14: CALCUL DE COUVERTURE PAR ÉCHANTILLON (10_COVERAGE)
========================================================================================
*/

process COMPUTE_BAM_COVERAGE {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/10_COVERAGE/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(precons_fasta), path(bam), path(bai)

    output:
    path "${sample_id}_${genotype}.cov", emit: sample_cov

    script:
    def sample_geno_id = "${sample_id}_${genotype}"
    """
    bedtools genomecov -ibam ${bam} -d -split | sed "s/\$/\t${sample_id}\tMINION/" > "${sample_geno_id}.cov"
    """
}
