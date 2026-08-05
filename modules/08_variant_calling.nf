/*
========================================================================================
    MODULE 08: VARIANT CALLING CLAIR3
========================================================================================
*/

process VARIANT_CALLING_CLAIR3 {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/08_VARIANT_CALLING/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(draft_fasta), path(bam), path(bai)

    output:
    tuple val(sample_id), val(genotype), path("*_clair3.vcf.gz"), optional: true, emit: vcf

    script:
    """
    run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${draft_fasta} \
        --threads=${task.cpus} \
        --platform="ont" \
        --model_path="${params.clair3_model_path}" \
        --output=.

    if [ -f "merge_output.vcf.gz" ]; then
        mv merge_output.vcf.gz ${sample_id}_${genotype}_clair3.vcf.gz
    fi
    """
}