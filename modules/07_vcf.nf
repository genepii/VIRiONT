nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 07: VARIANT CALLING CLAIR3 (07_VCF - COORDONNÉES CANONIQUES)
========================================================================================
*/
process CALL_VCF {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/07_VCF/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(canonical_ref), path(bam), path(bai)
    val clair3_model

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.vcf.gz"), path("${sample_id}_${genotype}.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    echo "=== Variant Calling Clair3 (${clair3_model}) pour ${sample_id} - Génotype ${genotype} ==="
    export CONDA_PREFIX="/opt/conda/envs/VIRiONT2_medaka"
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    # Indexation de la référence canonique
    samtools faidx "${canonical_ref}"

    # Récupération du modèle Clair3
    MODEL_PATH=\$(find \$CONDA_PREFIX -type d -name "${clair3_model}" 2>/dev/null | head -n 1)
    if [ -z "\$MODEL_PATH" ]; then
        MODEL_PATH="${clair3_model}"
    fi
    echo "Utilisation du modèle Clair3 : \$MODEL_PATH"

    # Exécution de Clair3
    run_clair3.sh \
        --bam_fn="${bam}" \
        --ref_fn="${canonical_ref}" \
        --threads=${task.cpus} \
        --platform="ont" \
        --model_path="\$MODEL_PATH" \
        --output=\$PWD/clair3_out \
        --include_all_ctgs

    # Récupération et indexation du VCF
    if [ -f "clair3_out/merge_output.vcf.gz" ]; then
        mv clair3_out/merge_output.vcf.gz "${sample_id}_${genotype}.vcf.gz"
        tabix -f -p vcf "${sample_id}_${genotype}.vcf.gz"
    elif [ -f "clair3_out/phased_merge_output.vcf.gz" ]; then
        mv clair3_out/phased_merge_output.vcf.gz "${sample_id}_${genotype}.vcf.gz"
        tabix -f -p vcf "${sample_id}_${genotype}.vcf.gz"
    else
        echo "❌ ERREUR : Aucun VCF généré par Clair3 dans clair3_out/"
        exit 1
    fi
    """
}