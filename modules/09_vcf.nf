nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 09: VARIANT CALLING CLAIR3 (09_VCF - COORDONNÉES CANONIQUES)
========================================================================================
*/
process CALL_VCF {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/09_VCF/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(canonical_ref), path(bam), path(bai)
    val clair3_model

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.vcf.gz"), path("${sample_id}_${genotype}.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    echo "=== Variant Calling Clair3 (${clair3_model}) pour ${sample_id} - Génotype ${genotype} ==="
    export CONDA_PREFIX="/opt/conda/envs/clair3"
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    # Indexation de la référence canonique
    samtools faidx "${canonical_ref}"

    # Récupération du modèle Clair3 embarqué
    if [ -d "\$CONDA_PREFIX/bin/models/${clair3_model}" ]; then
        MODEL_PATH="\$CONDA_PREFIX/bin/models/${clair3_model}"
    else
        MODEL_PATH="${clair3_model}"
    fi
    echo "Utilisation du modèle Clair3 : \$MODEL_PATH"

    # Exécution de Clair3 calibrée pour quasi-espèces virales et signaux sous-clonaux
    run_clair3.sh \
        --bam_fn="${bam}" \
        --ref_fn="${canonical_ref}" \
        --threads=${task.cpus} \
        --platform="ont" \
        --model_path="\$MODEL_PATH" \
        --output=\$PWD/clair3_out \
        --include_all_ctgs \
        --snp_min_af=0.001 \
        --indel_min_af=0.001 \
        --var_pct_full=0.01 \
        --ref_pct_full=0.01 \
        --qual=0 \
        --pileup_only \
        --haploid_sensitive \
        --print_ref_calls \
        --enable_variant_calling_at_sequence_head_and_tail

    # Récupération et indexation du VCF
    if [ -f "clair3_out/pileup.vcf.gz" ]; then
        mv clair3_out/pileup.vcf.gz "${sample_id}_${genotype}.vcf.gz"
        tabix -f -p vcf "${sample_id}_${genotype}.vcf.gz"
    elif [ -f "clair3_out/merge_output.vcf.gz" ]; then
        mv clair3_out/merge_output.vcf.gz "${sample_id}_${genotype}.vcf.gz"
        tabix -f -p vcf "${sample_id}_${genotype}.vcf.gz"
    else
        echo "ERREUR : Aucun VCF généré par Clair3 dans clair3_out/"
        exit 1
    fi
    """

    stub:
    """
    echo "=== [STUB] Variant Calling Clair3 pour ${sample_id} - Génotype ${genotype} ==="
    echo '##fileformat=VCFv4.2' | bgzip -c > "${sample_id}_${genotype}.vcf.gz"
    touch "${sample_id}_${genotype}.vcf.gz.tbi"
    """
}
