nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 09: VARIANT CALLING (SORTIE TEXTE CLAIR .VCF)
    - CALL_VCF           : Clair3 décompressé (alimente le QC)
    - RAW_PILEUP_COUNTS  : comptage brut des substitutions (bcftools .vcf clair)
    - RAW_INDEL_COUNTS   : comptage brut des indels (court-circuité .vcf clair)
========================================================================================
*/

process CALL_VCF {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/09_VCF/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(canonical_ref), path(bam), path(bai)
    val clair3_model

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.vcf"), emit: vcf

    script:
    """
    echo "=== Variant Calling Clair3 (${clair3_model}) pour ${sample_id} Génotype ${genotype} ==="
    export CONDA_PREFIX="/opt/conda/envs/clair3"
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    samtools faidx "${canonical_ref}"

    if [ -d "\$CONDA_PREFIX/bin/models/${clair3_model}" ]; then
        MODEL_PATH="\$CONDA_PREFIX/bin/models/${clair3_model}"
    else
        MODEL_PATH="${clair3_model}"
    fi
    echo "Utilisation du modèle Clair3 : \$MODEL_PATH"

    run_clair3.sh \\
        --bam_fn="${bam}" \\
        --ref_fn="${canonical_ref}" \\
        --threads=${task.cpus} \\
        --platform="ont" \\
        --model_path="\$MODEL_PATH" \\
        --output=\$PWD/clair3_out \\
        --include_all_ctgs \\
        --snp_min_af=0.005 \\
        --indel_min_af=0.01 \\
        --var_pct_full=0.005 \\
        --ref_pct_full=0.005 \\
        --min_coverage=10 \\
        --min_mq=5 \\
        --qual=0 \\
        --pileup_only \\
        --print_ref_calls \\
        --enable_variant_calling_at_sequence_head_and_tail

    if [ -f "clair3_out/pileup.vcf.gz" ]; then
        gzip -dc clair3_out/pileup.vcf.gz > "${sample_id}_${genotype}.vcf"
    elif [ -f "clair3_out/merge_output.vcf.gz" ]; then
        gzip -dc clair3_out/merge_output.vcf.gz > "${sample_id}_${genotype}.vcf"
    else
        echo "ERREUR : Aucun VCF généré par Clair3 dans clair3_out/"
        exit 1
    fi
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' > "${sample_id}_${genotype}.vcf"
    """
}

process RAW_PILEUP_COUNTS {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/09_VCF/${sample_id}/${genotype}" }, mode: 'copy'

    container "${projectDir}/singularity/clair3.sif"

    input:
    tuple val(sample_id), val(genotype), path(canonical_ref), path(bam), path(bai)

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.raw_counts.vcf"), emit: raw_counts

    script:
    """
    [ ! -f "${canonical_ref}.fai" ] && samtools faidx "${canonical_ref}"

    # -Ov génère directement un fichier texte VCF décompressé
    bcftools mpileup \\
        -f "${canonical_ref}" \\
        -a FORMAT/AD,FORMAT/DP \\
        -d 1000000 -L 1000000 \\
        -Q 0 -q 0 -B -A -I \\
        --no-version \\
        -Ov -o "${sample_id}_${genotype}.raw_counts.vcf" \\
        "${bam}"
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' > "${sample_id}_${genotype}.raw_counts.vcf"
    """
}

process RAW_INDEL_COUNTS {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/09_VCF/${sample_id}/${genotype}" }, mode: 'copy'

    container "${projectDir}/singularity/clair3.sif"

    input:
    tuple val(sample_id), val(genotype), path(canonical_ref), path(bam), path(bai)

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.raw_indels.vcf"), emit: raw_indels

    script:
    """
    echo "=== [BYPASS] RAW_INDEL_COUNTS temporairement ignoré pour ${sample_id} ==="

    # --- CODE ORIGINAL CONSERVÉ (-Ov pour texte brut) ---
    # [ ! -f "${canonical_ref}.fai" ] && samtools faidx "${canonical_ref}"
    # bcftools mpileup \\
    #     -f "${canonical_ref}" \\
    #     -a FORMAT/AD,FORMAT/DP \\
    #     -d 1000000 -L 1000000 --max-idepth 1000000 \\
    #     -m 1 -F 0.0001 \\
    #     -Q 0 -q 0 -B -A \\
    #     --no-version -Ou "${bam}" \\
    #   | bcftools view -v indels -Ov -o "${sample_id}_${genotype}.raw_indels.vcf"
    # ----------------------------------------------------

    echo '##fileformat=VCFv4.2' > "${sample_id}_${genotype}.raw_indels.vcf"
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' > "${sample_id}_${genotype}.raw_indels.vcf"
    """
}