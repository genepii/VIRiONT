nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 07: VARIANT CALLING (07_VCF - CLAIR3, PAR GÉNOTYPE)
========================================================================================
*/
process CALL_VCF {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/07_VCF/${sample_id}/${genotype}" }, mode: 'copy'
    input:
    tuple val(sample_id), val(genotype), path(final_consensus), path(bam), path(bai)
    val clair3_model
    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.vcf.gz"), path("${sample_id}_${genotype}.vcf.gz.tbi"), emit: vcf_tbi
    script:
    """
    echo "=== Variant Calling Clair3 (${clair3_model}) pour ${sample_id} - Génotype ${genotype} ==="
    export CONDA_PREFIX="/opt/conda/envs/VIRiONT2_medaka"
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    # Le header du consensus par génotype est déjà propre ("{sample_id}_Geno_{genotype}",
    # sans parenthèses), mais on garde ce nettoyage défensif au cas où.
    cp -L ${final_consensus} raw_ref.fasta
    sed -E '/^>/ s/[^A-Za-z0-9._>-]/_/g' raw_ref.fasta > ref.fasta
    samtools faidx ref.fasta

    samtools view -H ${bam} \\
        | awk 'BEGIN{OFS="\\t"} /^@SQ/{ for(i=1;i<=NF;i++){ if(\$i ~ /^SN:/){ split(\$i,a,"SN:"); name=a[2]; gsub(/[^A-Za-z0-9._-]/,"_",name); \$i="SN:"name } } } {print}' > new_header.sam
    samtools reheader new_header.sam ${bam} > ${sample_id}_${genotype}.reheadered.bam
    samtools index ${sample_id}_${genotype}.reheadered.bam

    MODEL_PATH=\$(find \$CONDA_PREFIX -type d -name "${clair3_model}" 2>/dev/null | head -n 1)
    if [ -z "\$MODEL_PATH" ]; then
        MODEL_PATH="${clair3_model}"
    fi
    echo "Utilisation du modèle Clair3 : \$MODEL_PATH"

    run_clair3.sh \\
        --bam_fn=${sample_id}_${genotype}.reheadered.bam \\
        --ref_fn=ref.fasta \\
        --threads=${task.cpus} \\
        --platform="ont" \\
        --model_path="\$MODEL_PATH" \\
        --output=\$PWD/clair3_out \\
        --include_all_ctgs

    if [ -f "clair3_out/merge_output.vcf.gz" ]; then
        mv clair3_out/merge_output.vcf.gz ${sample_id}_${genotype}.vcf.gz
        tabix -f -p vcf ${sample_id}_${genotype}.vcf.gz
    elif [ -f "clair3_out/phased_merge_output.vcf.gz" ]; then
        mv clair3_out/phased_merge_output.vcf.gz ${sample_id}_${genotype}.vcf.gz
        tabix -f -p vcf ${sample_id}_${genotype}.vcf.gz
    else
        echo "❌ ERREUR : Aucun VCF généré par Clair3 dans clair3_out/"
        exit 1
    fi
    """
}