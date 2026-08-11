nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 07: VARIANT CALLING (07_VCF - CLAIR3)
========================================================================================
*/
process CALL_VCF {
    tag "$sample_id"
    publishDir { "${params.outdir}/07_VCF/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(final_consensus), path(bam), path(bai)
    val clair3_model

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    echo "=== Variant Calling Clair3 (${clair3_model}) pour ${sample_id} ==="

    # 1. Activation de l'environnement Conda pour Clair3
    export CONDA_PREFIX="/opt/conda/envs/VIRiONT2_medaka"
    export PATH="\$CONDA_PREFIX/bin:\$PATH"

    # 2. Dé-symlinker + ASSAINIR les noms de contigs du consensus.
    #    amplicon_sorter produit des headers du type ">consensus_..._0(8845)" —
    #    les parenthèses cassent la commande shell interne de Clair3.
    #    On remplace tout caractère non alphanumérique/./_/- par "_", en préservant le ">".
    cp -L ${final_consensus} raw_ref.fasta
    sed -E '/^>/ s/[^A-Za-z0-9._>-]/_/g' raw_ref.fasta > ref.fasta
    samtools faidx ref.fasta

    # 3. Ré-header le BAM pour que ses noms de contigs (@SQ) correspondent
    #    exactement à ceux assainis dans ref.fasta (même transformation).
    samtools view -H ${bam} \\
        | awk 'BEGIN{OFS="\\t"} /^@SQ/{ for(i=1;i<=NF;i++){ if(\$i ~ /^SN:/){ split(\$i,a,"SN:"); name=a[2]; gsub(/[^A-Za-z0-9._-]/,"_",name); \$i="SN:"name } } } {print}' > new_header.sam
    samtools reheader new_header.sam ${bam} > ${sample_id}.reheadered.bam
    samtools index ${sample_id}.reheadered.bam

    # 4. Localisation automatique du chemin du modèle Clair3
    MODEL_PATH=\$(find \$CONDA_PREFIX -type d -name "${clair3_model}" 2>/dev/null | head -n 1)

    if [ -z "\$MODEL_PATH" ]; then
        MODEL_PATH="${clair3_model}"
    fi

    echo "Utilisation du modèle Clair3 : \$MODEL_PATH"

    # 5. Exécution de Clair3
    run_clair3.sh \\
        --bam_fn=${sample_id}.reheadered.bam \\
        --ref_fn=ref.fasta \\
        --threads=${task.cpus} \\
        --platform="ont" \\
        --model_path="\$MODEL_PATH" \\
        --output=\$PWD/clair3_out \\
        --include_all_ctgs

    # 6. Récupération et indexation du VCF final
    if [ -f "clair3_out/merge_output.vcf.gz" ]; then
        mv clair3_out/merge_output.vcf.gz ${sample_id}.vcf.gz
        tabix -f -p vcf ${sample_id}.vcf.gz
    elif [ -f "clair3_out/phased_merge_output.vcf.gz" ]; then
        mv clair3_out/phased_merge_output.vcf.gz ${sample_id}.vcf.gz
        tabix -f -p vcf ${sample_id}.vcf.gz
    else
        echo "❌ ERREUR : Aucun VCF généré par Clair3 dans clair3_out/"
        exit 1
    fi
    """
}