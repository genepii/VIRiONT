nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 06b: ALIGNEMENT PAR GÉNOTYPE VALIDÉ (06b_GENOTYPE_BAM)
    Réaligne TOUS les reads trimmed du barcode contre CHAQUE consensus de génotype
    validé individuellement (post-agrégation GENOTYPING). Nécessaire en cas de
    co-infection pour garantir un BAM/VCF/couverture strictement propre à chaque
    génotype, sans mélange avec les autres clusters bruts ou les autres génotypes.
========================================================================================
*/
process ALIGN_GENOTYPE_BAM {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/06b_GENOTYPE_BAM/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(trimmed_fastq), path(genotype_fasta)

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), emit: bam_bai
    tuple val(sample_id), val(genotype), path(genotype_fasta), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), emit: for_vcf

    script:
    """
    echo "=== Alignement minimap2 spécifique génotype ${genotype} pour ${sample_id} ==="
    samtools faidx ${genotype_fasta}
    minimap2 -ax map-ont -t ${task.cpus} ${genotype_fasta} ${trimmed_fastq} | \\
        samtools view -bS -F 4 - | \\
        samtools sort -o ${sample_id}_${genotype}.sorted.bam -
    samtools index ${sample_id}_${genotype}.sorted.bam
    """
}