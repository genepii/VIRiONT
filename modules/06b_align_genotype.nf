nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 06b: ALIGNEMENT CONTRE RÉFÉRENCE CANONIQUE DU GÉNOTYPE (06b_GENOTYPE_BAM)
    Aligne les reads contre la séquence canonique du génotype (HBV_genotype_wPrimer.fasta)
    pour que les coordonnées du VCF matchent exactement les tables mutation_GT*.csv.
========================================================================================
*/
process ALIGN_GENOTYPE_BAM {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/06b_GENOTYPE_BAM/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(trimmed_fastq), path(consensus_fasta)
    path hbv_primer_ref

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), emit: bam_bai
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_canonical_ref.fasta"), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), emit: for_vcf

    script:
    """
    # Extraction du code de génotype majeur (ex. E -> GTE, D3 -> GTD, A1 -> GTA)
    clean_gt=\$(echo "${genotype}" | sed -E 's/.*([A-J]).*/GT\\1/')
    echo "=== Alignement minimap2 contre référence canonique \${clean_gt} pour ${sample_id} ==="

    # Extraction du contig canonique correspondant depuis HBV_genotype_wPrimer.fasta
    samtools faidx "${hbv_primer_ref}" "\${clean_gt}" > "${sample_id}_${genotype}_canonical_ref.fasta" || \
    samtools faidx "${hbv_primer_ref}" "GTD" > "${sample_id}_${genotype}_canonical_ref.fasta"

    samtools faidx "${sample_id}_${genotype}_canonical_ref.fasta"

    # Alignement des reads trimmed
    minimap2 -ax map-ont -t ${task.cpus} "${sample_id}_${genotype}_canonical_ref.fasta" "${trimmed_fastq}" | \
        samtools view -bS -F 4 - | \
        samtools sort -o "${sample_id}_${genotype}.sorted.bam" -
    samtools index "${sample_id}_${genotype}.sorted.bam"
    """
}