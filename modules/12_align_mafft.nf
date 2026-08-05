/*
========================================================================================
    MODULE 12: ALIGNEMENT MAFFT & SÉLECTION DES RÉFÉRENCES
========================================================================================
*/

process ALIGN_MAFFT {
    publishDir "${params.outdir}/09_PHYLOGENY", mode: 'copy'

    input:
    path consensus_fasta
    path ref_seq
    path validated_tsvs

    output:
    path "aligned_consensus_and_refs.fasta", emit: aligned_fasta

    script:
    """
    # 1. Concaténation des séquences consensus du run avec les références
    cat ${consensus_fasta} ${ref_seq} > unaligned_all.fasta

    # 2. Alignement multiple avec MAFFT
    mafft --auto --thread ${task.cpus} unaligned_all.fasta > aligned_consensus_and_refs.fasta
    """
}