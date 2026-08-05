/*
========================================================================================
    MODULE 09: CONCATÉNATION DES SÉQUENCES CONSENSUS POLIES
========================================================================================
*/

process CONCAT_CONSENSUS {
    publishDir "${params.outdir}/09_PHYLOGENY", mode: 'copy'

    input:
    path polished_fastas

    output:
    path "all_consensus.fasta", emit: all_consensus

    script:
    """
    cat ${polished_fastas} > all_consensus.fasta
    """
}