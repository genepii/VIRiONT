nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 10: ALIGNEMENT MAFFT & ARBRE IQ-TREE (10_PHYLOGENY)
========================================================================================
*/
process PREPARE_TREE_REFS {
    publishDir "${params.outdir}/10_PHYLOGENY", mode: 'copy'
    input:
    path ref_fasta      // Base de référence complète (tous génotypes/sous-génotypes)
    path summary_tsv    // Génotypes validés du run
    output:
    path "tree_references.fasta", emit: filtered_refs
    script:
    """
    echo "=== Sélection des références pour l'arbre (squelette=${params.tree_skeleton_n}) ==="
    python3 ${projectDir}/bin/10_prepare_tree_refs.py "${ref_fasta}" "${summary_tsv}" tree_references.fasta ${params.tree_skeleton_n}
    """
}

process PHYLOGENY {
    publishDir "${params.outdir}/10_PHYLOGENY", mode: 'copy'
    input:
    path validated_consensus_fasta  // Fichier unifié du fil conducteur
    path ref_seq                    // Références filtrées (squelette + génotypes présents)
    output:
    path "aligned_consensus_and_refs.fasta" , emit: aligned_fasta
    path "*.treefile"                       , optional: true, emit: tree_file
    path "PHYLOGRAM_tree.pdf"               , optional: true, emit: tree_pdf
    script:
    """
    echo "=== 1. Concaténation des séquences validées du run et des références ==="
    cat "${validated_consensus_fasta}" "${ref_seq}" > unaligned_all.fasta
    echo "=== 2. Alignement Multiple MAFFT ==="
    mafft --auto --thread ${task.cpus} unaligned_all.fasta > aligned_consensus_and_refs.fasta
    echo "=== 3. Reconstruction de l'Arbre IQ-TREE ==="
    iqtree -s aligned_consensus_and_refs.fasta -m GTR+G -bb 1000 -nt ${task.cpus} -redo
    echo "=== 4. Rendu graphique PDF via R ==="
    if [ -f "aligned_consensus_and_refs.fasta.treefile" ]; then
        Rscript ${projectDir}/bin/10_plot_tree.R aligned_consensus_and_refs.fasta.treefile PHYLOGRAM_tree.pdf || true
    fi
    """
}