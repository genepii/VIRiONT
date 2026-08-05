/*
========================================================================================
    MODULE 13: ARBRE PHYLOGÉNÉTIQUE IQ-TREE & RENDU VISUEL R
========================================================================================
*/

process PHYLOGENY_IQTREE {
    publishDir "${params.outdir}/09_PHYLOGENY", mode: 'copy'

    input:
    path aligned_fasta

    output:
    path "*.treefile", optional: true, emit: tree_file
    path "PHYLOGRAM_tree.pdf", optional: true, emit: tree_pdf

    script:
    """
    # 1. Reconstruction phylogénétique IQ-TREE
    iqtree -s ${aligned_fasta} -m GTR+G -bb 1000 -nt AUTO

    # 2. Génération du rendu graphique PDF via le script R
    if [ -f "${aligned_fasta}.treefile" ]; then
        13_plot_tree.R ${aligned_fasta}.treefile PHYLOGRAM_tree.pdf
    fi
    """
}