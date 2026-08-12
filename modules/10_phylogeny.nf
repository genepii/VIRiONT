nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 10: MAFFT, IQ-TREE & PLOT PHYLOGÉNÉTIQUE (10_PHYLOGENY)
========================================================================================
*/
process PHYLOGENY {
    publishDir "${params.outdir}/10_PHYLOGENY", mode: 'copy'

    input:
    path genotyped_fastas  // Fichiers FASTA ré-orientés et isolés émis par GENOTYPING
    path ref_seq           // FASTA de référence

    output:
    path "aligned_consensus_and_refs.fasta" , emit: aligned_fasta
    path "*.treefile"                       , optional: true, emit: tree_file
    path "PHYLOGRAM_tree.pdf"                , optional: true, emit: tree_pdf

    script:
    """
    echo "=== 1. Concaténation des consensus orientés (5'->3') et des références ==="
    cat ${genotyped_fastas} ${ref_seq} > unaligned_all.fasta

    echo "=== 2. Alignement MAFFT ==="
    mafft --auto --thread ${task.cpus} unaligned_all.fasta > aligned_consensus_and_refs.fasta

    echo "=== 3. Reconstruction Phylogénétique IQ-TREE ==="
    iqtree -s aligned_consensus_and_refs.fasta -m GTR+G -bb 1000 -nt ${task.cpus} -redo

    echo "=== 4. Génération du rendu graphique PDF via R ==="
    if [ -f "aligned_consensus_and_refs.fasta.treefile" ]; then
        Rscript ${projectDir}/bin/11_plot_tree.R aligned_consensus_and_refs.fasta.treefile PHYLOGRAM_tree.pdf || true
    fi
    """
}