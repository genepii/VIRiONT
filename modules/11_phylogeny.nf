nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 11: ALIGNEMENT MAFFT & ARBRE IQ-TREE (11_PHYLOGENY)
========================================================================================
*/

process PREPARE_TREE_REFS {
    publishDir "${params.outdir}/11_PHYLOGENY", mode: 'copy'

    input:
    path ref_fasta
    path summary_tsv

    output:
    path "tree_references.fasta", emit: filtered_refs

    script:
    """
    echo "=== Sélection des références pour l'arbre ==="
    python3 ${projectDir}/bin/11_prepare_tree_refs.py "${ref_fasta}" "${summary_tsv}" tree_references.fasta ${params.tree_skeleton_n ?: 2}
    """
}

process PHYLOGENY {
    publishDir "${params.outdir}/11_PHYLOGENY", mode: 'copy'

    input:
    path validated_consensus_fasta
    path ref_seq

    output:
    path "aligned_consensus_and_refs.fasta"         , emit: aligned_fasta
    path "aligned_consensus_and_refs.fasta.treefile", emit: tree_file, optional: true

    script:
    """
    echo "=== 1. Concaténation des séquences ==="
    cat "${validated_consensus_fasta}" "${ref_seq}" > unaligned_all.fasta

    NUM_SEQS=\$(grep -c "^>" unaligned_all.fasta || true)
    echo "Nombre total de séquences : \${NUM_SEQS}"

    if [ "\${NUM_SEQS}" -lt 3 ]; then
        echo "Moins de 3 séquences. Arbre non calculable."
        touch aligned_consensus_and_refs.fasta
        exit 0
    fi

    echo "=== 2. Alignement Multiple MAFFT ==="
    mafft --auto --thread ${task.cpus} unaligned_all.fasta > aligned_consensus_and_refs.fasta

    echo "=== 3. Reconstruction de l'Arbre IQ-TREE ==="
    if [ "\${NUM_SEQS}" -ge 4 ]; then
        iqtree -s aligned_consensus_and_refs.fasta -m GTR+G -bb 1000 -nt ${task.cpus} -redo
    else
        iqtree -s aligned_consensus_and_refs.fasta -m GTR+G -nt ${task.cpus} -redo
    fi
    """
}

process PLOT_TREE {
    publishDir "${params.outdir}/11_PHYLOGENY", mode: 'copy'

    input:
    path tree_file

    output:
    path "PHYLOGRAM_tree.pdf", emit: tree_pdf, optional: true

    script:
    """
    echo "=== 4. Rendu graphique PDF via R ==="
    Rscript ${projectDir}/bin/11_plot_tree.R "${tree_file}" PHYLOGRAM_tree.pdf
    """
}