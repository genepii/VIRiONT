nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 05: MEDAKA POLISHING PAR GÉNOTYPE (05_CONSENSUS)
========================================================================================
*/
process MEDAKA_CONSENSUS {
    tag { "${sample_id}_${genotype}" }
    publishDir { "${params.outdir}/05_CONSENSUS/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(geno_fastq), path(ref_fasta)
    val medaka_model

    output:
    tuple val(sample_id), val(genotype), path("*_Geno_*.fasta"), emit: final_consensus

    script:
    """
    export KMP_DUPLICATE_LIB_OK=TRUE
    export OMP_NUM_THREADS=1

echo "=== Polissage Medaka pour ${sample_id} - Génotype ${genotype} ==="

    medaka_consensus \\
        -i "${geno_fastq}" \\
        -d "${ref_fasta}" \\
        -o medaka_out \\
        -m "${medaka_model}" \\
        -t ${task.cpus}

    awk -v header="${sample_id}_Geno_${genotype}" '
        /^>/ { print ">" header; next }
        { print }
    ' medaka_out/consensus.fasta > "${sample_id}_Geno_${genotype}.fasta"
    """

    stub:
    """
    echo -e ">${sample_id}_Geno_${genotype}\nACGT" > "${sample_id}_Geno_${genotype}.fasta"
    """
}

process COLLECT_CONSENSUS {
    publishDir "${params.outdir}/05_CONSENSUS", mode: 'copy'

    input:
    path consensus_files

    output:
    path "all_cons.fasta", emit: all_consensus

    script:
    """
    cat ${consensus_files} > all_cons.fasta
    """

    stub:
    """
    touch all_cons.fasta
    """
}