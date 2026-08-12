nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 05: MEDAKA POLISHING (05_CONSENSUS)
========================================================================================
*/
process MEDAKA_CONSENSUS {
    tag { sample_id }
    publishDir { "${params.outdir}/05_CONSENSUS/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq), path(preconsensus_fasta)
    val medaka_model

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: final_consensus

    script:
    """
    echo "=== Polissage Medaka pour ${sample_id} ==="

    medaka_consensus \
        -i "${trimmed_fastq}" \
        -d "${preconsensus_fasta}" \
        -o medaka_out \
        -m "${medaka_model}" \
        -t ${task.cpus}

    # Garantit que si medaka produit plusieurs contigs, chacun reçoit un ID unique (_c1, _c2...)
    awk -v sid="${sample_id}" '
    BEGIN { c = 0; }
    /^>/ {
        c++;
        if (c == 1) {
            print ">" sid;
        } else {
            print ">" sid "_c" c;
        }
        next;
    }
    { print }' medaka_out/consensus.fasta > "${sample_id}.fasta"
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
}