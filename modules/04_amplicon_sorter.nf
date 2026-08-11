nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 04: AMPLICON_SORTER (04_PRECONSENSUS)
========================================================================================
*/
process AMPLICON_SORTER {
    tag "$sample_id"
    container null
    
    // Utilisation des accolades {} pour évaluer $sample_id dynamiquement après le bloc input
    publishDir { "${params.outdir}/04_PRECONSENSUS/${sample_id}" }, mode: 'copy', saveAs: { filename -> filename.replace("${sample_id}_results/", "") }

    input:
    tuple val(sample_id), path(trimmed_fastq)
    val min_length
    val max_length

    output:
    tuple val(sample_id), path("${sample_id}_results/consensusfile.fasta"),                    emit: preconsensus_fasta, optional: true
    tuple val(sample_id), path("${sample_id}_results/results.csv"),                          emit: results_csv,         optional: true
    tuple val(sample_id), path("${sample_id}_results/results.txt"),                          emit: results_txt,         optional: true
    tuple val(sample_id), path("${sample_id}_results/${sample_id}_trimmed"),                 emit: sorter_details,      optional: true

    script:
    """
    echo "=== Amplicon_sorter (Pre-consensus) pour ${sample_id} [${min_length} bp - ${max_length} bp] ==="

    /usr/bin/singularity exec --no-home -B ${projectDir}:${projectDir} -B \$PWD:/tmp ${params.sif_main} bash -c '
        export MPLBACKEND=Agg
        export QT_QPA_PLATFORM=offscreen
        export TMPDIR=/tmp
        amplicon_sorter.py \
            -i ${trimmed_fastq} \
            -o ${sample_id}_results \
            -min ${min_length} \
            -max ${max_length} \
            -maxr ${params.max_reads} \
            -np ${task.cpus} < /dev/null
    '

    if [ ! -f "${sample_id}_results/consensusfile.fasta" ]; then
        echo "⚠️ WARNING : amplicon_sorter n'a produit aucun pré-consensus pour ${sample_id}."
    fi
    """
}