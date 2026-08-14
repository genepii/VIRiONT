nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 04: AMPLICON_SORTER (04_PRECONSENSUS)
========================================================================================
*/
process AMPLICON_SORTER {
    tag "$sample_id"
    container null

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

    # Comptage du nombre réel de reads dans le fichier trimmed (4 lignes par read en FASTQ)
    N_READS=\$(( \$(zcat -f ${trimmed_fastq} | wc -l) / 4 ))
    echo "Nombre de reads détectés pour ${sample_id} : \$N_READS"

    # Logique adaptative :
    # - Si le nombre total de reads est déjà <= max_reads, on les prend TOUS (-ar), pas de sous-échantillonnage.
    # - Sinon, on sous-échantillonne à max_reads, mais de façon ALEATOIRE (-ra) pour ne pas biaiser
    #   la sélection vers l'ordre du fichier (biais temporel/run potentiel).
    if [ "\$N_READS" -le "${params.max_reads}" ]; then
        READ_OPTS="-ar -maxr ${params.max_reads}"
        echo "--> Tous les reads seront utilisés (N_READS <= max_reads)."
    else
        READ_OPTS="-ra -maxr ${params.max_reads}"
        echo "--> Sous-échantillonnage aléatoire à ${params.max_reads} reads (N_READS > max_reads)."
    fi

    /usr/bin/singularity exec --no-home -B ${projectDir}:${projectDir} -B \$PWD:/tmp ${params.sif_main} bash -c "
        export MPLBACKEND=Agg
        export QT_QPA_PLATFORM=offscreen
        export TMPDIR=/tmp
        amplicon_sorter.py \\
            -i ${trimmed_fastq} \\
            -o ${sample_id}_results \\
            -min ${min_length} \\
            -max ${max_length} \\
            \$READ_OPTS \\
            -np ${task.cpus} < /dev/null
    "

    if [ ! -f "${sample_id}_results/consensusfile.fasta" ]; then
        echo "⚠️ WARNING : amplicon_sorter n'a produit aucun pré-consensus pour ${sample_id}."
    fi
    """
}