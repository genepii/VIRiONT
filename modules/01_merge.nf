nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 01: MERGE FASTQ
========================================================================================
*/

process MERGE_FASTQ {
    tag "$sample_id"
    publishDir path: { "${params.outdir}/01_MERGED" }, mode: 'copy'

    input:
    tuple val(sample_id), path(input_dir)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz"), emit: merged_fastq, optional: true

    script:
    """
    echo "=== Traitement de ${sample_id} dans : ${input_dir} ==="

    if [ ! -d "${input_dir}" ]; then
        echo "ℹ️ DOSSIER ABSENT : ${input_dir} n'existe pas. Échantillon/témoin ignoré."
        exit 0
    fi

    gz_count=\$(find -L "${input_dir}" -maxdepth 1 -type f -name "*.gz" 2>/dev/null | wc -l)
    fastq_count=\$(find -L "${input_dir}" -maxdepth 1 -type f -name "*.fastq" 2>/dev/null | wc -l)

    echo "Trouve : \$gz_count fichier(s) .gz et \$fastq_count fichier(s) .fastq"

    if [ "\$gz_count" -gt 0 ] && [ "\$fastq_count" -eq 0 ]; then
        find -L "${input_dir}" -maxdepth 1 -type f -name "*.gz" -exec cat {} + > ${sample_id}_merged.fastq.gz
    elif [ "\$fastq_count" -gt 0 ] && [ "\$gz_count" -eq 0 ]; then
        find -L "${input_dir}" -maxdepth 1 -type f -name "*.fastq" -exec cat {} + | gzip -c > ${sample_id}_merged.fastq.gz
    elif [ "\$gz_count" -gt 0 ] || [ "\$fastq_count" -gt 0 ]; then
        {
            find -L "${input_dir}" -maxdepth 1 -type f -name "*.fastq" -exec cat {} +
            find -L "${input_dir}" -maxdepth 1 -type f -name "*.gz" -exec gzip -dc {} +
        } | gzip -c > ${sample_id}_merged.fastq.gz
    else
        echo "ℹ️ TÉMOIN NÉGATIF / DOSSIER VIDE : Aucun fichier FASTQ dans ${input_dir}. Poursuite du pipeline."
        exit 0
    fi

    file_size=\$(stat -c%s "${sample_id}_merged.fastq.gz" 2>/dev/null || stat -f%z "${sample_id}_merged.fastq.gz")
    if [ "\$file_size" -lt 1000 ]; then
        echo "ℹ️ VOLUMÉTRIE INSUFFISANTE : Le fichier fait \${file_size} octets. Échantillon/témoin ignoré."
        rm -f ${sample_id}_merged.fastq.gz
        exit 0
    fi
    """
}