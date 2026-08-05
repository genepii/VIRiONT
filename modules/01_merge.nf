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
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz"), emit: merged_fastq

    script:
    """
    echo "=== Traitement de ${sample_id} dans : ${input_dir} ==="

    gz_count=\$(find -L "${input_dir}" -maxdepth 1 -type f -name "*.gz" | wc -l)
    fastq_count=\$(find -L "${input_dir}" -maxdepth 1 -type f -name "*.fastq" | wc -l)

    echo "Trouve : \$gz_count fichier(s) .gz et \$fastq_count fichier(s) .fastq"

    if [ "\$gz_count" -gt 0 ] && [ "\$fastq_count" -eq 0 ]; then
        find -L "${input_dir}" -maxdepth 1 -type f -name "*.gz" -exec cat {} + > ${sample_id}_merged.fastq.gz
    elif [ "\$fastq_count" -gt 0 ] && [ "\$gz_count" -eq 0 ]; then
        find -L "${input_dir}" -maxdepth 1 -type f -name "*.fastq" -exec cat {} + | gzip -c > ${sample_id}_merged.fastq.gz
    elif [ "\$gz_count" -gt 0 ] || [ "\$fastq_count" -gt 0 ]; then
        {
            find -L "${input_dir}" -maxdepth 1 -type f -name "*.fastq" -exec cat {} +
            find -L "${input_dir}" -maxdepth 1 -type f -name "*.gz" -exec zcat {} +
        } | gzip -c > ${sample_id}_merged.fastq.gz
    else
        echo "❌ ERREUR : Aucun fichier FASTQ ou FASTQ.GZ trouve dans ${input_dir}"
        exit 1
    fi

    # Seuil de sécurité abaissé à 10 Ko (filtre uniquement les fichiers header de 20 octets)
    file_size=\$(stat -c%s "${sample_id}_merged.fastq.gz")
    if [ "\$file_size" -lt 10000 ]; then
        echo "❌ ERREUR CRITIQUE : Le fichier généré est vide ou corrompu (\${file_size} octets)"
        exit 1
    fi
    """
}