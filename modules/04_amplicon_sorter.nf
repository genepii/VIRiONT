nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 04: CLUSTERING ET SELECTION (AMPLICON_SORTER)
========================================================================================
*/

process AMPLICON_SORTER {
    tag "$sample_id"
    publishDir path: { "${params.outdir}/04_AMPLICON_SORTER/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq)
    val min_length
    val max_length

    output:
    tuple val(sample_id), path("${sample_id}_*.fasta"), optional: true, emit: consensus_clusters

    script:
    def threads = task.cpus ?: 4
    """
    echo "=== Execution Amplicon_Sorter pour ${sample_id} [${min_length}-${max_length} bp] avec ${threads} threads ==="

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export PYTHONUNBUFFERED=1

    # 1. Décompression du FASTQ
    if [[ "${trimmed_fastq}" == *.gz ]]; then
        gunzip -c ${trimmed_fastq} > raw_input.fastq
    else
        cp ${trimmed_fastq} raw_input.fastq
    fi

    # 2. Comptage des reads disponibles
    NUM_READS=\$(awk 'NR%4==1' raw_input.fastq | wc -l)
    echo "Reads disponibles pour ${sample_id} : \$NUM_READS"

    if [ "\$NUM_READS" -ge 10 ]; then
        # Subsampling de sécurité à 5 000 reads max
        if [ "\$NUM_READS" -gt 5000 ]; then
            echo "⚡ Sous-échantillonnage à 5000 reads pour accélérer Amplicon_Sorter..."
            if command -v seqkit &> /dev/null; then
                seqkit sample -n 5000 raw_input.fastq -o input_reads.fastq
            else
                paste - - - - < raw_input.fastq | shuf -n 5000 | tr '\\t' '\\n' > input_reads.fastq
            fi
            rm -f raw_input.fastq
        else
            mv raw_input.fastq input_reads.fastq
        fi

        # 3. Lancement d'Amplicon_Sorter avec redirection explicite pour débloquer Singularity
        python3 \$(which amplicon_sorter.py) \
            -i input_reads.fastq \
            -o "${sample_id}" \
            -min ${min_length} \
            -max ${max_length} \
            -maxr 5000 \
            -ar \
            -ssg 85 \
            -ss 85 \
            -sc 92 \
            -np ${threads} > amplicon_sorter.log 2>&1 || true

        # 4. Conversion des fichiers .sorted / .group générés en FASTA si nécessaire
        for g in *.group *.sorted; do
            if [ -f "\$g" ]; then
                out_fa="\${g%.*}.fasta"
                awk 'NR%4==1{sub(/^@/," >");print} NR%4==2{print}' "\$g" > "\$out_fa"
            fi
        done

        # 5. Nettoyage
        rm -f *unique*.fasta input_reads.fastq
    else
        echo "⚠️ Trop peu de reads pour effectuer le clustering (\$NUM_READS reads)."
        rm -f raw_input.fastq
    fi
    """
}