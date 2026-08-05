nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 04: CLUSTERING ET SELECTION AMPLICON_SORTER
========================================================================================
*/

process AMPLICON_SORTER {
    tag { sample_id }
    publishDir path: { "${params.outdir}/04_AMPLICON_SORTER/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq)
    val min_length
    val max_length

    output:
    tuple val(sample_id), path("${sample_id}_*.fasta*"), optional: true, emit: consensus_clusters

    script:
    """
    echo "=== Execution Amplicon_Sorter pour ${sample_id} ==="

    # Restreindre le sur-threading des librairies C et forcer l'affichage direct des logs
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export PYTHONUNBUFFERED=1
    export TMPDIR=\$(pwd)

    # 1. Décompression locale du FastQ
    if [[ "${trimmed_fastq}" == *.gz ]]; then
        gzip -dc ${trimmed_fastq} > input_reads.fastq
    else
        cp ${trimmed_fastq} input_reads.fastq
    fi

    # 2. Lancement en mono-thread (-np 1) pour eviter le blocage multiprocessing sous Singularity
    amplicon_sorter.py \
        -i input_reads.fastq \
        -min ${min_length} \
        -max ${max_length} \
        -maxr 1000 \
        -sfq \
        -c \
        -np 1 \
        -o .

    rm -f input_reads.fastq

    shopt -s nullglob

    # 3. Récupération des consensus depuis les sous-dossiers (y compris .gz)
    find . -mindepth 2 -type f \\( -name "*.fasta*" -o -name "*.fa*" -o -name "*.group*" -o -name "*.sorted*" \\) -exec mv {} . \\;

    # 4. Conversion explicite en FASTA si des fichiers .group ou .sorted existent
    for g in *.group *.sorted; do
        if [ -f "\$g" ]; then
            out_fa="\${g%.*}.fasta"
            awk 'NR%4==1{sub(/^@/," >");print} NR%4==2{print}' "\$g" > "\$out_fa"
        fi
    done

    # 5. Normalisation des noms de fichiers pour Nextflow
    for f in *.fasta *.fasta.gz *.fa *.fa.gz; do
        if [ -f "\$f" ]; then
            if [[ "\$f" != "${sample_id}_"* ]]; then
                mv "\$f" "${sample_id}_\$f"
            fi
        fi
    done

    shopt -u nullglob
    """
}