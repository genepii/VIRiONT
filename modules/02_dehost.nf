nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 02: DEHOSTING (HOSTILE)
========================================================================================
*/

process DEHOST_HOSTILE {
    tag "$sample_id"
    publishDir "${params.outdir}/02_DEHOSTING", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_dehosted.fastq.gz"), emit: dehosted_fastq

    script:
    """
    echo "=== 02. Nettoyage du génome hôte (Hostile) pour ${sample_id} ==="

    INDEX_NAME="${params.hostile_index}"
    CACHE_DIR="${params.hostile_cache}"

    export HOSTILE_CACHE_DIR="\${CACHE_DIR}"
    mkdir -p "\${CACHE_DIR}"

    # Export des variables proxy CHU Lyon
    export http_proxy="http://di3224su:8888/"
    export https_proxy="http://di3224su:8888/"
    export HTTP_PROXY="http://di3224su:8888/"
    export HTTPS_PROXY="http://di3224su:8888/"

    # 1. Téléchargement automatique de l'index si absent du cache local
    if [ ! -d "\${CACHE_DIR}/\${INDEX_NAME}" ] && [ ! -f "\${CACHE_DIR}/\${INDEX_NAME}.mmi" ]; then
        echo " Index \${INDEX_NAME} absent du cache local. Téléchargement via proxy..."
        hostile index fetch --name "\${INDEX_NAME}" --aligner minimap2 || true
    fi

    # 2. Activation du mode hors-ligne si l'index est présent
    AIRPLANE=""
    if [ -d "\${CACHE_DIR}/\${INDEX_NAME}" ] || [ -f "\${CACHE_DIR}/\${INDEX_NAME}.mmi" ]; then
        AIRPLANE="--airplane"
    fi

    # 3. Nettoyage des reads
    hostile clean \
        --fastq1 ${merged_fastq} \
        --index "\${INDEX_NAME}" \
        --aligner minimap2 \
        --threads ${task.cpus} \
        \$AIRPLANE \
        -o .

    # 4. Normalisation de la sortie
    if [ -f "${sample_id}_merged.clean.fastq.gz" ]; then
        mv "${sample_id}_merged.clean.fastq.gz" "${sample_id}_dehosted.fastq.gz"
    else
        alt_file=\$(find . -maxdepth 1 -name "*.clean.fastq.gz" | head -n 1)
        if [ -n "\$alt_file" ]; then
            mv "\$alt_file" "${sample_id}_dehosted.fastq.gz"
        else
            echo "❌ Fichier Hostile non généré pour ${sample_id}"
            exit 1
        fi
    fi
    """
}