nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 05: ALIGNEMENT & CONSENSUS IUPAC PAR SOUS-TYPE
========================================================================================
*/

// module consensus


process ALIGN_SUBTYPE {
    tag { "${sample_id}_${genotype}" }

    input:
    tuple val(sample_id), val(genotype), path(geno_fastq), path(ref_fasta)

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), path(ref_fasta), emit: bam_ref

    script:
    """
    [ ! -f "${ref_fasta}.fai" ] && samtools faidx "${ref_fasta}"

    # Alignement avec tolérance aux délétions (-r 2000 -g 2000)
    minimap2 -ax map-ont -r 2000 -g 2000 --secondary=no -t ${task.cpus} "${ref_fasta}" "${geno_fastq}" | \
        samtools view -b -F 2308 | \
        samtools sort -@ ${task.cpus} -o "${sample_id}_${genotype}.sorted.bam"

    samtools index "${sample_id}_${genotype}.sorted.bam"
    """

    stub:
    """
    touch "${sample_id}_${genotype}.sorted.bam" "${sample_id}_${genotype}.sorted.bam.bai"
    """
}

process GENERATE_CONSENSUS {
    tag { "${sample_id}_${genotype}" }
    publishDir { "${params.outdir}/05_CONSENSUS/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_Geno_${genotype}.fasta"), emit: final_consensus

    script:
    def min_depth = params.consensus_min_depth ?: 10
    def min_iupac = params.iupac_min_freq != null ? params.iupac_min_freq : 0.08

    """
    echo "=== Consensus IUPAC multithreadé pour ${sample_id} (${genotype}) ==="

    [ ! -f "${ref_fasta}.fai" ] && samtools faidx "${ref_fasta}"

    CHR_NAME=\$(head -n 1 "${ref_fasta}.fai" | cut -f1)
    CHR_LEN=\$(head -n 1 "${ref_fasta}.fai" | cut -f2)

    STEP=\$(( (CHR_LEN + 3) / 4 ))

    for i in 0 1 2 3; do
        START=\$(( i * STEP + 1 ))
        END=\$(( (i + 1) * STEP ))
        [ \${END} -gt \${CHR_LEN} ] && END=\${CHR_LEN}

        if [ \${START} -le \${CHR_LEN} ]; then
            (
                samtools consensus \
                    -r "\${CHR_NAME}:\${START}-\${END}" \
                    -m simple \
                    -A \
                    -d ${min_depth} \
                    -H ${min_iupac} \
                    -f fasta \
                    "${bam}" | grep -v "^>" | tr -d '\\n' > "part_\${i}.txt"
            ) &
        fi
    done

    wait

    cat part_*.txt > full_seq.txt
    rm -f part_*.txt

    awk '{ gsub(/^N+/, "", \$0); gsub(/N+\$/, "", \$0); print \$0 }' full_seq.txt > seq_clean.txt
    rm -f full_seq.txt

    if [ ! -s seq_clean.txt ]; then
        echo "ERREUR : consensus vide pour ${sample_id} ${genotype}" >&2
        exit 1
    fi

    echo ">${sample_id}_Geno_${genotype}" > "${sample_id}_Geno_${genotype}.fasta"
    fold -w 70 seq_clean.txt >> "${sample_id}_Geno_${genotype}.fasta"
    rm -f seq_clean.txt
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