nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 04: COMPETITIVE ALIGNMENT & READ DISPATCH (04_PRECONSENSUS)
========================================================================================
*/

process COMPETITIVE_ALIGN {
    tag { "$sample_id" }
    publishDir { "${params.outdir}/04_PRECONSENSUS/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq)
    path ref_db
    val virus_name
    val tech_name

    output:
    tuple val(sample_id), path("reads_by_geno/*.fastq.gz"), path("reads_by_geno/*.ref.fasta"), emit: partitioned_reads, optional: true
    tuple val(sample_id), path("*_real_counts.tsv"), emit: real_counts
    tuple val(sample_id), path("*_competitive.bam"), path("*_competitive.bam.bai"), emit: bam_bai, optional: true

    script:
    """
    echo "=== Alignement competitif PRONAME pour ${sample_id} (${virus_name} - ${tech_name}) ==="
    echo "--> Reference retenue : ${ref_db.name}"

    # Alignement minimap2 avec filtrage MAPQ >= 10 direct dans le BAM
    minimap2 -ax map-ont -t ${task.cpus} --secondary=no "${ref_db}" "${trimmed_fastq}" | \\
        samtools view -b -F 2048 -q 10 | \\
        samtools sort -@ ${task.cpus} -o "${sample_id}_competitive.bam"

    samtools index "${sample_id}_competitive.bam"

    # Partitionnement et comptage
    04_dispatch_reads.py \\
        "${sample_id}_competitive.bam" \\
        "${ref_db}" \\
        "${sample_id}" \\
        "reads_by_geno" \\
        0.05 \\
        20
    """

    stub:
    """
    mkdir -p reads_by_geno
    touch "${sample_id}_competitive.bam" "${sample_id}_competitive.bam.bai"
    echo -e "cluster\\treads\\tpercentage\\nHBV_A2\\t1000\\t100.0" > "${sample_id}_real_counts.tsv"
    touch "reads_by_geno/${sample_id}_A2.fastq.gz"
    echo -e ">HBV_A2\\nACGT" > "reads_by_geno/${sample_id}_A2.ref.fasta"
    """
}