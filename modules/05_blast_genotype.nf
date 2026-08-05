nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 05: BLAST & GÉNOTYPAGE
========================================================================================
*/

process BLAST_GENOTYPE {
    tag "$sample_id"
    publishDir path: { "${params.outdir}/05_GENOTYPING/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fastas)
    path db_dir
    val ref_filename

    output:
    tuple val(sample_id), path("*_validated_genotypes.tsv"), emit: validated_tsv
    tuple val(sample_id), path("*_master_consensus.fasta"), optional: true, emit: master_fasta
    tuple val(sample_id), path("*.pdf"), optional: true, emit: report_pdf

    script:
    def db_name = ref_filename.replaceFirst(/\.[^.]+$/, "")
    """
    # 1. Alignement BLAST des clusters contre la base de référence
    blastn -db ${db_dir}/${db_name} \
           -query ${consensus_fastas} \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop stitle" \
           -out ${sample_id}_blast_results.tsv

    # 2. Agrégation et filtrage des génotypes via le script R
    5_generate_report.R ${sample_id}_blast_results.tsv ${params.mi_cutoff} ${sample_id}_validated_genotypes.tsv
    """
}