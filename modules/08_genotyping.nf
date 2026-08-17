nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 08: GÉNOTYPAGE, ORIENTATION & FILTRAGE CLINIQUE (08_GENOTYPING)
========================================================================================
*/
process GENOTYPING {
    publishDir "${params.outdir}/08_GENOTYPING", mode: 'copy', pattern: '{GLOBAL_Genotyping_Report.pdf,SUMMARY_Multi_Infection.tsv,validated_consensus_all.fasta,renamed_consensus/*.fasta}'
    publishDir "${params.outdir}/00_SUPDATA", mode: 'copy', pattern: 'DB/*'
    input:
    path 'consensus_inputs/*'  // FASTA consensus émis par Medaka
    path ref_fasta             // FASTA de référence
    path 'csvs/*'              // Comptages réels exhaustifs (samtools idxstats, ex-CSVs AmpliconSorter)
    output:
    path "GLOBAL_Genotyping_Report.pdf"      , emit: global_pdf
    path "SUMMARY_Multi_Infection.tsv"        , emit: summary_tsv
    path "validated_consensus_all.fasta"      , emit: validated_all_fasta
    path "renamed_consensus/*.fasta"          , emit: genotyped_fastas, optional: true
    path "DB/*"                               , emit: db_files
    script:
    def base_name = ref_fasta.baseName
    """
    mkdir -p DB
    if [ ! -f "DB/${base_name}.nhr" ]; then
        makeblastdb -in "${ref_fasta}" -dbtype nucl -out "DB/${base_name}"
    fi
    # Exécution du script Python de traitement unifié (avec ref_fasta pour les longueurs attendues)
    python3 ${projectDir}/bin/08_genotype_and_filter.py consensus_inputs csvs "DB/${base_name}" ${params.mi_cutoff} "${ref_fasta}"
    # Génération du rapport graphique PDF
    Rscript ${projectDir}/bin/08_generate_report.R SUMMARY_Multi_Infection.tsv ${params.mi_cutoff} GLOBAL_Genotyping_Report.pdf
    """
}