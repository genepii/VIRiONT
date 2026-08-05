/*
========================================================================================
    MODULE 16: SCREENING DES MUTATIONS VHB
========================================================================================
*/

process SEARCH_HBV_MUTATIONS {
    tag "$sample_id"
    publishDir path: { "${params.outdir}/10_MUTATIONS/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(vcf_gz)
    val virus_name
    path mutation_tables_dir

    output:
    tuple val(sample_id), path("*_vcf_variants.csv"), optional: true, emit: global_variants

    script:
    """
    if [ "${virus_name}" = "VHB" ]; then
        gunzip -c ${vcf_gz} > ${sample_id}_uncompressed.vcf
        
        16_search_mutation.R \
            ${sample_id}_uncompressed.vcf \
            ${mutation_tables_dir} \
            ${params.freq_min} \
            ${params.window_pos} \
            ${sample_id}_PC.csv ${sample_id}_BCP.csv ${sample_id}_DS.csv ${sample_id}_RT.csv \
            ${sample_id}_DPS1.csv ${sample_id}_DPS2.csv ${sample_id}_DHBx.csv ${sample_id}_C.csv \
            ${sample_id}_PC_F.csv ${sample_id}_BCP_F.csv ${sample_id}_DS_F.csv ${sample_id}_RT_F.csv \
            ${sample_id}_DPS1_F.csv ${sample_id}_DPS2_F.csv ${sample_id}_DHBx_F.csv ${sample_id}_C_F.csv \
            ${sample_id}_vcf_variants.csv
    else
        echo "Non VHB sample - skipping mutation search"
        touch ${sample_id}_vcf_variants.csv
    fi
    """
}