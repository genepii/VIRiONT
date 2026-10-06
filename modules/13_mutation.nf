nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 13: SCREENING DES MUTATIONS ET INDELS CLINIQUES VHB (13_MUTATION)
========================================================================================
*/

// =====================================================================================
// A. SUBSTITUTIONS (CODE D'ORIGINE SANCTUARISÉ)
// =====================================================================================
process SEARCH_HBV_MUTATIONS {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/13_MUTATION/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(vcf_gz)
    val virus_name
    path mutation_tables_dir
    path ref_fasta

    output:
    path "filtered/${sample_id}_${genotype}_raw_mutations.csv"   , emit: sample_raw_variants
    path "filtered/${sample_id}_${genotype}_major_mutations.csv" , emit: sample_major_variants

    script:
    def freq_cutoff = params.freq_min != null ? params.freq_min : 0.0
    def win_pos     = params.window_pos != null ? params.window_pos : 0

    """
    mkdir -p filtered

    HEADER="REFERENCE;REF_POS;REF;total_count;base_status;base;count;NUM_CODON;POS_TYPE;REF_AA;ALT_AA;Mutation_name;freq;warning;GENE;GENOTYPE"

    RAW_CSV="filtered/${sample_id}_${genotype}_raw_mutations.csv"
    MAJOR_CSV="filtered/${sample_id}_${genotype}_major_mutations.csv"

    echo "\$HEADER" > "\$RAW_CSV"
    echo "\$HEADER" > "\$MAJOR_CSV"

    if [ "${virus_name}" = "VHB" ] && [ -s "${vcf_gz}" ]; then
        zcat "${vcf_gz}" > "${sample_id}.vcf"

        if grep -v "^#" "${sample_id}.vcf" | grep -q '[^[:space:]]'; then
            Rscript ${projectDir}/bin/13_search_mutation.R \\
                --vcf "${sample_id}.vcf" \\
                --tables-dir "${mutation_tables_dir}" \\
                --ref-fasta "${ref_fasta}" \\
                --genotype "${genotype}" \\
                --freq-min ${freq_cutoff} \\
                --window-pos ${win_pos} \\
                --sample-id "${sample_id}" \\
                --out-variants "\$RAW_CSV"
        fi
    fi

    awk -F';' '\$5 == "majo"' "\$RAW_CSV" >> "\$MAJOR_CSV"
    """
}

process COLLECT_MUTATION_REPORTS {
    publishDir "${params.outdir}/13_MUTATION", mode: 'copy'

    input:
    path raw_variant_csvs

    output:
    path "ALL_SAMPLES_raw_variants.csv"   , emit: raw_variants_csv
    path "ALL_SAMPLES_major_variants.csv" , emit: major_variants_csv

    script:
    """
    HEADER="REFERENCE;REF_POS;REF;total_count;base_status;base;count;NUM_CODON;POS_TYPE;REF_AA;ALT_AA;Mutation_name;freq;warning;GENE;GENOTYPE"

    echo "\$HEADER" > ALL_SAMPLES_raw_variants.csv
    echo "\$HEADER" > ALL_SAMPLES_major_variants.csv

    for f in ${raw_variant_csvs}; do
        if [ -s "\$f" ]; then
            tail -n +2 "\$f" >> ALL_SAMPLES_raw_variants.csv
        fi
    done

    awk -F';' '\$5 == "majo"' ALL_SAMPLES_raw_variants.csv >> ALL_SAMPLES_major_variants.csv
    """
}

// =====================================================================================
// B. INDELS & FRAMESHIFTS (PROCESS DÉDIÉ ET ÉTANCHE)
// =====================================================================================
process SEARCH_HBV_INDELS {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/13_MUTATION/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(vcf_gz)
    val virus_name
    path mutation_tables_dir
    path ref_fasta

    output:
    path "filtered/${sample_id}_${genotype}_raw_indels.csv"   , emit: sample_raw_indels
    path "filtered/${sample_id}_${genotype}_major_indels.csv" , emit: sample_major_indels

    script:
    def freq_cutoff = params.freq_min != null ? params.freq_min : 0.0
    def win_pos     = params.window_pos != null ? params.window_pos : 0

    """
    mkdir -p filtered

    HEADER="REFERENCE;REF_POS;REF;total_count;base_status;base;count;NUM_CODON;POS_TYPE;REF_AA;ALT_AA;Mutation_name;freq;warning;GENE;GENOTYPE"

    RAW_CSV="filtered/${sample_id}_${genotype}_raw_indels.csv"
    MAJOR_CSV="filtered/${sample_id}_${genotype}_major_indels.csv"

    echo "\$HEADER" > "\$RAW_CSV"
    echo "\$HEADER" > "\$MAJOR_CSV"

    if [ "${virus_name}" = "VHB" ] && [ -s "${vcf_gz}" ]; then
        zcat "${vcf_gz}" > "${sample_id}.vcf"

        if grep -v "^#" "${sample_id}.vcf" | grep -q '[^[:space:]]'; then
            Rscript ${projectDir}/bin/13_search_indels.R \\
                --vcf "${sample_id}.vcf" \\
                --tables-dir "${mutation_tables_dir}" \\
                --ref-fasta "${ref_fasta}" \\
                --genotype "${genotype}" \\
                --freq-min ${freq_cutoff} \\
                --window-pos ${win_pos} \\
                --sample-id "${sample_id}" \\
                --out-indels "\$RAW_CSV"
        fi
    fi

    awk -F';' '\$5 == "majo"' "\$RAW_CSV" >> "\$MAJOR_CSV"
    """
}

process COLLECT_INDEL_REPORTS {
    publishDir "${params.outdir}/13_MUTATION", mode: 'copy'

    input:
    path raw_indel_csvs

    output:
    path "ALL_SAMPLES_raw_indels.csv"   , emit: raw_indels_csv
    path "ALL_SAMPLES_major_indels.csv" , emit: major_indels_csv

    script:
    """
    HEADER="REFERENCE;REF_POS;REF;total_count;base_status;base;count;NUM_CODON;POS_TYPE;REF_AA;ALT_AA;Mutation_name;freq;warning;GENE;GENOTYPE"

    echo "\$HEADER" > ALL_SAMPLES_raw_indels.csv
    echo "\$HEADER" > ALL_SAMPLES_major_indels.csv

    for f in ${raw_indel_csvs}; do
        if [ -s "\$f" ]; then
            tail -n +2 "\$f" >> ALL_SAMPLES_raw_indels.csv
        fi
    done

    awk -F';' '\$5 == "majo"' ALL_SAMPLES_raw_indels.csv >> ALL_SAMPLES_major_indels.csv
    """
}