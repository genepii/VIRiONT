nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 13: SCREENING DES MUTATIONS ET INDELS CLINIQUES VHB
    Sortie : 
      - 13_MUTATION/<barcode>/raw/ (archive bioinformatique)
      - 13_MUTATION/<barcode>/filtered/ (variants >= 1%, DP > 100 pour indels)
      - 13_MUTATION/<barcode>/glims_alert/ (alertes cliniques variants & indels)
      - 13_MUTATION/ALL_SAMPLES/ (rapports consolidés du run)
========================================================================================
*/

// =====================================================================================
// A. SUBSTITUTIONS
// =====================================================================================
process SEARCH_HBV_MUTATIONS {
    tag { "${sample_id}_${genotype}" }
    publishDir { "${params.outdir}/13_MUTATION/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(vcf_input)
    val virus_name
    path mutation_tables_dir
    path ref_fasta

    output:
    path "raw/${sample_id}_${genotype}_raw_mutations.csv"                  , emit: sample_raw_variants
    path "filtered/${sample_id}_${genotype}_filtered_mutations.csv"        , emit: sample_filtered_variants
    path "glims_alert/${sample_id}_${genotype}_glims_variants_alerts.csv"  , emit: sample_glims_variants_alerts

    script:
    def freq_filtered  = params.freq_min_filtered != null ? params.freq_min_filtered : 1.0
    def min_glims_mino = params.glims_min_minor_freq ?: 8.0
    def win_pos        = params.window_pos != null ? params.window_pos : 0

    """
    mkdir -p raw filtered glims_alert

    RAW_CSV="raw/${sample_id}_${genotype}_raw_mutations.csv"
    FILT_CSV="filtered/${sample_id}_${genotype}_filtered_mutations.csv"
    GLIMS_CSV="glims_alert/${sample_id}_${genotype}_glims_variants_alerts.csv"

    if [ "${virus_name}" = "VHB" ] && [ -s "${vcf_input}" ]; then
        if gzip -t "${vcf_input}" 2>/dev/null; then
            gzip -dc "${vcf_input}" > "${sample_id}.vcf"
        else
            cp "${vcf_input}" "${sample_id}.vcf"
        fi

        if grep -v "^#" "${sample_id}.vcf" | grep -q '[^[:space:]]'; then
            Rscript ${projectDir}/bin/13_search_mutation.R \\
                --vcf "${sample_id}.vcf" \\
                --tables-dir "${mutation_tables_dir}" \\
                --ref-fasta "${ref_fasta}" \\
                --genotype "${genotype}" \\
                --freq-min 0 \\
                --window-pos ${win_pos} \\
                --sample-id "${sample_id}" \\
                --out-variants "\$RAW_CSV"

            Rscript ${projectDir}/bin/13_search_mutation.R \\
                --vcf "${sample_id}.vcf" \\
                --tables-dir "${mutation_tables_dir}" \\
                --ref-fasta "${ref_fasta}" \\
                --genotype "${genotype}" \\
                --freq-min ${freq_filtered} \\
                --window-pos ${win_pos} \\
                --sample-id "${sample_id}" \\
                --out-variants "\$FILT_CSV"
        fi
    fi

    Rscript ${projectDir}/bin/13_extract_glims_alerts.R \\
        --input "\$FILT_CSV" \\
        --min-freq-minor ${min_glims_mino} \\
        --sample-id "${sample_id}" \\
        --out "\$GLIMS_CSV"
    """

    stub:
    """
    mkdir -p raw filtered glims_alert
    touch "raw/${sample_id}_${genotype}_raw_mutations.csv"
    touch "filtered/${sample_id}_${genotype}_filtered_mutations.csv"
    touch "glims_alert/${sample_id}_${genotype}_glims_variants_alerts.csv"
    """
}

process COLLECT_MUTATION_REPORTS {
    publishDir "${params.outdir}/13_MUTATION/ALL_SAMPLES", mode: 'copy'

    input:
    path raw_csvs
    path filtered_csvs
    path glims_csvs

    output:
    path "raw/ALL_SAMPLES_raw_variants.csv"                       , emit: raw_variants_csv
    path "filtered/ALL_SAMPLES_filtered_variants.csv"             , emit: filtered_variants_csv
    path "glims_alert/ALL_SAMPLES_glims_variants_alerts.csv"       , emit: glims_variants_alerts_csv

    script:
    """
    mkdir -p raw filtered glims_alert

    FIRST_RAW=\$(echo "${raw_csvs}" | awk '{print \$1}')
    if [ -n "\$FIRST_RAW" ] && [ -s "\$FIRST_RAW" ]; then
        head -n 1 "\$FIRST_RAW" > raw/ALL_SAMPLES_raw_variants.csv
        for f in ${raw_csvs}; do
            [ -s "\$f" ] && tail -n +2 "\$f" >> raw/ALL_SAMPLES_raw_variants.csv
        done
    else
        touch raw/ALL_SAMPLES_raw_variants.csv
    fi

    FIRST_FILT=\$(echo "${filtered_csvs}" | awk '{print \$1}')
    if [ -n "\$FIRST_FILT" ] && [ -s "\$FIRST_FILT" ]; then
        head -n 1 "\$FIRST_FILT" > filtered/ALL_SAMPLES_filtered_variants.csv
        for f in ${filtered_csvs}; do
            [ -s "\$f" ] && tail -n +2 "\$f" >> filtered/ALL_SAMPLES_filtered_variants.csv
        done
    else
        touch filtered/ALL_SAMPLES_filtered_variants.csv
    fi

    FIRST_GLIMS=\$(echo "${glims_csvs}" | awk '{print \$1}')
    if [ -n "\$FIRST_GLIMS" ] && [ -s "\$FIRST_GLIMS" ]; then
        head -n 1 "\$FIRST_GLIMS" > glims_alert/ALL_SAMPLES_glims_variants_alerts.csv
        for f in ${glims_csvs}; do
            [ -s "\$f" ] && tail -n +2 "\$f" >> glims_alert/ALL_SAMPLES_glims_variants_alerts.csv
        done
    else
        touch glims_alert/ALL_SAMPLES_glims_variants_alerts.csv
    fi
    """

    stub:
    """
    mkdir -p raw filtered glims_alert
    touch raw/ALL_SAMPLES_raw_variants.csv
    touch filtered/ALL_SAMPLES_filtered_variants.csv
    touch glims_alert/ALL_SAMPLES_glims_variants_alerts.csv
    """
}

// =====================================================================================
// B. INDELS & FRAMESHIFTS
// =====================================================================================
process SEARCH_HBV_INDELS {
    tag { "${sample_id}_${genotype}" }
    publishDir { "${params.outdir}/13_MUTATION/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(vcf_input)
    val virus_name
    path mutation_tables_dir
    path ref_fasta

    output:
    path "raw/${sample_id}_${genotype}_raw_indels.csv"                 , emit: sample_raw_indels
    path "filtered/${sample_id}_${genotype}_filtered_indels.csv"       , emit: sample_filtered_indels
    path "glims_alert/${sample_id}_${genotype}_glims_indels_alerts.csv", emit: sample_glims_indels_alerts

    script:
    def freq_filtered  = params.freq_min_filtered      != null ? params.freq_min_filtered      : 1.0
    def dp_filtered    = params.indel_min_dp_filtered  != null ? params.indel_min_dp_filtered  : 101
    def alt_filtered   = params.indel_min_alt_filtered != null ? params.indel_min_alt_filtered : 5
    def min_glims_mino = params.glims_min_minor_freq   ?: 8.0
    def win_pos        = params.window_pos != null ? params.window_pos : 0

    """
    mkdir -p raw filtered glims_alert

    RAW_CSV="raw/${sample_id}_${genotype}_raw_indels.csv"
    FILT_CSV="filtered/${sample_id}_${genotype}_filtered_indels.csv"
    GLIMS_CSV="glims_alert/${sample_id}_${genotype}_glims_indels_alerts.csv"

    # En-tête adapté et standardisé pour les INDELs
    INDEL_HEADER="REFERENCE;REF_POS;CODON;GENE;GENOTYPE;INDEL_TYPE;INDEL_SIZE_BP;CONSEQUENCE;REF_ALLELE;ALT_ALLELE;INDEL_NAME;FREQ;STATUS;ALT_READS;TOTAL_DEPTH;warning"
    echo "\$INDEL_HEADER" > "\$RAW_CSV"
    echo "\$INDEL_HEADER" > "\$FILT_CSV"

    if [ "${virus_name}" = "VHB" ] && [ -s "${vcf_input}" ]; then
        if gzip -t "${vcf_input}" 2>/dev/null; then
            gzip -dc "${vcf_input}" > "${sample_id}.vcf"
        else
            cp "${vcf_input}" "${sample_id}.vcf"
        fi

        if grep -v "^#" "${sample_id}.vcf" | grep -q '[^[:space:]]'; then
            Rscript ${projectDir}/bin/13_search_indels.R \\
                --vcf "${sample_id}.vcf" \\
                --tables-dir "${mutation_tables_dir}" \\
                --ref-fasta "${ref_fasta}" \\
                --genotype "${genotype}" \\
                --freq-min 0 --min-dp 1 --min-alt 1 \\
                --window-pos ${win_pos} \\
                --sample-id "${sample_id}" \\
                --out-indels "\$RAW_CSV"

            Rscript ${projectDir}/bin/13_search_indels.R \\
                --vcf "${sample_id}.vcf" \\
                --tables-dir "${mutation_tables_dir}" \\
                --ref-fasta "${ref_fasta}" \\
                --genotype "${genotype}" \\
                --freq-min ${freq_filtered} --min-dp ${dp_filtered} --min-alt ${alt_filtered} \\
                --window-pos ${win_pos} \\
                --sample-id "${sample_id}" \\
                --out-indels "\$FILT_CSV"
        fi
    fi

    # Extraction ciblée GLIMS pour les Indels
    Rscript ${projectDir}/bin/13_extract_glims_alerts.R \\
        --input "\$FILT_CSV" \\
        --min-freq-minor ${min_glims_mino} \\
        --sample-id "${sample_id}" \\
        --out "\$GLIMS_CSV"
    """

    stub:
    """
    mkdir -p raw filtered glims_alert
    touch "raw/${sample_id}_${genotype}_raw_indels.csv"
    touch "filtered/${sample_id}_${genotype}_filtered_indels.csv"
    touch "glims_alert/${sample_id}_${genotype}_glims_indels_alerts.csv"
    """
}

process COLLECT_INDEL_REPORTS {
    publishDir "${params.outdir}/13_MUTATION/ALL_SAMPLES", mode: 'copy'

    input:
    path raw_csvs
    path filtered_csvs
    path glims_csvs

    output:
    path "raw/ALL_SAMPLES_raw_indels.csv"                    , emit: raw_indels_csv
    path "filtered/ALL_SAMPLES_filtered_indels.csv"          , emit: filtered_indels_csv
    path "glims_alert/ALL_SAMPLES_glims_indels_alerts.csv"    , emit: glims_indels_alerts_csv

    script:
    """
    mkdir -p raw filtered glims_alert

    FIRST_RAW=\$(echo "${raw_csvs}" | awk '{print \$1}')
    if [ -n "\$FIRST_RAW" ] && [ -s "\$FIRST_RAW" ]; then
        head -n 1 "\$FIRST_RAW" > raw/ALL_SAMPLES_raw_indels.csv
        for f in ${raw_csvs}; do
            [ -s "\$f" ] && tail -n +2 "\$f" >> raw/ALL_SAMPLES_raw_indels.csv
        done
    else
        touch raw/ALL_SAMPLES_raw_indels.csv
    fi

    FIRST_FILT=\$(echo "${filtered_csvs}" | awk '{print \$1}')
    if [ -n "\$FIRST_FILT" ] && [ -s "\$FIRST_FILT" ]; then
        head -n 1 "\$FIRST_FILT" > filtered/ALL_SAMPLES_filtered_indels.csv
        for f in ${filtered_csvs}; do
            [ -s "\$f" ] && tail -n +2 "\$f" >> filtered/ALL_SAMPLES_filtered_indels.csv
        done
    else
        touch filtered/ALL_SAMPLES_filtered_indels.csv
    fi

    FIRST_GLIMS=\$(echo "${glims_csvs}" | awk '{print \$1}')
    if [ -n "\$FIRST_GLIMS" ] && [ -s "\$FIRST_GLIMS" ]; then
        head -n 1 "\$FIRST_GLIMS" > glims_alert/ALL_SAMPLES_glims_indels_alerts.csv
        for f in ${glims_csvs}; do
            [ -s "\$f" ] && tail -n +2 "\$f" >> glims_alert/ALL_SAMPLES_glims_indels_alerts.csv
        done
    else
        touch glims_alert/ALL_SAMPLES_glims_indels_alerts.csv
    fi
    """

    stub:
    """
    mkdir -p raw filtered glims_alert
    touch raw/ALL_SAMPLES_raw_indels.csv
    touch filtered/ALL_SAMPLES_filtered_indels.csv
    touch glims_alert/ALL_SAMPLES_glims_indels_alerts.csv
    """
}