nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 12: SCREENING DES MUTATIONS CLINIQUES VHB/VHD (12_MUTATION)
========================================================================================
*/

process SEARCH_HBV_MUTATIONS {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/12_MUTATION/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(vcf_gz)
    val virus_name
    path mutation_tables_dir

    output:
    path "filtered/${sample_id}_${genotype}_filtered_mutations.csv", emit: sample_variants
    path "all_results/*.csv"                                       , optional: true
    path "filtered/*.csv"                                          , optional: true

    script:
    """
    echo "=== Screening Clinique des Mutations pour ${sample_id} - Génotype ${genotype} ==="

    mkdir -p all_results filtered

    if [ "${virus_name}" = "VHB" ] && [ -s "${vcf_gz}" ]; then
        zcat "${vcf_gz}" > "${sample_id}.vcf"

        if grep -q -v "^#" "${sample_id}.vcf"; then
            Rscript ${projectDir}/bin/12_search_mutation.R \
                --vcf "${sample_id}.vcf" \
                --tables-dir "${mutation_tables_dir}" \
                --genotype "${genotype}" \
                --freq-min ${params.freq_min ?: 0.08} \
                --window-pos ${params.window_pos ?: 10} \
                --sample-id "${sample_id}" \
                --out-variants "filtered/${sample_id}_${genotype}_filtered_mutations.csv" || true

            [ -f "filtered/${sample_id}_${genotype}_filtered_mutations.csv" ] || touch "filtered/${sample_id}_${genotype}_filtered_mutations.csv"
        else
            touch "filtered/${sample_id}_${genotype}_filtered_mutations.csv"
        fi
    else
        touch "filtered/${sample_id}_${genotype}_filtered_mutations.csv"
    fi
    """
}

process COLLECT_MUTATION_REPORTS {
    publishDir "${params.outdir}/12_MUTATION", mode: 'copy'

    input:
    path variant_csvs

    output:
    path "ALL_SAMPLES_vcf_variants.csv", emit: global_variants_csv

    script:
    """
    echo "=== Consolidation globale du rapport de mutations ==="

    first=1
    for f in ${variant_csvs}; do
        if [ -s "\$f" ]; then
            if [ \$first -eq 1 ]; then
                cat "\$f" > ALL_SAMPLES_vcf_variants.csv
                first=0
            else
                tail -n +2 "\$f" >> ALL_SAMPLES_vcf_variants.csv
            fi
        fi
    done

    if [ \$first -eq 1 ]; then
        echo "REFERENCE;GENE;GENOTYPE;REF_POS;Position_EcoR1;REF;total_count;base_status;base;count;Mutation_name;freq;ALT;NUM_CODON;POS_TYPE;REF_AA;ALT_AA;warning" > ALL_SAMPLES_vcf_variants.csv
    fi
    """
}