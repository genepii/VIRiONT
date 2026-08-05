/*
========================================================================================
    MODULE 10: EXTRACTION INDIVIDUELLE DES MÉTRIQUES QC (08_QC_ANALYSIS)
========================================================================================
*/

process EXTRACT_SAMPLE_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/08_QC_ANALYSIS/temp", mode: 'copy'

    input:
    tuple val(sample_id), path(raw_fastq), path(dehost_fastq), path(trim_fastq), path(validated_tsv), path(bam), path(bai), path(vcf)

    output:
    path "${sample_id}_raw.csv", emit: raw_csv
    path "${sample_id}_dehost.csv", emit: dehost_csv
    path "${sample_id}_trimm.csv", emit: trimm_csv
    path "${sample_id}_geno.csv", emit: geno_csv

    script:
    """
    sample_raw="${sample_id}_raw.csv"
    sample_dehost="${sample_id}_dehost.csv"
    sample_trimm="${sample_id}_trimm.csv"
    sample_geno="${sample_id}_geno.csv"

    > "\$sample_raw"
    > "\$sample_dehost"
    > "\$sample_trimm"
    > "\$sample_geno"

    # 1. RAW
    if [ -f "${raw_fastq}" ] && [ -s "${raw_fastq}" ]; then
        gzip -dc ${raw_fastq} | awk -v sid="${sample_id}" '(NR%4==2){print length(\$0)";"sid";01_RAW;NONE;NA;NA;NA;NA"}' >> "\$sample_raw" || true
    fi

    # 2. DEHOSTING
    if [ -f "${dehost_fastq}" ] && [ -s "${dehost_fastq}" ]; then
        gzip -dc ${dehost_fastq} | awk -v sid="${sample_id}" '(NR%4==2){print length(\$0)";"sid";02_DEHOSTING;NONE;NA;NA;NA;NA"}' >> "\$sample_dehost" || true
    fi

    # 3. FILTERED_TRIMMED
    if [ -f "${trim_fastq}" ] && [ -s "${trim_fastq}" ]; then
        gzip -dc ${trim_fastq} | awk -v sid="${sample_id}" '(NR%4==2){print length(\$0)";"sid";03_FILTERED_TRIMMED;NONE;NA;NA;NA;NA"}' >> "\$sample_trimm" || true
    fi

    # Métriques d'alignement BAM et variants VCF
    mean_depth="NA"
    if [ -f "${bam}" ] && [ -s "${bam}" ]; then
        mean_depth=\$(samtools depth "${bam}" 2>/dev/null | awk '{sum+=\$3; cnt++} END {if(cnt>0) printf "%.1f", sum/cnt; else print "NA"}' || echo "NA")
    fi

    clair3_vars="NA"
    if [ -f "${vcf}" ] && [ -s "${vcf}" ]; then
        if grep -q -v "^#" "${vcf}"; then
            clair3_vars=\$(grep -v "^#" "${vcf}" | grep -cw "PASS" || true)
            if [ -z "\$clair3_vars" ] || [ "\$clair3_vars" -eq 0 ]; then
                clair3_vars=\$(grep -v -c "^#" "${vcf}" || echo "0")
            fi
        else
            clair3_vars=0
        fi
    fi

    # 4. GENOTYPING
    if [ -f "${validated_tsv}" ] && [ -f "${trim_fastq}" ] && [ -s "${trim_fastq}" ]; then
        temp_lengths=\$(mktemp)
        gzip -dc ${trim_fastq} | awk '(NR%4==2){print length(\$0)}' > "\$temp_lengths"

        grep -w "VALIDATED" "${validated_tsv}" | while IFS=\$'\\t' read -r s_id geno best_cluster total_reads ratio pident length status strand; do
            geno=\$(echo "\$geno" | tr -d '\\r\\n ')
            pident=\$(echo "\$pident" | tr -d '\\r\\n ')
            total_reads=\$(echo "\$total_reads" | tr -d '\\r\\n ')

            awk -v sid="${sample_id}" -v g="\$geno" -v p="\$pident" -v d="\$mean_depth" -v v="\$clair3_vars" -v a="\$total_reads" \
                '{print \$1";"sid";05_GENOTYPING;"g";"p";"d";"v";"a}' "\$temp_lengths" >> "\$sample_geno" || true
        done
        rm -f "\$temp_lengths"
    fi
    """
}