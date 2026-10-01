nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 14: DÉTECTION DES VARIANTS STRUCTURAUX / ÉPISSAGE VHB (14_SV_SPLICING)
========================================================================================
*/
process SV_SPLICING_VHB {
    tag "${sample_id}_${genotype}"
    publishDir path: { "${params.outdir}/14_SV_SPLICING/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), path(canonical_ref), path(bam), path(bai)

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_sv.vcf")          , emit: vcf
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_splicing_table.tsv"), emit: table
    path "${sample_id}_${genotype}_sniffles.log"                                           , emit: log

    script:
    """
    set -euo pipefail

    # Indexation de la référence si l'index n'est pas déjà présent
    if [ ! -f "${canonical_ref}.fai" ]; then
        samtools faidx "${canonical_ref}" 2>/dev/null || true
    fi

    # 1. Appel des variants structuraux avec Sniffles2 (seuil minimal rehaussé à 5 % et 50 reads)
    sniffles \\
        --input "${bam}" \\
        --vcf "${sample_id}_${genotype}_sv.vcf" \\
        --reference "${canonical_ref}" \\
        --mosaic \\
        --mosaic-af-min 0.05 \\
        --mosaic-af-max 1.0 \\
        --mosaic-include-germline \\
        --mosaic-qc-coverage-max-change-frac -1 \\
        --min-alignment-length 100 \\
        --minsupport 50 \\
        --minsvlen 150 \\
        --threads ${task.cpus} \\
        --mapq 15 2>&1 | tee "${sample_id}_${genotype}_sniffles.log"

    # 2. Extraction stricte avec filtre anti-bruit (SUPPORT >= 50 reads et VAF >= 5.0%)
    awk -F'\\t' '
    BEGIN {
        OFS="\\t";
        print "SAMPLE", "GENOTYPE", "CHROM", "POS_START", "POS_END", "SPLICE_LEN", "SUPPORT_READS", "TOTAL_DEPTH", "VAF_PERCENT", "STATUS"
    }
    !/^#/ && \$7 == "PASS" && \$8 ~ /SVTYPE=DEL/ && \$8 ~ /PRECISE/ {
        svlen = "NA"; end = "NA"; supp = "NA"; af = "0";
        if (match(\$8, /SVLEN=-?[0-9]+/)) { svlen = substr(\$8, RSTART+6, RLENGTH-6); gsub("-", "", svlen) }
        if (match(\$8, /END=[0-9]+/))      { end = substr(\$8, RSTART+4, RLENGTH-4) }
        if (match(\$8, /SUPPORT=[0-9]+/))  { supp = substr(\$8, RSTART+8, RLENGTH-8) }
        if (match(\$8, /AF=[0-9.]+/))      { af = substr(\$8, RSTART+3, RLENGTH-3) }

        split(\$9, fmt_keys, ":");
        split(\$10, fmt_vals, ":");
        dr = 0; dv = 0;
        for (i=1; i<=length(fmt_keys); i++) {
            if (fmt_keys[i] == "DR") dr = fmt_vals[i] + 0;
            if (fmt_keys[i] == "DV") dv = fmt_vals[i] + 0;
        }
        total_cov = dr + dv;
        af_pct = sprintf("%.2f", af * 100);

        status = "ALTERNATIVE_SPLICING";
        if (svlen >= 1200 && svlen <= 1300) {
            status = "SP1_CANONICAL";
        }

        # Filtre anti-bruit strict : VAF >= 5% et au moins 50 reads
        if ((supp + 0) >= 50 && (af + 0) >= 0.05) {
            print "${sample_id}", "${genotype}", \$1, \$2, end, svlen, supp, total_cov, af_pct, status
        }
    }' "${sample_id}_${genotype}_sv.vcf" > "${sample_id}_${genotype}_splicing_table.tsv"
    """

    stub:
    """
    touch "${sample_id}_${genotype}_sv.vcf"
    echo -e "SAMPLE\\tGENOTYPE\\tCHROM\\tPOS_START\\tPOS_END\\tSPLICE_LEN\\tSUPPORT_READS\\tTOTAL_DEPTH\\tVAF_PERCENT\\tSTATUS" > "${sample_id}_${genotype}_splicing_table.tsv"
    echo "=== Sniffles stub ===" > "${sample_id}_${genotype}_sniffles.log"
    """
}

/*
========================================================================================
    TABLEAU RÉCAPITULATIF GLOBAL DES ÉPISSAGES (À LA RACINE DE 14_SV_SPLICING)
========================================================================================
*/
process COLLECT_SPLICING_REPORTS {
    publishDir "${params.outdir}/14_SV_SPLICING", mode: 'copy'

    input:
    path(tables)

    output:
    path "global_splicing_summary.tsv", emit: global_summary

    script:
    """
    set -euo pipefail

    awk '
    FNR == 1 {
        if (NR == 1) {
            print \$0
        }
        next
    }
    {
        print \$0
    }' ${tables} > global_splicing_summary.tsv
    """

    stub:
    """
    echo -e "SAMPLE\\tGENOTYPE\\tCHROM\\tPOS_START\\tPOS_END\\tSPLICE_LEN\\tSUPPORT_READS\\tTOTAL_DEPTH\\tVAF_PERCENT\\tSTATUS" > global_splicing_summary.tsv
    """
}