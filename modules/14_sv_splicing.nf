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
    path splice_catalog

    output:
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_sv.vcf")          , emit: vcf
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_splicing_table.csv"), emit: table
    path "${sample_id}_${genotype}_sniffles.log"                                           , emit: log

    script:
    """
    set -euo pipefail

    # Indexation de la référence si l'index n'est pas déjà présent
    if [ ! -f "${canonical_ref}.fai" ]; then
        samtools faidx "${canonical_ref}" 2>/dev/null || true
    fi

    # 1. Conversion éphémère N -> D dédiée uniquement à Sniffles2
    if command -v samtools >/dev/null 2>&1; then
        samtools view -h "${bam}" | \
            awk 'BEGIN{FS=OFS="\\t"} /^@/{print; next} {gsub(/N/, "D", \$6); print}' | \
            samtools view -bS - > sniffles_input.bam
        samtools index sniffles_input.bam
    else
        python3 -c "
import pysam
with pysam.AlignmentFile('${bam}', 'rb') as infile:
    with pysam.AlignmentFile('sniffles_input.bam', 'wb', template=infile) as outfile:
        for read in infile:
            if read.cigartuples:
                read.cigartuples = [(2 if op == 3 else op, length) for op, length in read.cigartuples]
            outfile.write(read)
pysam.index('sniffles_input.bam')
"
    fi

    # 2. Appel des variants structuraux avec Sniffles2
    sniffles \
        --input sniffles_input.bam \
        --vcf "${sample_id}_${genotype}_sv.vcf" \
        --reference "${canonical_ref}" \
        --mosaic \
        --mosaic-af-min 0.05 \
        --mosaic-af-max 1.0 \
        --mosaic-include-germline \
        --mosaic-qc-coverage-max-change-frac -1 \
        --min-alignment-length 100 \
        --minsupport 50 \
        --minsvlen 150 \
        --threads ${task.cpus} \
        --mapq 15 2>&1 | tee "${sample_id}_${genotype}_sniffles.log"

    # Nettoyage immédiat du BAM temporaire
    rm -f sniffles_input.bam sniffles_input.bam.bai

    # 3. Extraction, croisement avec le catalogue et génération au format CSV (séparateur virgule)
    awk '
    BEGIN {
        FS="\\t";
        OFS=",";
        print "SAMPLE", "GENOTYPE", "CHROM", "POS_START", "POS_END", "SPLICE_LEN", "SUPPORT_READS", "TOTAL_DEPTH", "VAF_PERCENT", "TRANSCRIPT", "TYPE"
    }
    # Lecture du catalogue d épissage (1er fichier passé en argument)
    NR == FNR {
        if (\$1 != "NAME") {
            key = \$4 "_" \$5;
            sp_name[key] = \$1;
            sp_type[key] = \$7;
        }
        next;
    }
    # Traitement du fichier VCF Sniffles2
    !/^#/ && \$7 == "PASS" && \$8 ~ /SVTYPE=DEL/ && \$8 ~ /(^|;)PRECISE(;|\$)/ && \$8 !~ /IMPRECISE/ {
        svlen = "NA"; end = "NA"; supp = "NA"; af = "0";
        if (match(\$8, /SVLEN=-?[0-9]+/)) { svlen = substr(\$8, RSTART+6, RLENGTH-6); gsub("-", "", svlen) }
        if (match(\$8, /END=[0-9]+/))      { end = substr(\$8, RSTART+4, RLENGTH-4) }
        if (match(\$8, /SUPPORT=[0-9]+/))  { supp = substr(\$8, RSTART+8, RLENGTH-8) }
        if (match(\$8, /AF=[0-9.]+/))      { af = substr(\$8, RSTART+3, RLENGTH-3) }

        if ((svlen + 0) > 2200) next;

        split(\$9, fmt_keys, ":");
        split(\$10, fmt_vals, ":");
        dr = 0; dv = 0;
        for (i=1; i<=length(fmt_keys); i++) {
            if (fmt_keys[i] == "DR") dr = fmt_vals[i] + 0;
            if (fmt_keys[i] == "DV") dv = fmt_vals[i] + 0;
        }
        total_cov = dr + dv;
        af_pct = sprintf("%.2f", af * 100);

        vcf_key = \$2 "_" end;
        transcript = "UNKNOWN";
        tr_type = "NOVEL";

        if (vcf_key in sp_name) {
            transcript = sp_name[vcf_key];
            tr_type = sp_type[vcf_key];
        } else if ((svlen + 0) >= 1200 && (svlen + 0) <= 1300) {
            transcript = "SP1_HBV";
            tr_type = "CANONICAL";
        }

        # Filtre de confiance clinique (SUPPORT >= 50 et VAF >= 5.0%)
        if ((supp + 0) >= 50 && (af + 0) >= 0.05) {
            print "${sample_id}", "${genotype}", \$1, \$2, end, svlen, supp, total_cov, af_pct, transcript, tr_type
        }
    }' "${splice_catalog}" "${sample_id}_${genotype}_sv.vcf" > "${sample_id}_${genotype}_splicing_table.csv"
    """

    stub:
    """
    touch "${sample_id}_${genotype}_sv.vcf"
    echo "SAMPLE,GENOTYPE,CHROM,POS_START,POS_END,SPLICE_LEN,SUPPORT_READS,TOTAL_DEPTH,VAF_PERCENT,TRANSCRIPT,TYPE" > "${sample_id}_${genotype}_splicing_table.csv"
    echo "=== Sniffles stub ===" > "${sample_id}_${genotype}_sniffles.log"
    """
}

/*
========================================================================================
    TABLEAU RÉCAPITULATIF GLOBAL DES ÉPISSAGES (CSV)
========================================================================================
*/
process COLLECT_SPLICING_REPORTS {
    publishDir "${params.outdir}/14_SV_SPLICING", mode: 'copy'

    input:
    path(tables)

    output:
    path "global_splicing_summary.csv", emit: global_summary

    script:
    """
    set -euo pipefail

    awk '
    BEGIN { FS=","; OFS="," }
    FNR == 1 {
        if (NR == 1) {
            print \$0
        }
        next
    }
    {
        print \$0
    }' ${tables} > global_splicing_summary.csv
    """

    stub:
    """
    echo "SAMPLE,GENOTYPE,CHROM,POS_START,POS_END,SPLICE_LEN,SUPPORT_READS,TOTAL_DEPTH,VAF_PERCENT,TRANSCRIPT,TYPE" > global_splicing_summary.csv
    """
}