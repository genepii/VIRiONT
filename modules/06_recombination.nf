nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 06: DÉTECTION DES RECOMBINAISONS (jpHMM)
========================================================================================
*/

process SPLIT_VALIDATED_FASTA {
    label 'main_container'

    input:
    path validated_multi_fasta

    output:
    path "single_fastas/*.fasta", emit: single_fastas

    script:
    """
    mkdir -p single_fastas
    awk '
        /^>/ {
            # Nettoie le header : prend le premier token sans le ">" et sans les metadata (reads=...)
            sub(/^>/, "", \$1)
            sample_id = \$1
            outfile = "single_fastas/" sample_id ".fasta"
            print ">" sample_id > outfile
            next
        }
        {
            if (outfile != "") print \$0 >> outfile
        }
    ' "${validated_multi_fasta}"
    """
}

process RECOMBINATION_JPHMM {
    tag "${sample_id}"
    label 'main_container'

    publishDir path: { "${params.outdir}/06_recombination/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta)
    path canonical_ecori_ref

    output:
    tuple val(sample_id), path("jphmm_output"),                            emit: raw_dir
    tuple val(sample_id), path("${sample_id}_recombination_summary.tsv"),  emit: summary

    script:
    """
    mkdir -p jphmm_output

    # 1. Rotation circulaire sur la position 1 EcoRI canonique
    06_rotate_hbv_ecori.py "${consensus_fasta}" "${canonical_ecori_ref}" "${sample_id}_rotated.fasta"

    # 2. Inférence jpHMM en mode circulaire (-C) sur le génome recalé
    jpHMM \\
        -s "${sample_id}_rotated.fasta" \\
        -v HBV \\
        -C \\
        -P /opt/jphmm/priors/ \\
        -I /opt/jphmm/input/ \\
        -o jphmm_output/

    # 3. Parsing et formatage du diagnostic clinique
    06_parse_jphmm.py jphmm_output/ ${sample_id}
    """
}

process COLLECT_RECOMBINATION_SUMMARIES {
    label 'main_container'

    publishDir "${params.outdir}/06_recombination", mode: 'copy'

    input:
    path summaries

    output:
    path "recombination_global_summary.tsv", emit: global_summary

    script:
    """
    first=\$(echo "${summaries}" | tr ' ' '\\n' | head -n 1)
    head -n 1 "\$first" > recombination_global_summary.tsv

    for f in ${summaries}; do
        tail -n +2 "\$f"
    done | sort -V >> recombination_global_summary.tsv
    """
}