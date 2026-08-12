nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 00: INITIALISATION & MÉTADONNÉES DU RUN (00_INIT_LOGS)
========================================================================================
*/

process GENERATE_RUN_METADATA {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fastq_dir
    path target_ref
    val min_len
    val max_len

    output:
    path "fastq_content.txt", emit: fastq_log
    path "param_file.txt"   , emit: param_log

    script:
    """
    #!/usr/bin/env bash

    # ==========================================
    # 1. GENERATION DE FASTQ_CONTENT.TXT
    # ==========================================
    VALID_BARCODES=""
    NON_FASTQ_BARCODES=""
    EMPTY_BARCODES=""
    CORRUPTED_FILES=""

    for dir in \$(find -L "${fastq_dir}" -maxdepth 1 -mindepth 1 -type d | sort); do
        bname=\$(basename "\$dir")
        
        # Ignorer les dossiers hors barcodes
        if [[ ! "\$bname" =~ barcode[0-9]+ ]]; then
            continue
        fi

        # Vérification des fichiers
        all_files_count=\$(find "\$dir" -maxdepth 1 -type f | wc -l)
        fastq_count=\$(find "\$dir" -maxdepth 1 -type f \\( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \\) | wc -l)
        
        # Calcul de la taille du dossier en Ko
        dir_size_kb=\$(du -sk "\$dir" | awk '{print \$1}')

        if [ "\$all_files_count" -eq 0 ]; then
            EMPTY_BARCODES="\${EMPTY_BARCODES}\${bname}\\n"
        elif [ "\$fastq_count" -gt 0 ]; then
            # Si le dossier contient moins de 50 Ko de données FASTQ (témoins d'extraction vides / bruit)
            if [ "\$dir_size_kb" -lt 50 ]; then
                EMPTY_BARCODES="\${EMPTY_BARCODES}\${bname} (volumétrie négligeable : \${dir_size_kb} Ko)\\n"
            else
                VALID_BARCODES="\${VALID_BARCODES}\${bname}\\n"
                
                # Vérification de l'intégrité des fichiers gzip
                for gz in \$(find "\$dir" -maxdepth 1 -type f -name "*.gz"); do
                    if ! gzip -t "\$gz" 2>/dev/null; then
                        CORRUPTED_FILES="\${CORRUPTED_FILES}\${gz}\\n"
                    fi
                done
            fi
        else
            NON_FASTQ_BARCODES="\${NON_FASTQ_BARCODES}\${bname}\\n"
        fi
    done

    cat << EOF > fastq_content.txt
##########################
##### FASTQ ANALYSIS #####
##########################
barcode repository containing fastq/gz files and used for analysis:
\$(echo -e "\$VALID_BARCODES" | sed '/^\$/d')

##########################
barcode repository containing other files than fastq/gz and ignored for analysis:
\$(echo -e "\$NON_FASTQ_BARCODES" | sed '/^\$/d')

##########################
list of problematic files:
\$(echo -e "\$CORRUPTED_FILES" | sed '/^\$/d')

##########################
empty barcode repositories and ignored for analysis:
\$(echo -e "\$EMPTY_BARCODES" | sed '/^\$/d')
EOF

    # ==========================================
    # 2. GENERATION DE PARAM_FILE.TXT
    # ==========================================
    cat << EOF > param_file.txt
######################
#### PARAMS USED #####
######################
Data repository                         : ${params.fastq_dir}/
Result repository                       : ${params.outdir}/
Reference database used                 : ${target_ref.name}
Read minlength                          : ${min_len} bp
Read maxlength                          : ${max_len} bp
Quality filtering (Q-score)             : 10
Min reads threshold                     : ${params.min_reads}
Max reads subsampling                   : ${params.max_reads}
Multi-infection cutoff                  : ${params.mi_cutoff}%
Medaka model used                       : ${params.medaka_model}
Clair3 model used                       : ${params.clair3_model}
Hostile index used                      : ${params.hostile_index}
Nextflow engine version                 : ${workflow.nextflow.version}
Pipeline execution run name             : ${workflow.runName}
Command line executed                   : ${workflow.commandLine}
EOF
    """
}