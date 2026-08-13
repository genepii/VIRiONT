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
    # 1. RECUPERATION DU CSV POUR CORRESPONDANCE
    # ==========================================
    CSV_FILE=\$(find -L "${fastq_dir}" -maxdepth 1 -type f -name "*.csv" | head -n 1)

    # Déclaration des variables de résultat
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

        # Extraction du numéro de barcode pour rechercher dans le CSV
        num_barcode=\$(echo "\$bname" | grep -oE '[0-9]+' | sed 's/^0*//')
        formatted_num=\$(printf "%02d" "\$num_barcode")

        # Recherche du nom complet de l'échantillon dans le CSV s'il existe
        full_sample_name="barcode_\${formatted_num}"
        if [ -f "\$CSV_FILE" ]; then
            sample_match=\$(awk -F'[;,]' -v b="\${num_barcode}" -v fb="\${formatted_num}" '
                NR>1 {
                    sample=\$1; gsub(/^[ \t\r\n"]+|[ \t\r\n"]+\$/, "", sample);
                    barcode=\$2; gsub(/^[ \t\r\n"]+|[ \t\r\n"]+\$/, "", barcode);
                    if (barcode ~ b || barcode ~ fb) {
                        print "barcode_" fb "_" sample;
                        exit;
                    }
                }
            ' "\$CSV_FILE")
            
            if [ -n "\$sample_match" ]; then
                full_sample_name="\$sample_match"
            fi
        fi

        # Vérification des fichiers et de la volumétrie
        all_files_count=\$(find "\$dir" -maxdepth 1 -type f | wc -l)
        fastq_count=\$(find "\$dir" -maxdepth 1 -type f \\( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \\) | wc -l)
        dir_size_kb=\$(du -sk "\$dir" | awk '{print \$1}')

        if [ "\$all_files_count" -eq 0 ]; then
            EMPTY_BARCODES="\${EMPTY_BARCODES}\${full_sample_name}\\n"
        elif [ "\$fastq_count" -gt 0 ]; then
            if [ "\$dir_size_kb" -lt 50 ]; then
                EMPTY_BARCODES="\${EMPTY_BARCODES}\${full_sample_name} (volumétrie négligeable : \${dir_size_kb} Ko)\\n"
            else
                VALID_BARCODES="\${VALID_BARCODES}\${full_sample_name}\\n"
                
                # Vérification de l'intégrité gzip
                for gz in \$(find "\$dir" -maxdepth 1 -type f -name "*.gz"); do
                    if ! gzip -t "\$gz" 2>/dev/null; then
                        CORRUPTED_FILES="\${CORRUPTED_FILES}\${gz}\\n"
                    fi
                done
            fi
        else
            NON_FASTQ_BARCODES="\${NON_FASTQ_BARCODES}\${full_sample_name}\\n"
        fi
    done

    # Génération du fichier fastq_content.txt
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