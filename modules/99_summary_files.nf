nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 99: GÉNÉRATION DES FICHIERS DE CONFIGURATION ET RAPPORT RACINE
========================================================================================
*/

process GENERATE_SUMMARY_FILES {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path input_dir_path
    path ref_db_path

    output:
    path "fastq_content.txt", emit: fastq_content
    path "param_file.txt",   emit: param_file

    script:
    """
    echo "=== Génération de fastq_content.txt et param_file.txt ==="

    # 1. Génération de fastq_content.txt
    cat << 'EOF' > fastq_content.txt
##########################
##### FASTQ ANALYSIS #####
##########################
barcode repository containing fastq/gz files and used for analysis:
EOF

    if [ -d "${input_dir_path}" ]; then
        for dir in ${input_dir_path}/*; do
            if [ -d "\$dir" ]; then
                bname=\$(basename "\$dir")
                if find "\$dir" -maxdepth 1 -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" | grep -q .; then
                    echo "\$bname" >> fastq_content.txt
                fi
            fi
        done
    fi

    cat << 'EOF' >> fastq_content.txt
##########################
barcode repository containing other files than fastq/gz and ignored for analysis:
##########################
list of problematic files:
##########################
empty barcode repositories and ignored for analysis:
EOF

    # 2. Génération de param_file.txt sécurisé (sans évaluation Bash/Nextflow conflictuelle)
    cat << 'EOF' > param_file.txt
######################
#### PARAMS USED #####
######################
data repository: ${params.input}
result repository: ${params.outdir}
database used: ${params.orig_ref}
read minlength: ${params.min_length}
read maxlength: ${params.max_length}
quality filtering: 12
read 5' trimming length: NA
read 3' trimming length: NA
min coverage for consensus generation: NA
variant calling depth: NA
variant calling base quality threeshold: NA
multi-infection cutoff: ${params.mi_cutoff}
variant frequency: ${params.freq_min}
mutation research: ${params.window_pos}
read correction enabled: false
coverage correction: false
EOF
    """
}