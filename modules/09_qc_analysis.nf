nextflow.enable.dsl=2

/*
========================================================================================
    MODULE 09: CONTROL QUALITÉ CENTRALISÉ (09_QC_ANALYSIS)
========================================================================================
*/
process QC_ANALYSIS {
    publishDir "${params.outdir}/09_QC_ANALYSIS", mode: 'copy'

    input:
    path raw_fastqs
    path dehosted_fastqs
    path trimmed_fastqs
    path bams
    path bais
    path vcfs
    path summary_tsv

    output:
    path "QC_Metrics_Summary.tsv", emit: qc_table
    path "*.mosdepth.*"          , optional: true

    script:
    """
    echo "=== 1. Calcul des couvertures via Mosdepth ==="
    for bam in ${bams}; do
        if [ -f "\$bam" ]; then
            prefix=\$(basename "\$bam" .bam)
            mosdepth -n --fast-mode "\${prefix}_cov" "\$bam"
        fi
    done

    echo "=== 2. Génération de la table de synthèse QC ==="
    cat << 'EOF' > generate_qc.R
    # Extraction des longueurs de reads basées uniquement sur les séquences
    get_lengths <- function(fq_path) {
      if (!file.exists(fq_path) || file.size(fq_path) == 0) return(numeric(0))
      cmd <- paste0("zcat -f '", fq_path, "' | awk 'NR%4==2 {print length(\$0)}'")
      lens <- as.numeric(system(cmd, intern = TRUE))
      return(lens[!is.na(lens)])
    }

    # Lecture du résumé de génotypage
    summary_dt <- read.delim("${summary_tsv}", stringsAsFactors = FALSE)

    qc_list <- list()

    # Liste des fichiers bruts uniques
    raw_files <- list.files(".", pattern = ".*_merged\\\\.fastq(\\\\.gz)?", full.names = TRUE)

    for (rf in raw_files) {
      sid <- gsub("_merged\\\\.fastq(\\\\.gz)?", "", basename(rf))
      sid <- gsub("^\\\\./", "", sid)

      # 01_RAW
      lens_raw <- get_lengths(rf)
      if (length(lens_raw) > 0) {
        qc_list[[length(qc_list) + 1]] <- data.frame(
          sample = sid, step = "01_RAW", assignedref = "NONE",
          read_count = length(lens_raw), minlengthread = min(lens_raw),
          maxlengthread = max(lens_raw), meanread_length = round(mean(lens_raw), 1),
          medianread_length = as.character(round(median(lens_raw), 1)),
          pident_blast = "NA", mean_depth_coverage = "NA", clair3_variants = "NA", assigned_reads = "NA",
          stringsAsFactors = FALSE
        )
      }

      # 02_DEHOSTING
      dh_file <- list.files(".", pattern = paste0("^", sid, ".*dehosted.*fastq(\\\\.gz)?"), full.names = TRUE)[1]
      if (!is.na(dh_file) && file.exists(dh_file)) {
        lens_dh <- get_lengths(dh_file)
        if (length(lens_dh) > 0) {
          qc_list[[length(qc_list) + 1]] <- data.frame(
            sample = sid, step = "02_DEHOSTING", assignedref = "NONE",
            read_count = length(lens_dh), minlengthread = min(lens_dh),
            maxlengthread = max(lens_dh), meanread_length = round(mean(lens_dh), 1),
            medianread_length = as.character(round(median(lens_dh), 1)),
            pident_blast = "NA", mean_depth_coverage = "NA", clair3_variants = "NA", assigned_reads = "NA",
            stringsAsFactors = FALSE
          )
        }
      }

      # 03_FILTERED_TRIMMED
      tr_file <- list.files(".", pattern = paste0("^", sid, ".*trimmed.*fastq(\\\\.gz)?"), full.names = TRUE)[1]
      if (!is.na(tr_file) && file.exists(tr_file)) {
        lens_tr <- get_lengths(tr_file)
        if (length(lens_tr) > 0) {
          qc_list[[length(qc_list) + 1]] <- data.frame(
            sample = sid, step = "03_FILTERED_TRIMMED", assignedref = "NONE",
            read_count = length(lens_tr), minlengthread = min(lens_tr),
            maxlengthread = max(lens_tr), meanread_length = round(mean(lens_tr), 1),
            medianread_length = as.character(round(median(lens_tr), 1)),
            pident_blast = "NA", mean_depth_coverage = "NA", clair3_variants = "NA", assigned_reads = "NA",
            stringsAsFactors = FALSE
          )
        }
      }

      # 05_GENOTYPING
      sample_summary <- summary_dt[summary_dt\$sample == sid & summary_dt\$status == "VALIDATED", ]
      
      # Extraction de la couverture moyenne depuis le summary de Mosdepth
      mosdepth_summary <- list.files(".", pattern = paste0("^", sid, ".*\\\\.mosdepth\\\\.summary\\\\.txt\$"), full.names = TRUE)[1]
      mean_cov <- "NA"
      if (!is.na(mosdepth_summary) && file.exists(mosdepth_summary)) {
        ms_data <- read.delim(mosdepth_summary, stringsAsFactors = FALSE)
        total_row <- ms_data[ms_data\$chrom == "total", ]
        if (nrow(total_row) > 0) {
          mean_cov <- round(as.numeric(total_row\$mean[1]), 2)
        }
      }

      # Nombre de variants Clair3
      vcf_file <- list.files(".", pattern = paste0("^", sid, ".*\\\\.vcf(\\\\.gz)?\$"), full.names = TRUE)[1]
      n_var <- "NA"
      if (!is.na(vcf_file) && file.exists(vcf_file)) {
        vcf_cmd <- paste0("bcftools view -H '", vcf_file, "' | wc -l")
        n_var <- trimws(system(vcf_cmd, intern = TRUE))
      }

      if (nrow(sample_summary) > 0) {
        for (i in 1:nrow(sample_summary)) {
          qc_list[[length(qc_list) + 1]] <- data.frame(
            sample = sid,
            step = "05_GENOTYPING",
            assignedref = as.character(sample_summary\$genotype[i]),
            read_count = as.numeric(sample_summary\$total_reads[i]),
            minlengthread = 1000,
            maxlengthread = 5000,
            meanread_length = "NA",
            medianread_length = "NA",
            pident_blast = round(as.numeric(sample_summary\$pident[i]), 3),
            mean_depth_coverage = mean_cov,
            clair3_variants = n_var,
            assigned_reads = as.numeric(sample_summary\$total_reads[i]),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    final_df <- do.call(rbind, qc_list)
    write.table(final_df, "QC_Metrics_Summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    EOF

    Rscript generate_qc.R
    """
}