#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
========================================================================================
    VIRiONT_NF - PIPELINE COMPLET (FIL CONDUCTEUR CLINIQUE UNIFIÉ)
========================================================================================
*/
// Importation des modules
include { GENERATE_RUN_METADATA                          } from './modules/00_init_logs.nf'
include { MERGE_FASTQ                                    } from './modules/01_merge.nf'
include { DEHOST_HOSTILE                                 } from './modules/02_dehost.nf'
include { TRIM_CHOPPER                                   } from './modules/03_chopper.nf'
include { AMPLICON_SORTER                                } from './modules/04_amplicon_sorter.nf'
include { MEDAKA_CONSENSUS; COLLECT_CONSENSUS            } from './modules/05_consensus.nf'
include { ALIGN_BAM; COUNT_REAL_READS                    } from './modules/06_bam.nf'
include { CALL_VCF                                       } from './modules/07_vcf.nf'
include { GENOTYPING                                     } from './modules/08_genotyping.nf'
include { QC_ANALYSIS                                    } from './modules/09_qc_analysis.nf'
include { PHYLOGENY                                      } from './modules/10_phylogeny.nf'
include { COMPUTE_BAM_COVERAGE; PLOT_GLOBAL_COVERAGE     } from './modules/11_coverage.nf'
include { SEARCH_HBV_MUTATIONS; COLLECT_MUTATION_REPORTS  } from './modules/12_mutation.nf'
workflow {
    def fastq_path = file(params.fastq_dir)
    // =====================================================================================
    // A. Détection du CSV & Détermination des Paramètres Cliniques
    // =====================================================================================
    def csv_file = fastq_path.listFiles().find { file ->
        file.name.endsWith('.csv') && (
            file.name.contains("Sample_Sheet") ||
            file.name.contains("VHB") ||
            file.name.contains("VHD")
        )
    }
    if (!csv_file) {
        error "❌ Aucun fichier CSV valide trouvé dans : ${params.fastq_dir}/"
    }
    def csv_upper  = csv_file.name.toUpperCase()
    def virus_name = csv_upper.contains("VHD") ? "VHD" : "VHB"
    def tech_name  = csv_upper.contains("R0") ? "R0" : (csv_upper.contains("POL") ? "POL" : "WG")
    def min_length = 1000
    def max_length = 5000
    if (virus_name == "VHD" && tech_name == "R0") {
        min_length = 300
        max_length = 800
    } else if (virus_name == "VHD") {
        min_length = 1000
        max_length = 2000
    } else if (virus_name == "VHB" && tech_name == "POL") {
        min_length = 800
        max_length = 5000
    }
    def target_ref_file = null
    if (virus_name == "VHD") {
        target_ref_file = (tech_name == "R0") ? file("${projectDir}/ref/HDV_subtype_R0.fasta") : file("${projectDir}/ref/HDV_subtype_WG.fasta")
    } else {
        target_ref_file = (tech_name == "POL") ? file("${projectDir}/ref/HBV_subtype_POL.fasta") : file("${projectDir}/ref/HBV_subtype_WG.fasta")
    }
    // =====================================================================================
    // 00. INITIALISATION DES LOGS
    // =====================================================================================
    GENERATE_RUN_METADATA(
        fastq_path,
        target_ref_file,
        min_length,
        max_length
    )
    def csv_text = csv_file.text
    def separator = csv_text.contains(";") ? ';' : ','
    Channel
        .fromPath(csv_file)
        .splitCsv(header: true, sep: separator, quote: '"')
        .map { row ->
            def sample    = row.Sample ?: row.sample ?: row.Sample_ID ?: ""
            def component = row.Component ?: row.component ?: row.Barcode ?: ""
            def plate     = row.Plate_Position ?: row.plate_position ?: ""
            sample    = sample.toString().trim()
            component = component.toString().trim()
            plate     = plate.toString().trim()
            if (!sample || plate.contains("Plate") || sample.contains("Sample")) return null
            def num_match = (component =~ /[0-9]+/)
            if (!num_match) return null
            def num_barcode = num_match[0].toInteger()
            def formatted_barcode = String.format("%02d", num_barcode)
            def dir_path = file("${params.fastq_dir}/barcode${formatted_barcode}")
            if (!dir_path.exists()) {
                dir_path = file("${params.fastq_dir}/barcode${num_barcode}")
            }
            return tuple("barcode_${formatted_barcode}_${sample}", dir_path)
        }
        .filter { it != null }
        .set { samples_ch }
    // =====================================================================================
    // WORKFLOW PRINCIPAL
    // =====================================================================================
    // 01. Préparation Reads
    MERGE_FASTQ(samples_ch)
    DEHOST_HOSTILE(MERGE_FASTQ.out.merged_fastq)
    TRIM_CHOPPER(DEHOST_HOSTILE.out.dehosted_fastq, min_length, max_length)
    TRIM_CHOPPER.out.trimmed_fastq
        .filter { sample_id, fq -> fq.exists() && fq.size() > 100 }
        .set { valid_trimmed_ch }
    // 02. Clustering & Polishing
    AMPLICON_SORTER(valid_trimmed_ch, min_length, max_length)
    valid_trimmed_ch
        .join(AMPLICON_SORTER.out.preconsensus_fasta)
        .set { medaka_input_ch }
    MEDAKA_CONSENSUS(medaka_input_ch, params.medaka_model)
    COLLECT_CONSENSUS(
        MEDAKA_CONSENSUS.out.final_consensus.map { sample_id, fasta -> fasta }.collect()
    )
    // 03. Alignement BAM (sur TOUS les reads trimmed, contre les clusters/génotypes détectés)
    valid_trimmed_ch
        .join(MEDAKA_CONSENSUS.out.final_consensus)
        .set { bam_input_ch }
    ALIGN_BAM(bam_input_ch)

    // 03bis. Comptage exhaustif des reads réels par cluster (remplace le CSV sous-échantillonné
    // d'amplicon_sorter comme source de vérité pour les ratios cliniques)
    COUNT_REAL_READS(ALIGN_BAM.out.bam_bai)

    // 04. Variant Calling Clair3
    CALL_VCF(ALIGN_BAM.out.for_vcf, params.clair3_model)
    // 05. Hub de Génotypage & Création du Fil Conducteur
    COUNT_REAL_READS.out.real_counts
        .map { sample_id, tsv ->
            def target_dir = file("${workDir}/tmp_csvs")
            target_dir.mkdirs()
            def target_file = file("${target_dir}/${sample_id}_results.csv")
            tsv.copyTo(target_file)
            return target_file
        }
        .collect()
        .set { all_sorter_csvs_ch }
    GENOTYPING(
        MEDAKA_CONSENSUS.out.final_consensus.map { id, f -> f }.collect(),
        target_ref_file,
        all_sorter_csvs_ch
    )
    // =====================================================================================
    // 06. CONTRÔLE QUALITÉ CENTRALISÉ (Utilise validated_all_fasta)
    // =====================================================================================
    MERGE_FASTQ.out.merged_fastq.map { id, fq -> fq }.collect().set { all_raw }
    DEHOST_HOSTILE.out.dehosted_fastq.map { id, fq -> fq }.collect().set { all_dehosted }
    TRIM_CHOPPER.out.trimmed_fastq.map { id, fq -> fq }.collect().set { all_trimmed }
    ALIGN_BAM.out.bam_bai.map { id, bam, bai -> bam }.collect().set { all_bams }
    ALIGN_BAM.out.bam_bai.map { id, bam, bai -> bai }.collect().set { all_bais }
    CALL_VCF.out.map { it instanceof List ? it[1] : it }.collect().set { all_vcfs }
    QC_ANALYSIS(
        all_raw,
        all_dehosted,
        all_trimmed,
        all_bams,
        all_bais,
        all_vcfs,
        GENOTYPING.out.summary_tsv,
        GENOTYPING.out.validated_all_fasta,
        min_length,
        max_length
    )
    // =====================================================================================
    // 07. PHYLOGÉNIE (Utilise le même validated_all_fasta)
    // =====================================================================================
    PHYLOGENY(
        GENOTYPING.out.validated_all_fasta,
        target_ref_file
    )
    // =====================================================================================
    // 08. COUVERTURE GÉNOMIQUE
    // =====================================================================================
    COMPUTE_BAM_COVERAGE(ALIGN_BAM.out.bam_bai)
    PLOT_GLOBAL_COVERAGE(
        COMPUTE_BAM_COVERAGE.out.sample_cov.collect(),
        GENOTYPING.out.summary_tsv
    )
    // =====================================================================================
    // 09. SCREENING DES MUTATIONS
    // =====================================================================================
    GENOTYPING.out.summary_tsv
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.status == 'VALIDATED' || row.status == 'validated' }
        .map { row ->
            def sample_id = row.sample
            def geno_raw  = row.genotype ?: "GTD"
            def m = (geno_raw =~ /[A-I]/)
            def gt = m ? "GT${m[0]}" : "GTD"
            return tuple(sample_id, gt)
        }
        .unique()
        .combine(
            CALL_VCF.out.map { item -> tuple(item[0], item[1]) },
            by: 0
        )
        .map { sample_id, gt, vcf_gz ->
            return tuple(sample_id, gt, vcf_gz)
        }
        .set { vcf_genotyped_ch }
    SEARCH_HBV_MUTATIONS(
        vcf_genotyped_ch,
        virus_name,
        file(params.mutation_tables ?: "${projectDir}/mutation_table")
    )
    COLLECT_MUTATION_REPORTS(
        SEARCH_HBV_MUTATIONS.out.sample_variants.collect()
    )
}