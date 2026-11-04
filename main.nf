#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    VIRiONT_NF - PIPELINE COMPLET (VHB WG/POL & VHD WG/R0)
========================================================================================
*/

params.fastq_dir            = params.fastq_dir            ?: null
params.outdir               = params.outdir               ?: "results"
params.run_id               = params.run_id               ?: null
params.consensus_min_depth  = params.consensus_min_depth  ?: 10
params.iupac_min_freq       = params.iupac_min_freq != null ? params.iupac_min_freq : 0.08
params.clair3_model         = params.clair3_model         ?: "r1041_e82_400bps_sup_v4.2.0"
params.primers_adapters     = params.primers_adapters     ?: "${projectDir}/ref/porechop_adapters.txt"
params.mutation_tables      = params.mutation_tables      ?: "${projectDir}/mutation_table"
params.min_reads            = params.min_reads            ?: 20
params.max_reads            = params.max_reads            ?: 1000
params.mi_cutoff            = params.mi_cutoff            ?: 30.0
params.hostile_index        = params.hostile_index        ?: "human-t2t-hla"
params.hostile_cache        = params.hostile_cache        ?: "${projectDir}/hostile_cache"

// --- PARAMÈTRES MUTATIONS & GLIMS ---
params.freq_min_filtered    = params.freq_min_filtered    ?: 1.0
params.glims_min_minor_freq = params.glims_min_minor_freq ?: 8.0

// --- PARAMÈTRES FLAIR ---
params.min_flair_pct        = params.min_flair_pct        ?: 8.0
params.min_flair_size       = params.min_flair_size       ?: 150

include { GENERATE_RUN_METADATA                                     } from './modules/00_init_logs.nf'
include { MERGE_FASTQ                                               } from './modules/01_merge.nf'
include { DEHOST_HOSTILE                                            } from './modules/02_dehost.nf'
include { TRIM_CHOPPER                                              } from './modules/03_chopper.nf'
include { TRIM_PRIMERS                                              } from './modules/03b_Trimming.nf'
include { COMPETITIVE_ALIGN                                         } from './modules/04_competitive_align.nf'
include { ALIGN_SUBTYPE; GENERATE_CONSENSUS; COLLECT_CONSENSUS      } from './modules/05_consensus.nf'
include { RECOMBINATION_JPHMM; COLLECT_RECOMBINATION_SUMMARIES      } from './modules/06_recombination.nf'
include { CLINICAL_FILTER                                          } from './modules/07_clinical_filter.nf'
include { ALIGN_GENOTYPE_BAM                                       } from './modules/08_canonical_bam.nf'
include { CALL_VCF; RAW_PILEUP_COUNTS                                } from './modules/09_vcf.nf'
include { QC_MOSDEPTH; QC_ANALYSIS                                  } from './modules/10_qc_analysis.nf'
include { PREPARE_TREE_REFS; PHYLOGENY; PLOT_TREE                   } from './modules/11_phylogeny.nf'
include { COMPUTE_BAM_COVERAGE; PLOT_GLOBAL_COVERAGE                } from './modules/12_coverage.nf'
include { SEARCH_HBV_MUTATIONS; COLLECT_MUTATION_REPORTS; SEARCH_HBV_INDELS; COLLECT_INDEL_REPORTS } from './modules/13_mutation.nf'
include { FLAIR; MERGE_SUMMARIES                                    } from './modules/14_flair.nf'
include { FASTFINDER                                                } from './modules/15_fastfinder.nf'

workflow {
    if (!params.fastq_dir) {
        error "ERREUR : Le paramètre obligatoire '--fastq_dir' est manquant."
    }

    def fastq_path = file(params.fastq_dir)
    if (!fastq_path.exists()) {
        error "ERREUR : Le dossier spécifié n'existe pas : ${params.fastq_dir}"
    }

    def csv_file = fastq_path.listFiles()?.find { file ->
        file.name.toLowerCase().endsWith('.csv')
    }
    if (!csv_file) {
        error "ERREUR : Aucun fichier .csv trouvé dans : ${params.fastq_dir}/"
    }

    def csv_upper  = csv_file.name.toUpperCase()
    def virus_name = (csv_upper.contains("VHD") || csv_upper.contains("HDV")) ? "VHD" : "VHB"
    def tech_name  = "WG"

    if (csv_upper.contains("POL")) {
        tech_name = "POL"
    } else if (csv_upper.contains("R0")) {
        tech_name = "R0"
    }

    def raw_folder   = (fastq_path.name == 'fastq_pass') ? fastq_path.parent.name : fastq_path.name
    def clean_folder = raw_folder.replaceAll(/^ONT_run_/, '')

    def all_dates = ("${clean_folder}_${csv_file.name}" =~ /(?<![0-9])[0-9]{6}(?![0-9])/).findAll()
    def run_date  = all_dates ? all_dates[0] : "RUN"

    def run_id = params.run_id ?: "${run_date}_${virus_name}_${tech_name}"

    def min_length          = 1000
    def max_length          = 5000
    def ref_db_primary      = null
    def tree_ref            = null
    def canonical_align_ref = null

    def jphmm_ecori_ref     = file("${projectDir}/ref/HBV_jphmm_canonical_EcoRI.fasta")
    def primers_file        = file(params.primers_adapters, checkIfExists: true)
    def splice_catalog_file = file("${projectDir}/ref/hbv_splice_catalog.tsv")

    if (virus_name == "VHD") {
        if (tech_name == "R0") {
            min_length          = 300
            max_length          = 800
            tree_ref            = file("${projectDir}/ref/HDV_subtype_R0.fasta", checkIfExists: true)
            ref_db_primary      = tree_ref
            canonical_align_ref = tree_ref
        } else {
            min_length          = 1000
            max_length          = 2000
            tree_ref            = file("${projectDir}/ref/HDV_subtype_WG.fasta", checkIfExists: true)
            ref_db_primary      = tree_ref
            canonical_align_ref = tree_ref
        }
    } else if (tech_name == "POL") {
        min_length          = 800
        max_length          = 5000
        tree_ref            = file("${projectDir}/ref/HBV_subtype_POL.fasta", checkIfExists: true)
        ref_db_primary      = tree_ref
        canonical_align_ref = file("${projectDir}/ref/HBV_genotype_wPrimer.fasta", checkIfExists: true)
    } else {
        min_length          = 1000
        max_length          = 5000
        ref_db_primary      = file("${projectDir}/ref/HBV_subtype_WG.fasta", checkIfExists: true)
        tree_ref            = ref_db_primary
        canonical_align_ref = file("${projectDir}/ref/HBV_genotype_wPrimer.fasta", checkIfExists: true)
    }

    log.info "========================================================="
    log.info "SAMPLE SHEET RETENUE: ${csv_file.name}"
    log.info "IDENTIFIANT RUN     : ${run_id}"
    log.info "VIRUS               : ${virus_name}"
    log.info "PROTOCOLE           : ${tech_name}"
    log.info "LONGUEURS CHOPPER   : ${min_length} - ${max_length} bp"
    log.info "REF PRIMARY/UNIQUE  : ${ref_db_primary.name}"
    log.info "REF CANONIQUE       : ${canonical_align_ref.name}"
    log.info "ROGNAGE AMORCES     : Porechop_ABI (${primers_file.name})"
    log.info "CONSENSUS           : samtools consensus --ambig (freq>=${params.iupac_min_freq}, prof>=${params.consensus_min_depth})"
    log.info "REF ECORI JPHMM     : ${jphmm_ecori_ref.exists() ? jphmm_ecori_ref.name : 'absente'}"
    log.info "CATALOGUE ÉPISSAGE  : ${splice_catalog_file.exists() ? splice_catalog_file.name : 'absent'}"
    log.info "MODULE RECOMBINAISON: ${(virus_name == 'VHB' && tech_name == 'WG') ? 'ACTIF' : 'IGNORÉ'}"
    log.info "MODULE MUTATION     : ${virus_name == 'VHB' ? 'ACTIF' : 'IGNORÉ'}"
    log.info "SEUIL GLIMS MINO    : >= ${params.glims_min_minor_freq}%"
    log.info "MODULE ÉPISSAGE     : ${(virus_name == 'VHB' && tech_name == 'WG') ? 'ACTIF (FLAIR sur amplicons WG)' : 'IGNORÉ'}"
    log.info "MODULE FASTFINDER   : ACTIF"
    log.info "========================================================="

    GENERATE_RUN_METADATA(
        fastq_path,
        ref_db_primary.name,
        canonical_align_ref.name,
        min_length,
        max_length,
        virus_name,
        tech_name
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
            
            sample    = sample.toString().trim().replaceAll('\\s+', '_')
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

    MERGE_FASTQ(samples_ch)
    DEHOST_HOSTILE(MERGE_FASTQ.out.merged_fastq)
    TRIM_CHOPPER(DEHOST_HOSTILE.out.dehosted_fastq, min_length, max_length)

    TRIM_PRIMERS(TRIM_CHOPPER.out.trimmed_fastq, primers_file)

    TRIM_PRIMERS.out.fastq
        .filter { sample_id, fq -> fq.exists() && fq.size() > 100 }
        .set { valid_trimmed_ch }

    COMPETITIVE_ALIGN(
        valid_trimmed_ch,
        ref_db_primary,
        virus_name,
        tech_name
    )

    COMPETITIVE_ALIGN.out.partitioned_reads
        .flatMap { sample_id, fq_list, ref_list ->
            def fqs = fq_list instanceof List ? fq_list : [fq_list]
            def refs = ref_list instanceof List ? ref_list : [ref_list]
            def out = []
            def pattern = java.util.regex.Pattern.compile("^" + java.util.regex.Pattern.quote(sample_id) + "_(.+)\\.fastq\\.gz\$")
            fqs.each { fq ->
                def m = pattern.matcher(fq.name)
                if (m.find()) {
                    def geno = m.group(1)
                    def matching_ref = refs.find { it.name == "${sample_id}_${geno}.ref.fasta" }
                    if (matching_ref) {
                        out << tuple(sample_id, geno, fq, matching_ref)
                    }
                }
            }
            return out
        }
        .set { consensus_input_ch }

    // 05_CONSENSUS : Alignement sous-type + consensus IUPAC
    ALIGN_SUBTYPE(consensus_input_ch)
    GENERATE_CONSENSUS(ALIGN_SUBTYPE.out.bam_ref)

    COLLECT_CONSENSUS(
        GENERATE_CONSENSUS.out.final_consensus.map { sample_id, geno, fasta -> fasta }.collect().ifEmpty([])
    )

    CLINICAL_FILTER(
        GENERATE_CONSENSUS.out.final_consensus.map { sample_id, geno, fasta -> fasta }.collect().ifEmpty([]),
        valid_trimmed_ch.map { id, fq -> fq }.collect().ifEmpty([]),
        ref_db_primary,
        COMPETITIVE_ALIGN.out.real_counts.map { id, tsv -> tsv }.collect().ifEmpty([]),
        virus_name,
        tech_name
    )

    if (virus_name == "VHB" && tech_name == "WG") {
        CLINICAL_FILTER.out.validated_all_fasta
            .splitFasta(record: [id: true, seqString: true])
            .map { record ->
                def clean_id = record.id.tokenize(' ')[0]
                def fasta_file = file("${workDir}/tmp_validated_${clean_id}.fasta")
                fasta_file.text = ">${clean_id}\n${record.seqString}\n"
                tuple(clean_id, fasta_file)
            }
            .set { validated_jphmm_input_ch }

        RECOMBINATION_JPHMM(validated_jphmm_input_ch, jphmm_ecori_ref)

        COLLECT_RECOMBINATION_SUMMARIES(
            RECOMBINATION_JPHMM.out.summary.map { id, file -> file }.collect().ifEmpty([])
        )
    }

    CLINICAL_FILTER.out.genotyped_fastas
        .flatten()
        .map { f -> tuple(f.name, f) }
        .set { genotype_fasta_by_name }

    CLINICAL_FILTER.out.genotype_manifest
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.filename, row.sample_id, row.genotype, row.protocol) }
        .combine(genotype_fasta_by_name, by: 0)
        .map { filename, sample_id, genotype, protocol, fasta -> tuple(sample_id, genotype, protocol, fasta) }
        .set { genotype_ref_ch }

    genotype_ref_ch
        .combine(valid_trimmed_ch, by: 0)
        .map { sample_id, genotype, protocol, fasta, trimmed_fastq -> tuple(sample_id, genotype, protocol, trimmed_fastq, fasta) }
        .set { genotype_align_input_ch }

    ALIGN_GENOTYPE_BAM(
        genotype_align_input_ch,
        canonical_align_ref,
        virus_name
    )

    // Initialisation des canaux VCF et GLIMS
    def all_vcfs_ch = Channel.value([])
    def glims_alerts_file_ch = Channel.value(file("${projectDir}/assets/NO_FILE"))
    def glims_indels_file_ch = Channel.value(file("${projectDir}/assets/NO_FILE"))

    if (virus_name == "VHB") {
        CALL_VCF(ALIGN_GENOTYPE_BAM.out.for_vcf, params.clair3_model)
        RAW_PILEUP_COUNTS(ALIGN_GENOTYPE_BAM.out.for_vcf)

        all_vcfs_ch = CALL_VCF.out
            .map { sample_id, genotype, vcf -> vcf }
            .collect()
            .ifEmpty([])

        def raw_triplet_ch = RAW_PILEUP_COUNTS.out
            .map { sample_id, genotype, vcf -> tuple(sample_id, genotype, vcf) }

        // Substitutions
        SEARCH_HBV_MUTATIONS(
            raw_triplet_ch,
            virus_name,
            file(params.mutation_tables),
            canonical_align_ref
        )

        COLLECT_MUTATION_REPORTS(
            SEARCH_HBV_MUTATIONS.out.sample_raw_variants.collect().ifEmpty([]),
            SEARCH_HBV_MUTATIONS.out.sample_filtered_variants.collect().ifEmpty([]),
            SEARCH_HBV_MUTATIONS.out.sample_glims_variants_alerts.collect().ifEmpty([])
        )

        // Indels & Frameshifts
        SEARCH_HBV_INDELS(
            CALL_VCF.out,
            virus_name,
            file(params.mutation_tables),
            canonical_align_ref
        )

        COLLECT_INDEL_REPORTS(
            SEARCH_HBV_INDELS.out.sample_raw_indels.collect().ifEmpty([]),
            SEARCH_HBV_INDELS.out.sample_filtered_indels.collect().ifEmpty([]),
            SEARCH_HBV_INDELS.out.sample_glims_indels_alerts.collect().ifEmpty([])
        )

        glims_alerts_file_ch = COLLECT_MUTATION_REPORTS.out.glims_variants_alerts_csv
        glims_indels_file_ch = COLLECT_INDEL_REPORTS.out.glims_indels_alerts_csv
    }

    MERGE_FASTQ.out.merged_fastq.map { id, fq -> fq }.collect().ifEmpty([]).set { all_raw }
    DEHOST_HOSTILE.out.dehosted_fastq.map { id, fq -> fq }.collect().ifEmpty([]).set { all_dehosted }
    TRIM_PRIMERS.out.fastq.map { id, fq -> fq }.collect().ifEmpty([]).set { all_trimmed }

    ALIGN_GENOTYPE_BAM.out.bam_bai.map { sample_id, geno, bam, bai -> bam }.collect().ifEmpty([]).set { all_bams }
    ALIGN_GENOTYPE_BAM.out.bam_bai.map { sample_id, geno, bam, bai -> bai }.collect().ifEmpty([]).set { all_bais }

    QC_MOSDEPTH(
        all_bams,
        all_bais
    )

    QC_ANALYSIS(
        all_raw,
        all_dehosted,
        all_trimmed,
        all_bams,
        all_bais,
        all_vcfs_ch,
        CLINICAL_FILTER.out.summary_tsv,
        CLINICAL_FILTER.out.validated_all_fasta,
        QC_MOSDEPTH.out.cov_files.collect().ifEmpty([]),
        min_length,
        max_length,
        virus_name,
        tech_name
    )

    PREPARE_TREE_REFS(
        tree_ref,
        CLINICAL_FILTER.out.summary_tsv
    )
    PHYLOGENY(
        CLINICAL_FILTER.out.validated_all_fasta,
        PREPARE_TREE_REFS.out.filtered_refs
    )
    PLOT_TREE(PHYLOGENY.out.tree_file)

    COMPUTE_BAM_COVERAGE(ALIGN_GENOTYPE_BAM.out.cov_bam_bai)
    PLOT_GLOBAL_COVERAGE(
        COMPUTE_BAM_COVERAGE.out.sample_cov.collect().ifEmpty([]),
        virus_name,
        tech_name
    )

    // 14. Détection fine des transcrits épissés via FLAIR (VHB WG uniquement)
    if (virus_name == "VHB" && tech_name == "WG") {
        FLAIR(
            genotype_align_input_ch,
            splice_catalog_file
        )

        MERGE_SUMMARIES(
            FLAIR.out.summary_csv.collect().ifEmpty([])
        )
    }

    // 15. Export FastFinder / GLIMS (5 arguments : summary, consensus, substitutions, indels, run_id)
    FASTFINDER(
        CLINICAL_FILTER.out.summary_tsv,
        CLINICAL_FILTER.out.validated_all_fasta,
        glims_alerts_file_ch,
        glims_indels_file_ch,
        run_id
    )
}