#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    VIRiONT_NF - PIPELINE COMPLET (VHB WG/POL & VHD WG/R0)
========================================================================================
*/

// Initialisation propre des paramètres par défaut (évite les warnings Nextflow)
params.fastq_dir        = params.fastq_dir        ?: null
params.outdir           = params.outdir           ?: "results"
params.run_id           = params.run_id           ?: null
params.medaka_model     = params.medaka_model     ?: "r1041_e82_400bps_sup_v4.2.0"
params.clair3_model     = params.clair3_model     ?: "r1041_e82_400bps_sup_v4.2.0"
params.primers_adapters = params.primers_adapters ?: "${projectDir}/ref/porechop_adapters.txt"
params.mutation_tables  = params.mutation_tables  ?: "${projectDir}/mutation_table"
params.min_reads        = params.min_reads        ?: 20
params.max_reads        = params.max_reads        ?: 1000
params.mi_cutoff        = params.mi_cutoff        ?: 30.0
params.hostile_index    = params.hostile_index    ?: "human-t2t-hla"
params.hostile_cache    = params.hostile_cache    ?: "${projectDir}/hostile_cache"

include { GENERATE_RUN_METADATA                          } from './modules/00_init_logs.nf'
include { MERGE_FASTQ                                    } from './modules/01_merge.nf'
include { DEHOST_HOSTILE                                 } from './modules/02_dehost.nf'
include { TRIM_CHOPPER                                   } from './modules/03_chopper.nf'
include { TRIM_PRIMERS                                   } from './modules/03b_Trimming.nf'
include { COMPETITIVE_ALIGN                              } from './modules/04_competitive_align.nf'
include { MEDAKA_CONSENSUS; COLLECT_CONSENSUS            } from './modules/05_consensus.nf'
include { RECOMBINATION_JPHMM; COLLECT_RECOMBINATION_SUMMARIES } from './modules/06_recombination.nf'
include { CLINICAL_FILTER                                } from './modules/07_clinical_filter.nf'
include { ALIGN_GENOTYPE_BAM                             } from './modules/08_canonical_bam.nf'
include { CALL_VCF                                       } from './modules/09_vcf.nf'
include { QC_MOSDEPTH; QC_ANALYSIS                       } from './modules/10_qc_analysis.nf'
include { PREPARE_TREE_REFS; PHYLOGENY; PLOT_TREE        } from './modules/11_phylogeny.nf'
include { COMPUTE_BAM_COVERAGE; PLOT_GLOBAL_COVERAGE     } from './modules/12_coverage.nf'
include { SEARCH_HBV_MUTATIONS; COLLECT_MUTATION_REPORTS; SEARCH_HBV_INDELS; COLLECT_INDEL_REPORTS } from './modules/13_mutation.nf'
include { SV_SPLICING_VHB; COLLECT_SPLICING_REPORTS      } from './modules/14_sv_splicing.nf'
include { FASTFINDER                                     } from './modules/15_fastfinder.nf'

workflow {
    def fastq_path = file(params.fastq_dir)

    def csv_file = fastq_path.listFiles().find { file ->
        file.name.toLowerCase().endsWith('.csv')
    }
    if (!csv_file) {
        error "Aucun fichier .csv valide trouvé dans : ${params.fastq_dir}/"
    }

    def csv_upper  = csv_file.name.toUpperCase()
    def virus_name = (csv_upper.contains("VHD") || csv_upper.contains("HDV")) ? "VHD" : "VHB"
    def tech_name  = "WG"

    if (csv_upper.contains("POL")) {
        tech_name = "POL"
    } else if (csv_upper.contains("R0")) {
        tech_name = "R0"
    }

    // =====================================================================================
    // EXTRACTION STRICTE DE LA DATE À 6 CHIFFRES (ex: 260931) SANS LE NUMÉRO DE RUN
    // =====================================================================================
    def raw_folder   = (fastq_path.name == 'fastq_pass') ? fastq_path.parent.name : fastq_path.name
    def clean_folder = raw_folder.replaceAll(/^ONT_run_/, '')

    def all_dates = ("${clean_folder}_${csv_file.name}" =~ /(?<![0-9])[0-9]{6}(?![0-9])/).findAll()
    def run_date  = all_dates ? all_dates[0] : (clean_folder =~ /[0-9]{6}/ ? (clean_folder =~ /[0-9]{6}/)[0] : "RUN")

    def run_id = params.run_id ?: "${run_date}_${virus_name}_${tech_name}"

    // =====================================================================================
    // BASES DE DONNÉES DE RÉFÉRENCE SELON LA MATRICE BIOLOGIQUE
    // =====================================================================================
    def min_length          = 1000
    def max_length          = 5000
    def ref_db_primary      = null
    def tree_ref            = null
    def canonical_align_ref = null

    def jphmm_ecori_ref = file("${projectDir}/ref/HBV_jphmm_canonical_EcoRI.fasta")
    def primers_file    = file(params.primers_adapters, checkIfExists: true)

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
    log.info "IDENTIFIANT RUN      : ${run_id}"
    log.info "VIRUS                : ${virus_name}"
    log.info "PROTOCOLE            : ${tech_name}"
    log.info "LONGUEURS CHOPPER    : ${min_length} - ${max_length} bp"
    log.info "REF PRIMARY/UNIQUE   : ${ref_db_primary.name}"
    log.info "REF CANONIQUE        : ${canonical_align_ref.name}"
    log.info "ROGNAGE AMORCES      : Porechop_ABI (${primers_file.name})"
    log.info "REF ECORI JPHMM      : ${jphmm_ecori_ref.exists() ? jphmm_ecori_ref.name : 'absente (ref/)'}"
    log.info "MODULE RECOMBINAISON : ${(virus_name == 'VHB' && tech_name == 'WG') ? 'ACTIF (jpHMM sur isolats validés)' : 'IGNORÉ (Non-applicable)'}"
    log.info "MODULE MUTATION      : ${virus_name == 'VHB' ? 'ACTIF (Sortie: 13_MUTATION)' : 'IGNORÉ (Aucune sortie)'}"
    log.info "MODULE SV SPLICING   : ${(virus_name == 'VHB' && tech_name == 'WG') ? 'ACTIF (Sniffles2 sur amplicons WG)' : 'IGNORÉ (Non-applicable)'}"
    log.info "MODULE FASTFINDER    : ACTIF (Sortie: 15_FASTFINDER)"
    log.info "========================================================="

    // 00. Génération des métadonnées du run
    GENERATE_RUN_METADATA(
        fastq_path,
        ref_db_primary.name,
        canonical_align_ref.name,
        min_length,
        max_length,
        virus_name,
        tech_name
    )

    // Parsing de la Sample Sheet
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
    
    // 01. Préparation, filtrage qualité (Chopper) et Amorces (Porechop_ABI)
    MERGE_FASTQ(samples_ch)
    DEHOST_HOSTILE(MERGE_FASTQ.out.merged_fastq)
    TRIM_CHOPPER(DEHOST_HOSTILE.out.dehosted_fastq, min_length, max_length)

    // Étape 03b : Rognage des amorces
    TRIM_PRIMERS(TRIM_CHOPPER.out.trimmed_fastq, primers_file)

    TRIM_PRIMERS.out.fastq
        .filter { sample_id, fq -> fq.exists() && fq.size() > 100 }
        .set { valid_trimmed_ch }

    // 02. Alignement compétitif
    COMPETITIVE_ALIGN(
        valid_trimmed_ch,
        ref_db_primary,
        virus_name,
        tech_name
    )

    // Préparation du canal pour Medaka
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
        .set { medaka_input_ch }

    // 03. Polissage Medaka
    MEDAKA_CONSENSUS(medaka_input_ch, params.medaka_model)

    COLLECT_CONSENSUS(
        MEDAKA_CONSENSUS.out.final_consensus.map { sample_id, geno, fasta -> fasta }.collect().ifEmpty([])
    )

    // 04. Filtrage clinique et validation des isolats (avec tech_name)
    CLINICAL_FILTER(
        MEDAKA_CONSENSUS.out.final_consensus.map { sample_id, geno, fasta -> fasta }.collect().ifEmpty([]) ,
        valid_trimmed_ch.map { id, fq -> fq }.collect().ifEmpty([]),
        ref_db_primary,
        COMPETITIVE_ALIGN.out.real_counts.map { id, tsv -> tsv }.collect().ifEmpty([]),
        virus_name,
        tech_name
    )

    // 04bis. Détection des Recombinaisons (jpHMM UNIQUEMENT sur VHB WG)
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

    // 05. Alignement sur référence canonique
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

    // 06. Variant calling et screening des mutations (Substitutions ET Indels)
    def all_vcfs_ch = Channel.value([])

    if (virus_name == "VHB") {
        CALL_VCF(ALIGN_GENOTYPE_BAM.out.for_vcf, params.clair3_model)
        
        all_vcfs_ch = CALL_VCF.out.vcf_tbi
            .map { sample_id, genotype, vcf, tbi -> vcf }
            .collect()
            .ifEmpty([])

        def vcf_triplet_ch = CALL_VCF.out.vcf_tbi
            .map { sample_id, genotype, vcf, tbi -> tuple(sample_id, genotype, vcf) }

        // 1. Substitutions (script R original sanctuarisé)
        SEARCH_HBV_MUTATIONS(
            vcf_triplet_ch,
            virus_name,
            file(params.mutation_tables),
            canonical_align_ref
        )

        COLLECT_MUTATION_REPORTS(
            SEARCH_HBV_MUTATIONS.out.sample_raw_variants.collect().ifEmpty([])
        )

        // 2. Indels & Frameshifts (script R dédié et étanche)
        SEARCH_HBV_INDELS(
            vcf_triplet_ch,
            virus_name,
            file(params.mutation_tables),
            canonical_align_ref
        )

        COLLECT_INDEL_REPORTS(
            SEARCH_HBV_INDELS.out.sample_raw_indels.collect().ifEmpty([])
        )
    }

    // 07. Contrôle Qualité (QC)
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

    // 08. Analyse Phylogénétique
    PREPARE_TREE_REFS(
        tree_ref,
        CLINICAL_FILTER.out.summary_tsv
    )
    PHYLOGENY(
        CLINICAL_FILTER.out.validated_all_fasta,
        PREPARE_TREE_REFS.out.filtered_refs
    )
    PLOT_TREE(PHYLOGENY.out.tree_file)

    // 09. Profil de couverture
    COMPUTE_BAM_COVERAGE(ALIGN_GENOTYPE_BAM.out.bam_bai)
    PLOT_GLOBAL_COVERAGE(
        COMPUTE_BAM_COVERAGE.out.sample_cov.collect().ifEmpty([]),
        virus_name,
        tech_name
    )

    // 10. Détection des variants structuraux / épissage (VHB WG uniquement)
    if (virus_name == "VHB" && tech_name == "WG") {
        SV_SPLICING_VHB(ALIGN_GENOTYPE_BAM.out.for_vcf)

        COLLECT_SPLICING_REPORTS(
            SV_SPLICING_VHB.out.table.map { sample_id, geno, tsv -> tsv }.collect().ifEmpty([])
        )
    }

    // 11. Export FastFinder / GLIMS
    FASTFINDER(
        CLINICAL_FILTER.out.summary_tsv,
        CLINICAL_FILTER.out.validated_all_fasta,
        run_id
    )
}