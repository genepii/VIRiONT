#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    VIRiONT_NF - PIPELINE COMPLET (16 MODULES DSL2)
========================================================================================
*/

// Importation de tous les modules
include { MERGE_FASTQ            } from './modules/01_merge.nf'
include { DEHOST_HOSTILE         } from './modules/02_dehost.nf'
include { TRIM_CHOPPER           } from './modules/03_chopper.nf'
include { AMPLICON_SORTER        } from './modules/04_amplicon_sorter.nf'
include { BLAST_GENOTYPE         } from './modules/05_blast_genotype.nf'
include { GLOBAL_REPORT          } from './modules/06_global_report.nf'
include { POLISH_MEDAKA          } from './modules/07_polish_medaka.nf'
include { VARIANT_CALLING_CLAIR3 } from './modules/08_variant_calling.nf'
include { CONCAT_CONSENSUS       } from './modules/09_concat_consensus.nf'
include { EXTRACT_SAMPLE_QC      } from './modules/10_extract_qc.nf'
include { COMPUTE_QC_METRICS     } from './modules/11_compute_metrics.nf'
include { ALIGN_MAFFT            } from './modules/12_align_mafft.nf'
include { PHYLOGENY_IQTREE       } from './modules/13_phylogeny_iqtree.nf'
include { COMPUTE_BAM_COVERAGE   } from './modules/14_bam_coverage.nf'
include { PLOT_GLOBAL_COVERAGE   } from './modules/15_plot_coverage.nf'
include { SEARCH_HBV_MUTATIONS   } from './modules/16_mutation_screening.nf'

process PREPARE_REF {
    tag "${orig_ref.name}"
    publishDir "${params.outdir}/00_SUPDATA", mode: 'copy'

    input:
    path orig_ref

    output:
    path "DB", emit: db_dir
    path "REFSEQ/${orig_ref.name}", emit: ref_seq

    script:
    """
    mkdir -p REFSEQ DB
    cp ${orig_ref} REFSEQ/
    makeblastdb -in ${orig_ref} -out DB/${orig_ref.baseName} -input_type fasta -dbtype nucl > /dev/null
    """
}

workflow {

    // A. Détection du CSV
    def fastq_path = file(params.fastq_dir)
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

    // B. Détermination du Virus, de la Techno et des seuils de longueur
    def csv_upper = csv_file.name.toUpperCase()
    def virus_name = csv_upper.contains("VHD") ? "VHD" : "VHB"
    def tech_name  = csv_upper.contains("R0") ? "R0" : (csv_upper.contains("POL") ? "POL" : "WG")

    def ref_filename = ""
    def min_length   = 1000
    def max_length   = 5000

    if (virus_name == "VHD" && tech_name == "R0") {
        ref_filename = "HDV_subtype_R0.fasta"
        min_length   = 300
        max_length   = 800
    } else if (virus_name == "VHD") {
        ref_filename = "HDV_subtype.fasta"
        min_length   = 1000
        max_length   = 2000
    } else if (virus_name == "VHB" && (tech_name == "R0" || tech_name == "POL")) {
        ref_filename = "HBV_subtype_R0.fasta"
        min_length   = 800
        max_length   = 5000
    } else {
        ref_filename = "HBV_subtype.fasta"
        min_length   = 1000
        max_length   = 5000
    }

    // C. Préparation de la Référence
    def orig_ref_file = file("${params.ref_dir}/${ref_filename}")
    PREPARE_REF(orig_ref_file)

    // D. Parsing CSV
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
            
            // Passage en chemin absolu complet
            def dir_path = file("${params.fastq_dir}/barcode${num_barcode}")
            
            return tuple("barcode_${formatted_barcode}_${sample}", dir_path)
        }
        .filter { it != null }
        .set { samples_ch }

    // =====================================================================================
    // WORKFLOW PRINCIPAL
    // =====================================================================================

    // 01. Fusion & Dehosting
    MERGE_FASTQ(samples_ch)
    DEHOST_HOSTILE(MERGE_FASTQ.out.merged_fastq)

    // 02. Chopper, Amplicon_Sorter, BLAST & Rapport Global
    TRIM_CHOPPER(DEHOST_HOSTILE.out.dehosted_fastq, min_length, max_length)

    // Sécurisation : filtrage des FASTQ non vides pour débloquer l'instanciation
    TRIM_CHOPPER.out.trimmed_fastq
        .filter { sample_id, fq -> fq.exists() && fq.size() > 100 }
        .set { valid_trimmed_ch }

    AMPLICON_SORTER(valid_trimmed_ch, min_length, max_length)
    
    BLAST_GENOTYPE(
        AMPLICON_SORTER.out.consensus_clusters,
        PREPARE_REF.out.db_dir,
        ref_filename
    )
    GLOBAL_REPORT(BLAST_GENOTYPE.out.validated_tsv.map { it[1] }.collect())

    // 03. Polissage Medaka, Variant Calling Clair3 & Concaténation FASTA
    BLAST_GENOTYPE.out.master_fasta
        .flatMap { sample_id, fasta_files ->
            def list = fasta_files instanceof List ? fasta_files : [fasta_files]
            list.collect { f ->
                def geno = f.name.replaceAll("^${sample_id}_", "").replaceAll("_master_consensus\\.fasta\$", "")
                return tuple(sample_id, geno, f)
            }
        }
        .combine(TRIM_CHOPPER.out.trimmed_fastq, by: 0)
        .map { sample_id, geno, master_fasta, trimmed_fastq ->
            return tuple(sample_id, geno, "cluster", master_fasta, trimmed_fastq)
        }
        .set { medaka_input_ch }

    POLISH_MEDAKA(medaka_input_ch, PREPARE_REF.out.ref_seq)
    VARIANT_CALLING_CLAIR3(POLISH_MEDAKA.out.bam_tuple)
    CONCAT_CONSENSUS(POLISH_MEDAKA.out.final_consensus.map { it[1] }.collect())

    // 04. Métrologie QC
    MERGE_FASTQ.out.merged_fastq
        .join(DEHOST_HOSTILE.out.dehosted_fastq)
        .join(TRIM_CHOPPER.out.trimmed_fastq)
        .join(BLAST_GENOTYPE.out.validated_tsv)
        .join(POLISH_MEDAKA.out.bam_tuple.map { sample_id, geno, precons, bam, bai -> tuple(sample_id, bam, bai) })
        .join(VARIANT_CALLING_CLAIR3.out.vcf.map { sample_id, geno, vcf -> tuple(sample_id, vcf) })
        .map { sample_id, raw_fq, dehost_fq, trim_fq, val_tsv, bam, bai, vcf ->
            return tuple(sample_id, raw_fq, dehost_fq, trim_fq, val_tsv, bam, bai, vcf)
        }
        .set { qc_input_ch }

    EXTRACT_SAMPLE_QC(qc_input_ch)
    COMPUTE_QC_METRICS(
        EXTRACT_SAMPLE_QC.out.raw_csv.collect(),
        EXTRACT_SAMPLE_QC.out.dehost_csv.collect(),
        EXTRACT_SAMPLE_QC.out.trimm_csv.collect(),
        EXTRACT_SAMPLE_QC.out.geno_csv.collect()
    )

    // 05. Phylogénie (MAFFT & IQ-TREE)
    ALIGN_MAFFT(
        CONCAT_CONSENSUS.out.all_consensus, 
        PREPARE_REF.out.ref_seq,
        BLAST_GENOTYPE.out.validated_tsv.map { it[1] }.collect()
    )
    PHYLOGENY_IQTREE(ALIGN_MAFFT.out.aligned_fasta)

    // 06. Profils de Couverture BAM
    COMPUTE_BAM_COVERAGE(POLISH_MEDAKA.out.bam_tuple)
    PLOT_GLOBAL_COVERAGE(COMPUTE_BAM_COVERAGE.out.sample_cov.collect())

    // 07. Screening des Mutations Cliniques
    def mutation_tables_path = file(params.mutation_tables_dir)
    SEARCH_HBV_MUTATIONS(
        VARIANT_CALLING_CLAIR3.out.vcf,
        virus_name,
        mutation_tables_path
    )
}