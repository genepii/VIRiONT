nextflow.enable.dsl=2
/*
========================================================================================
    MODULE 08: ALIGNEMENT CONTRE RÉFÉRENCE CANONIQUE DU GÉNOTYPE (08_CANONICAL_BAM)
========================================================================================
*/
process ALIGN_GENOTYPE_BAM {
    tag "${sample_id}_${genotype}"
    publishDir { "${params.outdir}/08_CANONICAL_BAM/${sample_id}/${genotype}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(genotype), val(protocol), path(trimmed_fastq), path(consensus_fasta)
    path canonical_ref
    val virus_name

    output:
    // 1. Pour QC Mosdepth & QC Analysis (map-ont STRICTEMENT INTACT)
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), emit: bam_bai

    // 2. Pour Clair3 / VCF Calling & Mutations (map-ont STRICTEMENT INTACT)
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_canonical_ref.fasta"), path("${sample_id}_${genotype}.sorted.bam"), path("${sample_id}_${genotype}.sorted.bam.bai"), emit: for_vcf

    // 3. Pour Module 12 Couverture (splice)
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}.coverage.sorted.bam"), path("${sample_id}_${genotype}.coverage.sorted.bam.bai"), emit: cov_bam_bai

    // 4. Pour Sniffles2 / Module 14 UNIQUEMENT (splice)
    tuple val(sample_id), val(genotype), path("${sample_id}_${genotype}_canonical_ref.fasta"), path("${sample_id}_${genotype}.coverage.sorted.bam"), path("${sample_id}_${genotype}.coverage.sorted.bam.bai"), emit: for_splicing

    script:
    """
    samtools faidx "${canonical_ref}"
    if [ "${virus_name}" = "VHB" ]; then
        # Extraction de la PREMIÈRE lettre A-J trouvée dans le génotype (identique
        # à la logique R de 13_search_mutation.R : gsub non-lettres + substr(1,1)).
        clean_gt=\$(echo "${genotype}" | grep -oE '[A-Ja-j]' | head -n1 | tr 'a-z' 'A-Z')
        if [ -z "\$clean_gt" ]; then
            echo "ERREUR : Impossible d'extraire une lettre de génotype A-J depuis '${genotype}' pour ${sample_id}." >&2
            exit 1
        fi
        clean_gt="GT\${clean_gt}"
        echo "=== Alignement VHB minimap2 contre référence canonique \${clean_gt} pour ${sample_id} (génotype brut: ${genotype}) ==="
        if ! samtools faidx "${canonical_ref}" "\${clean_gt}" > "${sample_id}_${genotype}_canonical_ref.fasta" 2>/dev/null; then
            echo "ERREUR : Contig '\${clean_gt}' introuvable dans ${canonical_ref} pour ${sample_id}/${genotype}. Pas de fallback silencieux — vérifier les headers de la référence et la valeur du génotype." >&2
            exit 1
        fi
    else
        echo "=== Alignement VHD minimap2 contre référence ${genotype} pour ${sample_id} ==="
        if ! samtools faidx "${canonical_ref}" "${genotype}" > "${sample_id}_${genotype}_canonical_ref.fasta" 2>/dev/null; then
            echo "WARNING : Contig '${genotype}' introuvable, utilisation du premier contig disponible de ${canonical_ref} pour ${sample_id}." >&2
            samtools faidx "${canonical_ref}" \$(head -n 1 "${canonical_ref}.fai" | cut -f1) > "${sample_id}_${genotype}_canonical_ref.fasta"
        fi
    fi
    samtools faidx "${sample_id}_${genotype}_canonical_ref.fasta"

    # --- 1. BAM CANONIQUE (Clair3, QC Mosdepth, Mutations) - map-ont STRICTEMENT INTACT ---
    echo "=== 1/2 Alignement standard map-ont (Clair3 / Mutations / QC) ==="
    minimap2 -ax map-ont -t ${task.cpus} "${sample_id}_${genotype}_canonical_ref.fasta" "${trimmed_fastq}" | \
        samtools view -bS -F 4 - | \
        samtools sort -o "${sample_id}_${genotype}.sorted.bam" -
    samtools index "${sample_id}_${genotype}.sorted.bam"

    # --- 2. BAM COUVERTURE & SNIFFLES2 DÉDIÉ (Module 12 & Module 14) ---
    echo "=== 2/2 Alignement spliced (-ax splice) pour couverture et Sniffles2 (${virus_name} ${protocol}) ==="
    minimap2 -ax splice -t ${task.cpus} "${sample_id}_${genotype}_canonical_ref.fasta" "${trimmed_fastq}" | \
        samtools view -bS -F 4 - | \
        samtools sort -o "${sample_id}_${genotype}.coverage.sorted.bam" -
    samtools index "${sample_id}_${genotype}.coverage.sorted.bam"
    """
}