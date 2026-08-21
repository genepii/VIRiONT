#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- ifelse(length(args) >= 1, args[1], "validated_consensus_all.fasta")
if (!file.exists(fasta_file)) {
  cat("Aucun fichier consensus validé trouvé pour la matrice.\n")
  quit(save = "no", status = 0)
}

# Lecture multi-FASTA
lines <- readLines(fasta_file)
seqs <- list()
current_hdr <- ""
current_seq <- ""
for (line in lines) {
  if (grepl("^>", line)) {
    if (current_hdr != "") seqs[[current_hdr]] <- toupper(current_seq)
    current_hdr <- sub("^>", "", line)
    current_hdr <- sub("\\s.*$", "", current_hdr)
    current_seq <- ""
  } else {
    current_seq <- paste0(current_seq, gsub("\\s+", "", line))
  }
}
if (current_hdr != "") seqs[[current_hdr]] <- toupper(current_seq)

samples <- sort(names(seqs))
n <- length(samples)
if (n < 2) {
  cat("Moins de 2 échantillons validés : matrice pairwise non générée.\n")
  quit(save = "no", status = 0)
}

tmpdir <- tempdir()
write_fasta <- function(id, seq, path) writeLines(c(paste0(">", id), seq), path)

# Calcule le nombre de mismatches RÉELS entre 2 séquences via alignement local
# blastn (gère nativement indels/décalages), en sommant tous les HSPs (au cas où
# un grand gap interne — consensus partiel — sépare l'alignement en plusieurs blocs,

blast_mismatches <- function(seq1, seq2) {
  q_path <- file.path(tmpdir, "q.fasta")
  s_path <- file.path(tmpdir, "s.fasta")
  write_fasta("q", seq1, q_path)
  write_fasta("s", seq2, s_path)
  cmd <- sprintf(
    "blastn -query %s -subject %s -outfmt '6 pident length mismatch gapopen qstart qend' -evalue 1e-10",
    shQuote(q_path), shQuote(s_path)
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  if (length(out) == 0) return(list(mismatches = NA, coverage_pct = 0))

  total_mismatch <- 0
  intervals <- list()
  for (line in out) {
    parts <- strsplit(line, "\t")[[1]]
    mismatch <- as.numeric(parts[3])
    qstart   <- as.numeric(parts[5])
    qend     <- as.numeric(parts[6])
    if (qstart > qend) { tmp <- qstart; qstart <- qend; qend <- tmp }
    intervals[[length(intervals) + 1]] <- c(qstart, qend)
    total_mismatch <- total_mismatch + mismatch
  }

  # Fusion des intervalles pour calculer la couverture réelle (sans double-compte
  # en cas de HSPs chevauchants)
  intervals <- intervals[order(sapply(intervals, `[`, 1))]
  merged <- list(intervals[[1]])
  if (length(intervals) > 1) {
    for (k in 2:length(intervals)) {
      last <- merged[[length(merged)]]
      cur <- intervals[[k]]
      if (cur[1] <= last[2] + 1) {
        merged[[length(merged)]] <- c(last[1], max(last[2], cur[2]))
      } else {
        merged[[length(merged) + 1]] <- cur
      }
    }
  }
  covered <- sum(sapply(merged, function(x) x[2] - x[1] + 1))
  min_len <- min(nchar(seq1), nchar(seq2))
  coverage_pct <- if (min_len > 0) round((covered / min_len) * 100, 1) else 0

  list(mismatches = total_mismatch, coverage_pct = coverage_pct)
}

matrix_rows <- list()
plot_data <- data.frame()
for (i in 1:n) {
  for (j in 1:n) {
    s1_id <- samples[i]
    s2_id <- samples[j]
    if (i == j) {
      mismatches <- 0
      coverage_pct <- 100
    } else {
      res <- blast_mismatches(seqs[[s1_id]], seqs[[s2_id]])
      mismatches <- res$mismatches
      coverage_pct <- res$coverage_pct
    }
    if (j >= i) {
      plot_data <- rbind(plot_data, data.frame(
        comp1 = s1_id, comp2 = s2_id, count = mismatches, stringsAsFactors = FALSE
      ))
    }
    matrix_rows[[length(matrix_rows) + 1]] <- data.frame(
      comp1 = s1_id, comp2 = s2_id, count = mismatches,
      coverage_pct = coverage_pct, stringsAsFactors = FALSE
    )
  }
}

write.table(do.call(rbind, matrix_rows), "matrix_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

plot_data$comp1 <- factor(plot_data$comp1, levels = samples)
plot_data$comp2 <- factor(plot_data$comp2, levels = rev(samples))

p <- ggplot(plot_data, aes(x = comp1, y = comp2, fill = count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(count), "NA", count)), color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#2b83ba", "#abdda4", "#fdae61", "#d7191c"),
    trans = "sqrt",
    na.value = "grey80",
    breaks = c(0, 10, 50, 100, 300, 700),
    name = "Mismatches"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("matrix_comp.pdf", plot = p, width = 9, height = 9)
cat("Matrice pairwise (alignement blastn) générée avec succès.\n")