#!/usr/bin/env Rscript
# bioxflow/rnaseq — gene-level summarisation and dual-quantification concordance.
#
# 1. Imports Salmon transcript estimates with tximport (keeps average transcript
#    lengths so DESeq2 can correct for differential isoform usage).
# 2. Builds the STAR gene-count matrix, choosing the strand column from the data.
# 3. Compares the two independent quantifications sample by sample and checks that
#    Salmon and STAR agree on library strandedness.

suppressPackageStartupMessages({
  library(tximport)
  library(ggplot2)
})
if (capabilities("cairo")) options(bitmapType = "cairo")

parse_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- 1
  while (i <= length(a)) {
    key <- sub("^--", "", a[i])
    if (!key %in% names(defaults)) stop("unknown argument: ", a[i])
    defaults[[key]] <- a[i + 1]
    i <- i + 2
  }
  defaults
}
opt <- parse_args(list(samples = NULL, tx2gene = NULL, salmon = "salmon", star = "star", threshold = "0.9"))
threshold <- as.numeric(opt$threshold)

samples <- read.delim(opt$samples, stringsAsFactors = FALSE, colClasses = "character")
tx2gene <- read.delim(opt$tx2gene, stringsAsFactors = FALSE, colClasses = "character")
gene_names <- tx2gene$gene_name[match(unique(tx2gene$gene_id), tx2gene$gene_id)]
names(gene_names) <- unique(tx2gene$gene_id)

write_matrix <- function(m, file, digits = 3) {
  out <- data.frame(gene_id = rownames(m), gene_name = gene_names[rownames(m)],
                    round(m, digits), check.names = FALSE)
  write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------------------------------------------------------------- Salmon ----
files <- file.path(opt$salmon, samples$sample, "quant.sf")
names(files) <- samples$sample
missing <- files[!file.exists(files)]
if (length(missing)) stop("missing Salmon results: ", paste(missing, collapse = ", "))

txi_tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreAfterBar = TRUE)
known  <- rownames(txi_tx$counts) %in% tx2gene$transcript_id
if (!any(known)) stop("None of the Salmon transcript IDs are present in the GTF. ",
                      "Check that --transcript_fasta and --gtf come from the same annotation release.")
if (mean(known) < 0.9)
  warning(sprintf("%.1f%% of quantified transcripts are absent from the GTF and were dropped", 100 * mean(!known)))
txi <- summarizeToGene(txi_tx, tx2gene[, c("transcript_id", "gene_id")], ignoreAfterBar = TRUE)
saveRDS(txi, "txi.rds")

write_matrix(txi$counts,    "salmon.gene_counts.tsv")
write_matrix(txi$abundance, "salmon.gene_tpm.tsv")
write_matrix(txi$length,    "salmon.gene_lengths.tsv")
tx_out <- function(m, file) write.table(data.frame(transcript_id = rownames(m), round(m, 3), check.names = FALSE),
                                        file, sep = "\t", quote = FALSE, row.names = FALSE)
tx_out(txi_tx$counts,    "salmon.transcript_counts.tsv")
tx_out(txi_tx$abundance, "salmon.transcript_tpm.tsv")

# Strandedness from Salmon's own mapping-orientation tallies. The tallies are used rather
# than Salmon's "expected_format" call, which stays at its unstranded default when a
# library is too small for Salmon's automatic detection to trigger.
salmon_lib <- lapply(samples$sample, function(s) {
  f <- file.path(opt$salmon, s, "lib_format_counts.json")
  if (!file.exists(f)) return(list(call = NA_character_, frac = NA_real_))
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  n <- function(key) {
    m <- regmatches(txt, regexec(paste0('"', key, '" *: *([0-9]+)'), txt))[[1]]
    if (length(m) > 1) as.numeric(m[2]) else 0
  }
  rev <- sum(vapply(c("ISR", "MSR", "OSR", "SR"), n, 0))
  fwd <- sum(vapply(c("ISF", "MSF", "OSF", "SF"), n, 0))
  if (rev + fwd == 0) return(list(call = "unstranded", frac = NA_real_))
  fr <- rev / (rev + fwd)
  list(call = if (fr >= 0.8) "reverse" else if (fr <= 0.2) "forward" else "unstranded", frac = round(fr, 4))
})
salmon_strand <- vapply(salmon_lib, function(x) x$call, "")

strand <- data.frame(sample = samples$sample, salmon_strandedness = salmon_strand,
                     salmon_fraction_reverse = vapply(salmon_lib, function(x) x$frac, 0),
                     star_strandedness = NA_character_, star_fraction_reverse = NA_real_,
                     agreement = NA_character_, stringsAsFactors = FALSE)

# ------------------------------------------------------------------ STAR ----
star_files <- file.path(opt$star, paste0(samples$sample, ".ReadsPerGene.out.tab"))
have_star  <- dir.exists(opt$star) && all(file.exists(star_files))

if (have_star) {
  star_list <- lapply(seq_along(star_files), function(i) {
    d <- read.delim(star_files[i], header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("gene_id", "unstranded", "forward", "reverse"))
    d <- d[!grepl("^N_", d$gene_id), ]
    fr <- sum(d$reverse) / max(1, sum(d$forward) + sum(d$reverse))
    call <- if (fr >= 0.8) "reverse" else if (fr <= 0.2) "forward" else "unstranded"
    strand$star_strandedness[i]     <<- call
    strand$star_fraction_reverse[i] <<- round(fr, 4)
    setNames(d[[call]], d$gene_id)
  })
  genes <- names(star_list[[1]])
  star  <- do.call(cbind, lapply(star_list, function(v) v[genes]))
  dimnames(star) <- list(genes, samples$sample)
  write_matrix(star, "star.gene_counts.tsv", digits = 0)
  strand$agreement <- ifelse(strand$salmon_strandedness == strand$star_strandedness, "yes", "NO")

  # ----------------------------------------------------- concordance ----
  common <- intersect(rownames(txi$counts), rownames(star))
  conc <- do.call(rbind, lapply(samples$sample, function(s) {
    a <- txi$counts[common, s]; b <- star[common, s]
    keep <- pmax(a, b) >= 10
    if (sum(keep) < 10) keep <- pmax(a, b) > 0
    la <- log2(a[keep] + 1); lb <- log2(b[keep] + 1)
    data.frame(sample = s,
               genes_compared = sum(keep),
               spearman = round(suppressWarnings(cor(a[keep], b[keep], method = "spearman")), 4),
               pearson_log2 = round(suppressWarnings(cor(la, lb)), 4),
               pct_within_2fold = round(100 * mean(abs(la - lb) <= 1), 2),
               salmon_assigned = round(sum(a)),
               star_assigned = sum(b),
               stringsAsFactors = FALSE)
  }))
  conc$strand_agreement <- strand$agreement
  conc$status <- ifelse(!is.na(conc$spearman) & conc$spearman >= threshold & conc$strand_agreement == "yes", "PASS", "REVIEW")
  write.table(conc, "concordance_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  long <- do.call(rbind, lapply(samples$sample, function(s)
    data.frame(sample = s, salmon = txi$counts[common, s], star = star[common, s])))
  long <- long[pmax(long$salmon, long$star) > 0, ]
  lab <- setNames(sprintf("%s  (rho = %.3f)", conc$sample, conc$spearman), conc$sample)
  p <- ggplot(long, aes(log2(star + 1), log2(salmon + 1))) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    geom_point(alpha = 0.25, size = 0.5, colour = "#0e7c86") +
    facet_wrap(~sample, labeller = as_labeller(lab)) +
    labs(x = "STAR gene counts, log2(n + 1)", y = "Salmon gene counts, log2(n + 1)",
         title = "Dual-quantification concordance",
         subtitle = "Each point is a gene; two independent methods should agree along the diagonal") +
    theme_bw(base_size = 11)
  ncol_f <- ceiling(sqrt(nrow(samples)))
  ggsave("concordance_scatter.png", p, width = min(16, 3.2 * ncol_f + 1),
         height = min(16, 3 * ceiling(nrow(samples) / ncol_f) + 1.2), dpi = 120)
  file.copy("concordance_scatter.png", "bioxflow_concordance_scatter_mqc.png")

  header <- c(
    "# id: 'bioxflow_concordance'",
    "# section_name: 'Dual-quantification concordance'",
    sprintf("# description: 'Gene-level agreement between two independent quantification routes: Salmon selective alignment and STAR genome alignment. Samples with Spearman correlation below %s, or where the two methods disagree on library strandedness, are marked REVIEW.'", threshold),
    "# plot_type: 'table'",
    "# pconfig:",
    "#     id: 'bioxflow_concordance_table'",
    "#     title: 'Salmon vs STAR concordance'")
  con <- file("bioxflow_concordance_mqc.tsv", "w")
  writeLines(header, con)
  write.table(conc, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
} else {
  message("STAR counts not provided: concordance check skipped (Salmon-only mode)")
}

write.table(strand, "strandedness.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
header <- c(
  "# id: 'bioxflow_strandedness'",
  "# section_name: 'Library strandedness'",
  "# description: 'Strandedness inferred from the data. Salmon infers it from read-to-transcript mappings; STAR from sense versus antisense gene counts.'",
  "# plot_type: 'table'",
  "# pconfig:",
  "#     id: 'bioxflow_strandedness_table'",
  "#     title: 'Inferred strandedness'")
con <- file("bioxflow_strandedness_mqc.tsv", "w")
writeLines(header, con)
write.table(strand[, colSums(!is.na(strand)) > 0, drop = FALSE], con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
