#!/usr/bin/env Rscript
# bioxflow/rnaseq — sample-level QC and differential expression with DESeq2.
#
# Counts come from tximport, so DESeq2 uses per-gene, per-sample average transcript
# length offsets. Each contrast reports Wald-test p-values (BH adjusted) together
# with apeglm-shrunken log2 fold changes for ranking and plotting.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
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
opt <- parse_args(list(txi = NULL, samples = NULL, tx2gene = NULL, contrasts = NULL, reference = NULL,
                       alpha = "0.05", lfc_threshold = "0", min_count = "10", outdir = "deseq2"))
alpha <- as.numeric(opt$alpha)
lfc_threshold <- as.numeric(opt$lfc_threshold)
min_count <- as.numeric(opt$min_count)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(...) file.path(opt$outdir, ...)

txi     <- readRDS(opt$txi)
samples <- read.delim(opt$samples, stringsAsFactors = FALSE, colClasses = "character")
tx2gene <- read.delim(opt$tx2gene, stringsAsFactors = FALSE, colClasses = "character")
gene_name <- function(ids) tx2gene$gene_name[match(ids, tx2gene$gene_id)]

stopifnot(identical(colnames(txi$counts), samples$sample))
conditions <- unique(samples$condition)
reference  <- if (is.null(opt$reference) || opt$reference %in% c("", "null")) conditions[1] else opt$reference
if (!reference %in% conditions) stop("reference condition '", reference, "' is not in the samplesheet")

samples$condition <- relevel(factor(samples$condition, levels = conditions), ref = reference)
use_batch <- length(setdiff(unique(samples$batch), "NA")) > 1
if (use_batch) samples$batch <- factor(samples$batch)
design <- if (use_batch) ~ batch + condition else ~ condition
rownames(samples) <- samples$sample

notes <- character()
note <- function(...) { msg <- sprintf(...); message(msg); notes <<- c(notes, msg) }
note("Design: %s; reference condition: %s", paste(deparse(design), collapse = ""), reference)

multi_condition <- length(conditions) > 1
dds <- DESeqDataSetFromTximport(txi, colData = samples,
                                design = if (multi_condition) design else ~ 1)

# keep genes with at least `min_count` reads in as many samples as the smallest group
smallest <- min(table(samples$condition))
keep <- rowSums(counts(dds) >= min_count) >= smallest
note("Genes kept after low-count filter (>= %d reads in >= %d samples): %d of %d",
     min_count, smallest, sum(keep), length(keep))
dds <- dds[keep, ]

fit_deseq <- function(d) {
  tryCatch(DESeq(d, quiet = TRUE), error = function(e) {
    note("Parametric dispersion fit failed (%s); using fitType = 'mean'", conditionMessage(e))
    tryCatch(DESeq(d, fitType = "mean", quiet = TRUE), error = function(e2) {
      note("Mean dispersion fit failed; using gene-wise dispersion estimates")
      d <- estimateSizeFactors(d)
      d <- estimateDispersionsGeneEst(d)
      dispersions(d) <- mcols(d)$dispGeneEst
      nbinomWaldTest(d)
    })
  })
}

# ------------------------------------------------------------ sample QC ----
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(norm), gene_name = gene_name(rownames(norm)), round(norm, 2), check.names = FALSE),
            out("normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

vsd <- tryCatch(varianceStabilizingTransformation(dds, blind = TRUE), error = function(e) {
  note("VST failed (%s); using log2(normalised count + 1) for QC plots", conditionMessage(e))
  normTransform(dds)
})
vmat <- assay(vsd)
write.table(data.frame(gene_id = rownames(vmat), gene_name = gene_name(rownames(vmat)), round(vmat, 4), check.names = FALSE),
            out("transformed_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (ncol(vmat) >= 3 && nrow(vmat) >= 3) {
  top <- head(order(apply(vmat, 1, var), decreasing = TRUE), 500)
  pca <- prcomp(t(vmat[top, , drop = FALSE]))
  pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pd  <- data.frame(sample = colnames(vmat), condition = samples$condition, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  write.table(pd, out("pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  p <- ggplot(pd, aes(PC1, PC2, colour = condition, label = sample)) +
    geom_point(size = 3.5) +
    labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
         title = "Principal component analysis", subtitle = "500 most variable genes, variance-stabilised counts") +
    theme_bw(base_size = 12)
  if (requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE)
  ggsave(out("pca.png"), p, width = 7, height = 5.2, dpi = 130)
  file.copy(out("pca.png"), "bioxflow_pca_mqc.png")

  dmat <- as.matrix(dist(t(vmat)))
  ann  <- data.frame(condition = samples$condition, row.names = colnames(vmat))
  pheatmap(dmat, annotation_col = ann, main = "Sample-to-sample distances",
           color = colorRampPalette(c("#0e3b43", "#5fb4bd", "#f4fbfc"))(60),
           filename = out("sample_distances.png"), width = 7, height = 6)
  file.copy(out("sample_distances.png"), "bioxflow_sample_distances_mqc.png")
}

# ------------------------------------------------ differential expression ----
contrasts <- NULL
if (multi_condition) {
  if (!is.null(opt$contrasts)) {
    contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE, strip.white = TRUE)
    if (!all(c("treatment", "control") %in% names(contrasts))) stop("contrasts file needs columns: treatment,control")
    bad <- setdiff(c(contrasts$treatment, contrasts$control), conditions)
    if (length(bad)) stop("contrast conditions not in samplesheet: ", paste(bad, collapse = ", "))
  } else {
    contrasts <- data.frame(treatment = setdiff(conditions, reference), control = reference, stringsAsFactors = FALSE)
  }
}

summary_rows <- list()
if (is.null(contrasts)) {
  note("Only one condition in the samplesheet: differential expression skipped")
} else if (min(table(samples$condition)) < 2) {
  note("Differential expression skipped: every condition needs at least 2 replicates (3 or more recommended)")
} else {
  if (min(table(samples$condition)) < 3)
    note("Fewer than 3 replicates in at least one condition: statistical power is limited")
  dds <- fit_deseq(dds)

  for (k in seq_len(nrow(contrasts))) {
    trt <- contrasts$treatment[k]; ctl <- contrasts$control[k]
    name <- paste0(trt, "_vs_", ctl)
    res <- results(dds, contrast = c("condition", trt, ctl), alpha = alpha)

    # shrunken effect sizes: refit with the contrast's control as reference level
    shr_type <- "apeglm"
    shr <- tryCatch({
      d2 <- dds
      d2$condition <- relevel(d2$condition, ref = ctl)
      d2 <- nbinomWaldTest(d2, quiet = TRUE)
      lfcShrink(d2, coef = paste0("condition_", trt, "_vs_", ctl), type = "apeglm", quiet = TRUE)
    }, error = function(e) {
      note("%s: apeglm shrinkage failed (%s); using DESeq2 'normal' shrinkage", name, conditionMessage(e))
      shr_type <<- "normal"
      lfcShrink(dds, contrast = c("condition", trt, ctl), type = "normal", quiet = TRUE)
    })

    tab <- data.frame(
      gene_id = rownames(res), gene_name = gene_name(rownames(res)),
      baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
      stat = res$stat, pvalue = res$pvalue, padj = res$padj,
      log2FoldChange_shrunk = shr$log2FoldChange, lfcSE_shrunk = shr$lfcSE,
      stringsAsFactors = FALSE)
    sig <- !is.na(tab$padj) & tab$padj < alpha & abs(tab$log2FoldChange_shrunk) >= lfc_threshold
    tab$significant <- ifelse(sig, ifelse(tab$log2FoldChange > 0, "up", "down"), "no")
    tab <- tab[order(tab$padj, -abs(tab$log2FoldChange_shrunk)), ]
    write.table(tab, out(paste0(name, ".deseq2_results.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

    cols <- c(up = "#c2410c", down = "#0e7c86", no = "grey75")
    pd <- tab[!is.na(tab$padj), ]
    v <- ggplot(pd, aes(log2FoldChange_shrunk, -log10(padj), colour = significant)) +
      geom_point(size = 1, alpha = 0.7) +
      scale_colour_manual(values = cols, breaks = c("up", "down", "no")) +
      geom_hline(yintercept = -log10(alpha), linetype = 2, colour = "grey50") +
      labs(x = sprintf("log2 fold change (%s shrunk)", shr_type), y = "-log10 adjusted p-value",
           title = sprintf("%s vs %s", trt, ctl),
           subtitle = sprintf("%d up, %d down at padj < %s", sum(tab$significant == "up"), sum(tab$significant == "down"), alpha)) +
      theme_bw(base_size = 12)
    if (requireNamespace("ggrepel", quietly = TRUE) && any(pd$significant != "no"))
      v <- v + ggrepel::geom_text_repel(data = head(pd[pd$significant != "no", ], 15), aes(label = gene_name),
                                        size = 2.8, show.legend = FALSE, max.overlaps = 30)
    ggsave(out(paste0(name, ".volcano.png")), v, width = 7, height = 5.5, dpi = 130)
    file.copy(out(paste0(name, ".volcano.png")), paste0("bioxflow_volcano_", name, "_mqc.png"))

    m <- ggplot(tab[tab$baseMean > 0, ], aes(baseMean, log2FoldChange_shrunk, colour = significant)) +
      geom_point(size = 1, alpha = 0.7) + scale_x_log10() +
      scale_colour_manual(values = cols, breaks = c("up", "down", "no")) +
      geom_hline(yintercept = 0, colour = "grey40") +
      labs(x = "mean of normalised counts", y = sprintf("log2 fold change (%s shrunk)", shr_type),
           title = sprintf("MA plot: %s vs %s", trt, ctl)) +
      theme_bw(base_size = 12)
    ggsave(out(paste0(name, ".ma.png")), m, width = 7, height = 5, dpi = 130)

    summary_rows[[name]] <- data.frame(contrast = name, genes_tested = sum(!is.na(tab$padj)),
                                       up = sum(tab$significant == "up"), down = sum(tab$significant == "down"),
                                       alpha = alpha, lfc_threshold = lfc_threshold, shrinkage = shr_type,
                                       stringsAsFactors = FALSE)
  }
}

if (length(summary_rows)) {
  smry <- do.call(rbind, summary_rows)
  write.table(smry, out("de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  header <- c(
    "# id: 'bioxflow_de_summary'",
    "# section_name: 'Differential expression summary'",
    "# description: 'DESeq2 Wald test with Benjamini-Hochberg correction. Fold changes are apeglm-shrunken estimates.'",
    "# plot_type: 'table'",
    "# pconfig:",
    "#     id: 'bioxflow_de_summary_table'",
    "#     title: 'Differentially expressed genes'")
  con <- file("bioxflow_de_summary_mqc.tsv", "w")
  writeLines(header, con)
  write.table(smry, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}

writeLines(notes, out("analysis_notes.txt"))
writeLines(capture.output(sessionInfo()), out("R_sessionInfo.txt"))
if (!length(list.files(".", pattern = "_mqc\\."))) writeLines(notes, "bioxflow_de_notes_mqc.txt")
