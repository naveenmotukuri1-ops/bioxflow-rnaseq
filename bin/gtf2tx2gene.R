#!/usr/bin/env Rscript
# Build a transcript -> gene map (transcript_id, gene_id, gene_name) from a GTF file.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: gtf2tx2gene.R <genes.gtf> <tx2gene.tsv>")

gtf <- read.delim(args[1], header = FALSE, comment.char = "#", quote = "",
                  stringsAsFactors = FALSE, colClasses = "character")
if (ncol(gtf) < 9) stop("GTF file must have 9 tab-separated columns")

feature <- if (any(gtf$V3 == "transcript")) "transcript" else "exon"
attr <- unique(gtf$V9[gtf$V3 == feature])

get_attr <- function(key, x) {
  pattern <- paste0('^(.*; *)?', key, ' "([^"]*)".*$')
  ifelse(grepl(pattern, x), sub(pattern, "\\2", x), NA_character_)
}

tx <- data.frame(
  transcript_id = get_attr("transcript_id", attr),
  gene_id       = get_attr("gene_id", attr),
  gene_name     = get_attr("gene_name", attr),
  stringsAsFactors = FALSE
)
tx <- tx[!is.na(tx$transcript_id) & !is.na(tx$gene_id), ]
tx$gene_name[is.na(tx$gene_name) | tx$gene_name == ""] <- tx$gene_id[is.na(tx$gene_name) | tx$gene_name == ""]
tx <- tx[!duplicated(tx$transcript_id), ]
if (nrow(tx) == 0) stop("No transcript_id / gene_id attributes found in the GTF")

write.table(tx, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("tx2gene: %d transcripts, %d genes", nrow(tx), length(unique(tx$gene_id))))
