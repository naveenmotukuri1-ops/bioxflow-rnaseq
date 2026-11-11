# Output

```
results/
├── multiqc/multiqc_report.html        start here
├── concordance/
├── quantification/
├── differential_expression/
├── salmon/<sample>/
├── star/
├── fastp/  fastqc/
├── reference/
└── pipeline_info/
```

## multiqc/

`multiqc_report.html` gathers everything in one page: the concordance table and scatter plot,
inferred strandedness, PCA, sample distances, the differential expression summary, and the
FastQC, fastp, Salmon, STAR and samtools statistics.

## concordance/

The quality gate that is specific to this pipeline.

`concordance_summary.tsv`, one row per sample:

| Column | Meaning |
|---|---|
| `genes_compared` | Genes with at least 10 reads by either method |
| `spearman` | Rank correlation of Salmon and STAR gene counts |
| `pearson_log2` | Pearson correlation of log2(count + 1) |
| `pct_within_2fold` | Percentage of genes whose two counts differ by less than 2-fold |
| `salmon_assigned`, `star_assigned` | Reads assigned to genes by each method |
| `strand_agreement` | Whether the two methods infer the same library strandedness |
| `status` | `PASS`, or `REVIEW` if `spearman` is below `--concordance_threshold` or the strandedness calls disagree |

`concordance_scatter.png` plots the two counts per gene for every sample.
`strandedness.tsv` lists the strandedness inferred by Salmon and by STAR.

How to read it: Salmon and STAR use different algorithms and different evidence. When they
agree, the gene counts do not depend on the choice of method. A sample marked `REVIEW`
usually points to a wrong or mismatched annotation, contamination, heavy degradation or a
sample swap, and should be looked at before trusting its differential expression results.
Some disagreement is expected for gene families with near-identical copies, where STAR
discards multi-mapping reads and Salmon assigns them probabilistically.

## quantification/

| File | Content |
|---|---|
| `salmon.gene_counts.tsv` | Estimated read counts per gene (used for DESeq2) |
| `salmon.gene_tpm.tsv` | Transcripts per million per gene |
| `salmon.gene_lengths.tsv` | Effective gene length per sample |
| `salmon.transcript_counts.tsv`, `salmon.transcript_tpm.tsv` | Transcript-level estimates |
| `star.gene_counts.tsv` | Uniquely mapped read counts per gene from STAR |
| `tx2gene.tsv` | Transcript to gene map built from the GTF |

## differential_expression/

| File | Content |
|---|---|
| `<treatment>_vs_<control>.deseq2_results.tsv` | One row per tested gene, sorted by adjusted p-value |
| `<treatment>_vs_<control>.volcano.png`, `.ma.png` | Volcano and MA plots |
| `de_summary.tsv` | Number of up- and down-regulated genes per contrast |
| `normalized_counts.tsv` | DESeq2 normalised counts |
| `transformed_counts.tsv` | Variance-stabilised counts for clustering and heatmaps |
| `pca.png`, `pca.tsv`, `sample_distances.png` | Sample-level QC |
| `analysis_notes.txt` | Model design, filters and any fallback that was used |
| `R_sessionInfo.txt` | Exact R package versions |

Columns of the results table:

| Column | Meaning |
|---|---|
| `baseMean` | Mean normalised count across samples |
| `log2FoldChange`, `lfcSE` | Maximum-likelihood estimate and its standard error |
| `stat`, `pvalue`, `padj` | Wald statistic, p-value, Benjamini-Hochberg adjusted p-value |
| `log2FoldChange_shrunk`, `lfcSE_shrunk` | apeglm-shrunken effect size; use this for ranking and plots |
| `significant` | `up`, `down` or `no` at the chosen `--alpha` and `--lfc_threshold` |

A positive fold change means higher expression in the treatment than in the control.

## salmon/, star/, fastp/, fastqc/

Per-sample tool outputs: Salmon `quant.sf` and logs, STAR sorted BAM with index, splice
junctions and mapping statistics, fastp trimming reports, FastQC reports on the raw reads.

## reference/

The Salmon and STAR indices built for this run. Reuse them with `--salmon_index` and
`--star_index` to save time.

## pipeline_info/

Nextflow execution report, timeline and trace, plus the versions of all tools used.
