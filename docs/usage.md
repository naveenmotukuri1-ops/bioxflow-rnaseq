# Usage

## 1. Install

You need [Nextflow](https://www.nextflow.io) 25.04 or newer (Java 17+) and one of
Docker, Singularity/Apptainer, or conda/mamba.

## 2. Samplesheet

A CSV file with one row per FASTQ pair:

```csv
sample,fastq_1,fastq_2,condition,batch
ctrl_1,/data/ctrl_1_R1.fastq.gz,/data/ctrl_1_R2.fastq.gz,control,day1
ctrl_2,/data/ctrl_2_R1.fastq.gz,/data/ctrl_2_R2.fastq.gz,control,day2
ctrl_3,/data/ctrl_3_R1.fastq.gz,/data/ctrl_3_R2.fastq.gz,control,day3
drug_1,/data/drug_1_R1.fastq.gz,/data/drug_1_R2.fastq.gz,treated,day1
drug_2,/data/drug_2_R1.fastq.gz,/data/drug_2_R2.fastq.gz,treated,day2
drug_3,/data/drug_3_R1.fastq.gz,/data/drug_3_R2.fastq.gz,treated,day3
```

| Column      | Required | Meaning |
|-------------|----------|---------|
| `sample`    | yes | Unique library name. Rows with the same name are treated as extra sequencing runs of one library and are concatenated. |
| `fastq_1`   | yes | Read 1 (gzipped FASTQ). Absolute path, URL, or a path relative to the launch directory or to the samplesheet. |
| `fastq_2`   | no  | Read 2. Leave empty for single-end data. |
| `condition` | yes | Experimental group, used for differential expression. Start with a letter; letters, numbers, `_` and `.` only. |
| `batch`     | no  | Batch or pairing variable. If present with more than one level the model becomes `~ batch + condition`. |

Strandedness is not a column: it is inferred from the data by two methods and reported.

## 3. Run

```bash
nextflow run bioxflow/rnaseq -profile docker \
    --input samplesheet.csv \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --gtf Homo_sapiens.GRCh38.112.gtf.gz \
    --outdir results
```

Use a genome FASTA and GTF from the **same source and release** (Ensembl or GENCODE).
For GENCODE files add `--gencode`.

### Profiles

| Profile | What it does |
|---|---|
| `docker`, `singularity`, `apptainer` | Run every tool from a pinned container image |
| `conda`, `mamba` | Create pinned conda environments instead |
| `arm` | Add to `docker` on Apple Silicon to run the amd64 images under emulation |
| `laptop` | Salmon-only mode without genome decoys, for machines with about 16 GB RAM |
| `test` | Run the bundled yeast test data (about 5 minutes) |

### Memory

| Mode | Human / mouse genome |
|---|---|
| Default (Salmon with decoys + STAR) | about 36 GB RAM |
| `--skip_alignment` (Salmon with decoys) | about 20 GB RAM |
| `-profile laptop` (Salmon, no decoys) | about 12 GB RAM |

In Salmon-only modes the concordance check is not available, because it needs the second
quantification from STAR.

## 4. Parameters

Run `nextflow run bioxflow/rnaseq --help` for the full list.

### Differential expression

* `--reference_condition control` sets the baseline. Default: the first condition in the samplesheet.
* By default every other condition is compared with the reference. For specific comparisons
  give `--contrasts contrasts.csv`:

  ```csv
  treatment,control
  drugA,vehicle
  drugB,vehicle
  drugB,drugA
  ```
* `--alpha 0.05` is the adjusted p-value cut-off; `--lfc_threshold 1` additionally requires
  an absolute shrunken log2 fold change of at least 1.
* Each condition needs at least 2 replicates to be tested; 3 or more are recommended.

### Tool arguments

Every tool's extra arguments are set in [`conf/modules.config`](../conf/modules.config) through
`ext.args`. Override them in your own config file and pass it with `-c my.config`:

```groovy
process {
    withName: 'FASTP'        { ext.args = '--length_required 36 --trim_poly_g' }
    withName: 'SALMON_QUANT' { ext.args = '--validateMappings --seqBias --numBootstraps 30' }
}
```

### Running on a cluster

Set the upper limits to the size of your largest node and choose an executor:

```groovy
// my.config
process.executor = 'slurm'
params.max_cpus   = 32
params.max_memory = '120.GB'
```

## 5. Reproducing a run

Every tool version is pinned. `results/pipeline_info/` records the software versions, the
execution trace and a full run report, and `differential_expression/R_sessionInfo.txt`
records the R package versions. Pin the pipeline itself with `-r <tag or commit>`.

## 6. Troubleshooting

**`there is no package called 'GenomeInfoDbData'` with `-profile conda` or `mamba`.**
This Bioconductor data package is downloaded by a conda post-install step that sometimes
fails silently. Remove the R environment from the conda cache (`work/conda/` by default, or
`$NXF_CONDA_CACHEDIR`) and re-run with `-resume` so it is created again. The container
profiles are not affected.
