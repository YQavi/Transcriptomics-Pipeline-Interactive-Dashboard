# Bulk RNA-seq of Cancer Associated Fibroblasts in Cervical Cancer: Pipeline and Interactive Dashboard

[![CI](https://github.com/YQavi/bulk-rnaseq-cervical-caf/actions/workflows/ci.yml/badge.svg)](https://github.com/YQavi/bulk-rnaseq-cervical-caf/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

An end to end bulk RNA-seq project, from raw SRA reads to an interactive R Shiny dashboard, built on public dataset **GSE347724** (cancer associated fibroblasts, CAFs, from late stage vs early stage FIGO cervical cancer). The pipeline is a set of numbered scripts, one per step, that runs on a single laptop and is configured through one settings file and one sample sheet so it can be reused on other datasets.

> **Live dashboard:** `<https://yqavi.shinyapps.io/cervical-caf-rnaseq-dashboard/>`
> **Author:** Yusuf Qavi, M.S. Bioinformatics, Georgetown University

<!-- TODO: add a 10 second GIF of the click-to-boxplot interaction at docs/img/dashboard_demo.gif -->
![Dashboard demo](docs/img/dashboard_demo.gif)

---

## Why this project

Differential expression results usually end up as a static volcano plot and a spreadsheet. That format hides the question a biologist actually asks next: *"What does this particular gene look like across my samples?"* This project closes that gap. The pipeline produces the statistics, and the Shiny dashboard lets a reader move from the genome wide view to a single gene in one click.

It also demonstrates the practical side of the field: working with a dataset that has no accompanying paper, a very small sample size, and a laptop instead of a cluster, and being explicit about what the data can and cannot support.

## Pipeline at a glance

```mermaid
flowchart LR
    A[SRA accessions] --> B[prefetch + fasterq-dump]
    B --> C[fastp QC and trimming]
    C --> D[HISAT2 alignment]
    D --> E[samtools sort + index]
    E --> F[featureCounts]
    F --> G[DESeq2]
    G --> H[Shiny dashboard]
```

Each step is its own script in `scripts/`. Every step skips work that is already done, so the whole chain can be re-run safely, or any single step can be run on its own.

## Design decisions

| Step | Tool | Why this choice |
|---|---|---|
| Download | `prefetch` + `fasterq-dump --split-files` | Fetches the compact `.sra` once, then converts to paired FASTQ locally. Temporary files go to fast internal storage and outputs are gzipped to save space. |
| QC and trimming | `fastp` | One pass handles adapter detection (automatic for paired end), quality filtering (Q20) and produces per sample HTML and JSON reports. |
| Alignment | `HISAT2` | A splice aware aligner is required for RNA-seq. HISAT2 needs far less memory than STAR for the human genome, which matters on a laptop. |
| Sorting | `samtools sort` (piped) | Alignments stream straight from HISAT2 into the sorter, so the large intermediate SAM file is never written to disk. |
| Quantification | `featureCounts` | Fast gene level counting. Strandedness (`-s`) is exposed in `config.sh` because the wrong value can silently discard most reads. |
| Statistics | `DESeq2` | RNA-seq counts are overdispersed integers, so a negative binomial model with median-of-ratios normalization is more appropriate than a per gene t-test, especially with few replicates. |
| Effect sizes | `apeglm` shrinkage | Low count genes produce wildly inflated fold changes. Shrinkage gives stable log2 fold changes for ranking and plotting. The original Wald p-values are kept alongside. |
| Visualization | VST counts | Variance stabilized values are used for PCA and sample distance heatmaps, while normalized counts are used for per gene plots. |

### Reproducibility

- **Numbered scripts, resumable runs.** Each step writes to a temporary name and renames on success, then skips itself next time, so an interrupted run never leaves a half-written file that looks finished.
- **One environment file** (`environment.yml`) pins the command line tools and the R packages.
- **Config driven:** library type, strandedness, reference and thread counts live in `config.sh`. Samples and groups live in `config/samples.tsv`. The DESeq2 design, contrast and organism are one settings block at the top of `scripts/deseq2.R`. No step needs editing to change dataset.
- **Safety checks:** missing inputs, a missing HISAT2 index and failures on either side of a pipe are reported per sample. The DESeq2 step verifies that sample IDs and contrast levels exist and warns when any group has fewer than two replicates.

## The Shiny dashboard

Static figures answer one question. An interactive app lets a reader ask many, which is why the results are delivered as a Shiny application rather than a PDF.

| Tab | What it does |
|---|---|
| **Sample QC and PCA** | Interactive PCA and a sample to sample distance heatmap to check that replicates cluster before trusting any gene level result. |
| **Run QC** | HISAT2 alignment rate and featureCounts assignment rate per sample, with a warning if assignment is low (a classic sign of the wrong strandedness setting). |
| **P-value Histogram** | Distribution of raw p-values, the standard diagnostic for whether a test is calibrated and whether there is any signal above the null. |
| **Volcano Plot** | WebGL scatter of every gene. Fold change and p-value cutoffs are live sliders, with the p-value cutoff moving in steps of 0.001. Top genes are labeled and any gene can be highlighted by name. |
| **MA Plot** | Mean expression against fold change, to separate real effects from low count noise. Clicking a point behaves like the volcano plot. |
| **Heatmap** | Row z-scored VST expression of the top genes passing the current cutoffs. |
| **Click to inspect** | Clicking any gene on the volcano or MA plot instantly draws that gene's expression across all samples in a side panel. |
| **Gene Inspector** | Search by gene symbol (mapped from Ensembl IDs), see log2FC, p-value and adjusted p-value, and download the per sample plot as a PNG. |
| **Pathways (GSEA)** | fgsea on the Wald statistic against MSigDB Hallmark plus custom CAF marker sets, with a bar plot of top enriched sets and a downloadable table. |
| **Results Table** | Searchable, sortable table with a selector for up regulated, down regulated or not significant genes, and CSV download. |
| **Methods and Caveats** | The statistical approach and the single sample limitation, kept inside the app so it travels with every screenshot and link. |

Design notes:
- **Reactive by construction.** Changing a threshold recomputes classification, colors, cutoff lines and the table in one pass, with no reruns of R scripts.
- **Dataset agnostic.** Group names, colors, titles and the comparison label are read from the pipeline output, so the same app works for any two group contrast.
- **Light to host.** The pipeline writes a slim data object (no DESeq2 classes, no local paths) so the app deploys with a small package footprint.
- **Colorblind safe.** Group and significance colors use the Okabe-Ito palette.
- **Honest by default.** When a group has a single sample, the sidebar shows a warning that p-values are exploratory, and it reports how many genes would pass the current p-value cutoff by chance alone.

Run it locally from the repository root:

```r
shiny::runApp("app")
```

Deploy to shinyapps.io with `Rscript deploy.R` after step 5 has produced `app/data/demo_deseq2_data.RData` (see the comments in that file for one time account setup). A smoke test lives in `tests/` and can be run with `Rscript tests/testthat.R` (needs `shinytest2`).

For reproducible R packages, run `renv::init()` then `renv::snapshot()` once and commit the resulting `renv.lock`.

## Experimental design and statistics

The design has one late stage and two early stage samples, so the statistical thinking is part of the deliverable, not an afterthought.

- **Pre-specified plan.** [`docs/analysis_plan.md`](docs/analysis_plan.md) fixes the estimand, thresholds, multiple testing rule and sensitivity analyses before the results are examined, and lists any later change.
- **Why replication, not the test, matters.** `analysis/simulate_design_power.R` simulates counts with known differential genes and compares a per gene t-test with DESeq2 from 1 vs 2 up to 10 vs 10 samples, reporting false positive rate, power and false discovery proportion. Run it with `make simulate`.
- **Calibration is tested.** The automated tests check, on simulated data, that fold changes have the right direction, that strong effects are detected, and that about 5% of null genes have p < 0.05.
- **Multiple testing is visible.** The dashboard shows the number of genes expected to pass the cutoff by chance and a p-value histogram.

## Working within hardware limits

This was developed on an Apple M1 laptop with the project stored on an external drive. The pipeline contains several deliberate workarounds that generalize to any resource limited setup.

| Constraint | Workaround |
|---|---|
| External drive is slow for small random writes | Scratch space (`SCRATCH_DIR`) points to internal storage for SRA dumps and sort temp files. Only final outputs touch the external drive. |
| Large intermediate files | HISAT2 output is piped directly into `samtools sort`, so no uncompressed SAM is written. FASTQs are gzipped immediately. |
| Limited RAM | HISAT2 is used instead of STAR. Sort memory is capped per thread (`SAMTOOLS_MEM_PER_THREAD`) so peak usage is predictable. |
| Building an index is memory heavy | A prebuilt HISAT2 index is used, and step 3 stops immediately with a clear message if it is missing. |
| Limited cores | Thread counts per tool are settings in `config.sh`. Steps run one after another, so tools never compete for the same cores. |
| Laptop sleep and drive unmount mid run | Long commands run under `caffeinate` on macOS so the machine stays awake, and missing inputs fail loudly at the start instead of deep in a run. Finished samples are skipped, so a rerun resumes where it stopped. |
| Silent failures in piped commands | Scripts use `pipefail` and check the exit status of both HISAT2 and `samtools sort`, so a failure on either side marks that sample as failed. |
| Apple Silicon package gaps | If a bioconda package has no `osx-arm64` build, create the environment with `CONDA_SUBDIR=osx-64`. |
| Spaces in the project path (common on external drives) | Every script quotes its paths. Create the conda environment in its normal location, not inside the project folder, because some Bioconductor data packages run install scripts that break on spaces. |

Free disk space is the usual real bottleneck: BAMs and sort scratch can need two to three times the size of the input FASTQs.

## Reproduce it

```bash
# 1. Create the environment (in your normal conda location)
conda env create -f environment.yml && conda activate txnomics

# 2. Edit config.sh (paths, paired/single end, strandedness) and config/samples.tsv
#    Put the HISAT2 index and a matching GTF where config.sh expects them

# 3. Run every step, or run them one at a time
caffeinate -i bash scripts/run_all.sh      # drop "caffeinate -i" on Linux
bash scripts/04_count_featurecounts.sh     # e.g. a single step
bash scripts/run_all.sh 5                  # start from step 5 (DESeq2) onward

# 4. Launch the dashboard
Rscript -e 'shiny::runApp("app")'
```

| Step | Script | Output |
|---|---|---|
| 1 | `scripts/01_download_sra.sh` | `<sample>_1.fastq.gz`, `_2.fastq.gz` |
| 2 | `scripts/02_trim_fastp.sh` | `clean_fastq/`, `reports/` (fastp HTML and JSON) |
| 3 | `scripts/03_align_hisat2.sh` | `aligned_bam/<sample>.sorted.bam` plus HISAT2 summaries |
| 4 | `scripts/04_count_featurecounts.sh` | `gene_counts.txt` and its `.summary` |
| 5 | `scripts/05_run_deseq2.sh` (runs `scripts/deseq2.R`) | `results/deseq2/` and `app/data/demo_deseq2_data.RData` |

Before running on a new dataset, confirm **paired vs single end**, **strandedness**, and that the **GTF matches the genome build** of the index. Step 4 prints the assigned read percentage for each sample, which is the quickest check that the strandedness setting is right.

## Repository structure

```
bulk-rnaseq-cervical-caf/
├── README.md
├── LICENSE                   # MIT
├── CITATION.cff
├── Makefile                  # make deseq2 | app | test | simulate | deploy
├── config.sh                 # paths, library type, threads (shared by the shell steps)
├── environment.yml           # conda environment: tools + R packages
├── deploy.R                  # shinyapps.io deployment
├── config/
│   ├── samples.tsv           # sample, condition, accession, optional covariates
│   └── gene_sets/            # custom GMT gene sets for pathway analysis
├── scripts/
│   ├── 01_download_sra.sh
│   ├── 02_trim_fastp.sh
│   ├── 03_align_hisat2.sh
│   ├── 04_count_featurecounts.sh
│   ├── 05_run_deseq2.sh      # wrapper around deseq2.R
│   ├── deseq2.R              # DE, GSEA, run QC, outputs, static QC plots
│   └── run_all.sh            # runs steps 1 to 5 in order
├── analysis/
│   └── simulate_design_power.R   # t-test vs DESeq2 across replicate numbers
├── app/
│   ├── app.R                 # Shiny dashboard
│   └── data/                 # slim demo_deseq2_data.RData for the hosted app
├── tests/
│   ├── test_shell_steps.sh   # stubbed test of steps 1 to 5
│   ├── ci_environment.yml
│   └── testthat/             # DESeq2 step on simulated data, dashboard smoke test
├── docs/
│   ├── WRITEUP.md  analysis_plan.md  decisions.md  data_dictionary.md  troubleshooting.md
│   └── img/                  # screenshots and demo GIF used in this README
├── .github/workflows/ci.yml  # continuous integration
└── results/                  # generated outputs (not tracked)
```

## Testing

```bash
make test          # shell steps (stubbed tools) and the DESeq2 step on simulated data
make test-shell    # runs in seconds, needs no bioinformatics tools
```

The shell test replaces the real tools with stubs and checks paths with spaces, skipping of finished work, failure detection on both sides of the HISAT2 and samtools pipe, strandedness and single-end options, and step 5 inputs. The R test generates simulated featureCounts output and checks the outputs, the structure of the dashboard data object, the statistics, and that bad inputs fail loudly. GitHub Actions runs both on every push.

## Documentation

| Document | What it covers |
|---|---|
| [`docs/WRITEUP.md`](docs/WRITEUP.md) | Framework for the full write-up: tool by tool methods and reasoning, design and statistics, results scaffolding, a defensibility checklist and references |
| [`docs/analysis_plan.md`](docs/analysis_plan.md) | Pre-specified analysis plan and a reusable template |
| [`docs/decisions.md`](docs/decisions.md) | Design decisions and trade-offs, including choices that were reversed |
| [`docs/data_dictionary.md`](docs/data_dictionary.md) | Every input and output file and column |
| [`docs/troubleshooting.md`](docs/troubleshooting.md) | Problems hit while building this and their fixes |

## Limitations and next steps

- The custom CAF gene sets are illustrative marker lists and should be replaced or checked against published signatures before any claim is made from them.
- Validate candidate genes against an independent cervical cancer or CAF dataset.
- Add a donor level covariate to the design if more samples become available.

## Data and acknowledgments

Data: NCBI GEO accession [GSE347724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347724). This dataset has no associated publication, so no original analysis is being reproduced here. Tools: SRA Toolkit, fastp, HISAT2, samtools, Subread (featureCounts), DESeq2, apeglm, fgsea, R Shiny, plotly.

## License and citation

Released under the [MIT License](LICENSE). If you use the pipeline or dashboard, please cite it using [`CITATION.cff`](CITATION.cff).
