# Data dictionary

## Inputs

| File | Description |
|---|---|
| `config/samples.tsv` | Tab separated. Columns: `sample` (matches the BAM file name), `condition` (group label used by the design), `accession` (SRA run, falls back to `sample`). Extra columns are allowed and can be used as covariates in the design. |
| `config.sh` | Paths, library type (paired, strandedness) and thread counts for the shell steps. |
| `config/gene_sets/*.gmt` | Custom gene sets for pathway analysis, one per line: name, description, then genes. |
| `<DATA_DIR>/gene_counts.txt` | featureCounts output: a `#` command line, a header, then one row per gene (Geneid, Chr, Start, End, Strand, Length, one count column per BAM). |
| `<DATA_DIR>/gene_counts.txt.summary` | featureCounts assignment summary (Assigned, Unassigned_*). Optional for step 5. |
| `<DATA_DIR>/aligned_bam/<sample>_summary.txt` | HISAT2 alignment summary. Optional for step 5. |

## Outputs of step 5 (`results/deseq2/`)

| File | Description |
|---|---|
| `app_data.RData` | Slim object loaded by the dashboard (see below). Copied to `app/data/demo_deseq2_data.RData` for deployment. |
| `deseq2_data.RData` | Full object: `dds`, `res_df`, `norm_counts`, `vsd`, `colData`, `contrast_info`. |
| `deseq2_all_results.csv` | One row per tested gene, ranked by p-value (columns below). |
| `normalized_counts.csv` | Counts divided by the DESeq2 size factors (genes by samples). Used for per gene plots. |
| `vst_counts.csv` | Variance stabilized values (genes by samples). Used for PCA, distances and heatmaps. |
| `top_candidate_genes.csv` | Genes with p < `candidate_pvalue` and |log2FC| >= `candidate_lfc`. For exploration only. |
| `gsea_results.csv` | fgsea output (columns below). Empty if pathway analysis was skipped. |
| `run_qc.csv` | Per sample QC (columns below). Absent if no QC inputs were found. |
| `qc_pca_volcano.pdf` | Static PCA and volcano plot. |

## `deseq2_all_results.csv` and `res_df`

| Column | Meaning |
|---|---|
| `gene_id` | Identifier from the annotation (for example Ensembl). |
| `gene_name` | Unique display label: gene symbol, with the ID appended when a symbol is duplicated. |
| `symbol` | Plain gene symbol, used for pathway analysis. |
| `baseMean` | Mean of normalized counts across all samples. |
| `log2FoldChange` | Shrunken (apeglm) log2 fold change, numerator over denominator (late vs early). Positive means higher in late stage. |
| `lfcSE` | Posterior standard deviation of the shrunken log2 fold change. |
| `pvalue` | Wald test p-value of the unshrunken coefficient. NA for all-zero genes or Cook's distance outliers. |
| `padj` | Benjamini-Hochberg adjusted p-value. NA where `pvalue` is NA or the gene was removed by independent filtering. |
| `stat` | Wald statistic (unshrunken log2FC divided by its standard error). Used to rank genes for pathway analysis. |

## `run_qc.csv` and `qc_df`

| Column | Meaning |
|---|---|
| `sample`, `group` | Sample ID and its group. |
| `total_reads` | Input reads (read pairs for paired-end data) from the HISAT2 summary. |
| `alignment_rate` | HISAT2 overall alignment rate, percent. |
| `detected_genes` | Genes with at least one count in this sample. |
| `assigned_reads`, `assigned_pct` | featureCounts assigned reads and the percentage of all counted reads. |

## `gsea_results.csv` and `gsea_df`

| Column | Meaning |
|---|---|
| `pathway` | Gene set name. |
| `pval`, `padj` | fgsea p-value and Benjamini-Hochberg adjusted p-value. |
| `ES`, `NES` | Enrichment score and its normalized value. Positive NES means enriched among genes higher in the numerator group. |
| `size` | Number of genes of the set found in the ranking. |
| `leadingEdge` | Genes driving the enrichment, comma separated. |
| `source` | `MSigDB` or `Custom`. |

## Objects inside `app_data.RData`

`res_df`, `norm_counts`, `vst_mat`, `gsea_df`, `qc_df` (may be NULL), `contrast_info` (list: variable, numerator, denominator, design, title) and `colData` (sample table with the group variables). A unit test fails if any of these names changes.
