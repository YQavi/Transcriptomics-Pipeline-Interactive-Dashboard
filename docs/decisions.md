# Design decisions

Short records of the choices that shaped the project, including the ones I reversed. Each states the context, the choice, and the trade-off.

## 1. Plain numbered scripts instead of a workflow manager
- **Context.** I first built the pipeline in Snakemake. For three samples on one laptop it added a layer to learn and debug (for example, conda environments created inside a project path with a space broke a Bioconductor install).
- **Decision.** Five numbered scripts that share one settings file, each idempotent: finished work is skipped and files are renamed into place only on success.
- **Trade-off.** Less automatic dependency tracking and no cluster execution. If the project grew to dozens of samples or needed a scheduler, Snakemake or Nextflow would be the right move; the script boundaries map directly onto rules.

## 2. HISAT2 and featureCounts, not STAR or Salmon
- **Decision.** Alignment based counting with HISAT2 and featureCounts.
- **Why.** Low memory on a laptop, BAMs that can be inspected in a genome browser, and an interpretable alignment rate and assignment rate for QC.
- **Trade-off.** Salmon (alignment free) is faster and often at least as accurate for gene level counts; STAR is more sensitive for novel junctions. I would choose differently for a large cohort.

## 3. DESeq2 with apeglm, not edgeR or limma-voom
- **Decision.** DESeq2 for testing and apeglm for effect sizes.
- **Why.** A well documented negative binomial model with automatic handling of small designs, and shrinkage that fits the ranking and plotting use case.
- **Trade-off.** limma-voom becomes attractive with larger designs and complex covariates; edgeR gives similar results to DESeq2 here.

## 4. Raw p-value thresholds in the dashboard
- **Decision.** The volcano coloring uses raw p-value and shrunken log2FC with an adjustable p-value step of 0.001, and the dashboard shows the number of expected false positives.
- **Why.** With one late sample almost nothing will pass an adjusted threshold, but the dashboard should still allow exploration.
- **Trade-off.** Raw p-values overstate evidence. Every view carries a caveat and padj is always shown in the tables.

## 5. Wald statistic as the pathway ranking metric
- **Decision.** Rank genes by the Wald statistic, not by log2FC or p-value alone.
- **Why.** It combines effect size and uncertainty, so noisy low-count genes do not dominate the ranking.

## 6. A slim data object for the dashboard
- **Decision.** The pipeline writes a separate file of plain data frames and matrices (no DESeq2 objects, no local paths) that the app loads.
- **Why.** Small, hostable on a free tier, and the app does not need DESeq2 installed.
- **Trade-off.** The full DESeq2 object is saved separately for analysis work, so there are two files to keep consistent; a test checks the slim file's contents.

## 7. No normal tissue in the comparison
- **Decision.** Do not merge an external normal cervix dataset.
- **Why.** Bulk normal tissue has a different cell composition than purified fibroblasts, and a separate study confounds condition with batch.

## 8. Colorblind safe palette
- **Decision.** Okabe-Ito colors for groups and significance classes, with shape or labels used where meaning matters so color is never the only cue.
