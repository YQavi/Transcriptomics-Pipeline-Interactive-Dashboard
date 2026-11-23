# Analysis plan (pre-specified)

Commit this file **before** looking at results, then record any deviation in the "Changes" table at the bottom. It is a lightweight pre-registration. The template section is reusable for other projects.

Status key: **[done]** implemented and tested, **[planned]** specified but not yet run.

## Template (copy for a new project)

| Item | Question to answer |
|---|---|
| Hypothesis or question | What exactly is being asked? Exploratory or confirmatory? |
| Estimand | What quantity is being estimated (for example, log2 fold change late vs early)? |
| Unit of replication | What is an independent observation (donor, section, cell line)? |
| Sample size and power | How many units, and what effect could they detect? |
| Confounders and batches | What else differs between groups? What is balanced? |
| Primary analysis | Which test or model, with which settings? |
| Multiple testing | Which error rate is controlled, and how? |
| Thresholds | Cutoffs fixed in advance, and for what purpose? |
| Controls | Positive and negative controls |
| Sensitivity analyses | What would change the conclusion? |
| Validation | Independent data, simulation or orthogonal assay |
| Reporting | What is reported regardless of outcome? |

## This project

| Item | Plan | Status |
|---|---|---|
| Question | Which genes and pathways differ between CAFs from late and early stage FIGO cervical cancer? Exploratory, hypothesis generating. | n/a |
| Estimand | log2 fold change late vs early (shrunken for ranking), with the Wald test p-value for the unshrunken coefficient. | [done] |
| Unit of replication | One CAF sample per row of `samples.tsv` (probably one donor each). n = 2 early, 1 late. | [done] |
| Power | Report the replicates per group needed for 80% power at a 2-fold change for low, medium and high expressed genes, using RNASeqPower or PROPER. | [planned] |
| Confounders | Donor, passage and batch are confounded with stage and cannot be separated with this design. Stated as a limitation. | [done] |
| Primary analysis | DESeq2, design `~ condition`, reference `early_stage`, Wald test, apeglm shrinkage, genes with fewer than 10 total counts removed. | [done] |
| Multiple testing | Benjamini-Hochberg with DESeq2 independent filtering; padj reported for every gene. The dashboard thresholds on raw p-value for exploration and shows the expected number of false positives. | [done] / expected-false-positive counter [planned] |
| Thresholds | Dashboard defaults: raw p < 0.05, |log2FC| >= 2, adjustable in 0.001 steps. Static candidate table: p < 0.05 and |log2FC| >= 1.5. Used only to rank and explore, never as evidence of significance. | [done] |
| Positive controls | CAF markers (ACTA2, FAP, PDGFRB, COL1A1) expressed in all samples. | [planned] |
| Negative controls | Simulated data with known truth: the fraction of null genes with p < 0.05 should be near 5% (unit test); designs compared by simulation. | [done] (unit test), [planned] (full simulation run) |
| Pathway analysis | fgsea on the Wald statistic with MSigDB Hallmark plus custom CAF sets; sets of 5 to 500 genes; BH adjustment. Results interpreted as hypotheses. | [done] |
| Sensitivity analyses | (1) leave one early sample out and correlate fold changes; (2) shrunken vs unshrunken; (3) threshold sweep. | [planned] |
| Validation | None possible within this dataset. Next step is an independent dataset with replicates, compared at gene set level; qPCR on a few candidates. | [planned] |
| Reporting | All QC, the p-value histogram, the number of genes at each cutoff, and the single-sample limitation are reported whether or not anything looks interesting. | [done] |

## Changes after seeing the data

| Date | Change | Reason |
|---|---|---|
| | | |
