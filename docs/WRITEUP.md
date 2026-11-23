# Write-up framework: Bulk RNA-seq of cancer associated fibroblasts (GSE347724)

**How to use this document.** Sections without tags are finished explanations of the methods and are meant to be kept. Anything tagged **[RESULT]** needs a number, figure or table from your own run, and **[FILL]** needs a fact you must look up (for example from the GEO sample pages). Delete this box and the tags when you publish. Aim for a final length of 3,000 to 4,500 words plus figures. Write for a reader who knows biology but has not used these tools.

The goal of the write-up is not to show that the pipeline ran. It is to show that you understand what each tool assumes, where it fails, and how you checked it.

---

## 1. Title and abstract

**Title.** Bulk RNA-seq of cancer associated fibroblasts from early and late stage cervical cancer: an end to end pipeline and interactive dashboard.

**Abstract template (150 words).**
> Cancer associated fibroblasts (CAFs) shape the tumor microenvironment, but public data on how they differ by disease stage is sparse. I processed the three CAF samples of GEO series GSE347724 (one late stage and two early stage FIGO cervical cancer) from raw reads to differential expression using fastp, HISAT2, featureCounts, DESeq2 and apeglm, and tested pathway enrichment with fgsea. Because the design has a single late stage sample, I treat all inference as exploratory, quantify what that design can and cannot support by simulation, and deliver results through an R Shiny dashboard that links genome wide views to single gene evidence. **[RESULT: one sentence on the main finding, for example the number of genes at a stated cutoff and the top pathway direction.]** The repository is fully scripted, tested and runs on a laptop.

**Key points (3 to 5 bullets).**
- What the data can support: **[RESULT]**.
- The most defensible biological observation: **[RESULT]**.
- What the simulation shows about 1 vs 2 designs: **[RESULT from analysis/simulate_design_power.R]**.
- What I would need to make a real claim: at least three independent donors per stage.

---

## 2. Background and question

- **CAFs.** Fibroblasts that live in tumor stroma. They are heterogeneous: contractile, matrix producing "myCAF" like states and inflammatory "iCAF" like states have been described in several cancers. They remodel extracellular matrix, secrete cytokines and influence immune infiltration and drug penetration. **[FILL: two or three citations on CAFs in cervical cancer specifically.]**
- **FIGO stage.** The International Federation of Gynecology and Obstetrics staging describes how far the cancer has spread. "Early" and "late" here are the labels given in the GEO sample records. **[FILL: the exact FIGO stages from the GEO records, and whether the samples came from different patients.]**
- **Question.** Which genes and pathways differ between CAFs from late and early stage disease? I frame this as a hypothesis generating question, not a test of a specific hypothesis.
- **Why the lack of a paper matters.** There is no published analysis to compare against, so the pipeline must justify itself with internal checks (sample QC, positive controls, simulation) rather than by reproducing a known result.

---

## 3. Data

| Run accession | GEO sample | Label in GEO | Group used here | Patient / donor | Sequencing |
|---|---|---|---|---|---|
| SRR40723926 | **[FILL]** | CAFs from late stage FIGO cervical cancer-1 | late_stage | **[FILL]** | **[FILL: platform, read length]** |
| SRR40723927 | **[FILL]** | CAFs from early stage FIGO cervical cancer (replicate) | early_stage | **[FILL]** | **[FILL]** |
| SRR40723928 | **[FILL]** | CAFs from early stage FIGO cervical cancer (replicate) | early_stage | **[FILL]** | **[FILL]** |

Confirm which accession belongs to which GEO sample before publishing; the mapping in `config/samples.tsv` is the one thing the code cannot verify. Also record **[FILL]**: paired or single end, library prep kit, strandedness, and whether the cells were cultured and for how many passages (passage number is a technical factor that can masquerade as biology).

---

## 4. Experimental design and statistical analysis plan

This section is where the project shows statistical judgment. The full pre-specified plan is in [`analysis_plan.md`](analysis_plan.md); the reasoning is below.

### 4.1 Unit of replication and what the data can support

The unit of replication is the **independent biological sample** (ideally one donor), not the number of reads or genes. Sequencing deeper adds precision about each sample but does not add replicates.

The design has two early stage samples and one late stage sample. DESeq2 will run on this because it estimates each gene's dispersion by borrowing information across thousands of genes (section 5.6). The test therefore *works mechanically*. What it cannot do is tell you how much the late stage group varies from person to person, because there is only one late stage sample. The reported p-values describe how surprising the late sample looks against the early samples and the shared dispersion trend, **not** how likely the difference is to generalize to other patients.

### 4.2 Confounding

Stage is almost certainly confounded with donor, and possibly with culture conditions, passage and batch. With one sample in one group these cannot be separated, whatever the statistics. State this plainly in the limitations, and do not use the word "significant" without "exploratory".

### 4.3 Power and sample size (what a valid design would need)

Power for RNA-seq depends on the fold change you want to detect, the biological variability (dispersion), the mean expression of the genes, the sequencing depth and the number of replicates. Report a small table of the replicates per group needed for 80% power at a 2-fold change, using a published tool (**[FILL]**: `RNASeqPower`, `PROPER` or `ssizeRNA`), for a low, medium and high expressed gene and for a typical human dispersion (biological coefficient of variation of about 0.4 for patient samples, lower for cell lines). The point to make is that for patient derived samples the answer is usually in the range of a handful of replicates per group, not one.

### 4.4 Multiple testing

About 20,000 genes are tested. At a raw p < 0.05, roughly 1,000 genes are expected to pass **by chance alone** if nothing is truly different. The dashboard therefore shows the expected number of false positives next to the cutoff, and reports Benjamini-Hochberg adjusted p-values (padj), which control the expected proportion of false discoveries among the genes called. State which threshold each figure uses.

### 4.5 Simulation: how the test choice and replication change the answer

`analysis/simulate_design_power.R` simulates negative binomial counts with 10% truly differential genes and compares a per gene t-test (the quick method in the original project blueprint) with DESeq2 from 1 vs 2 up to 10 vs 10 samples. Report, for each design:

- the false positive rate among truly null genes at p < 0.05 (should be about 5%),
- power at BH padj < 0.05, and
- the false discovery proportion among genes called.

**[RESULT: paste design_simulation.png and a 3 sentence reading.]** The expected message: with very few replicates both approaches are unreliable, DESeq2 is better calibrated because it models count variance and shares information across genes, and no test rescues a design with one sample in a group.

### 4.6 Planned sensitivity analyses

1. **Leave one early sample out.** Compare late vs early-1 and late vs early-2 separately and correlate the two sets of log2 fold changes. A low correlation means the "early" group is not internally consistent and the contrast is unreliable. **[RESULT: Spearman correlation.]**
2. **Shrunken vs unshrunken fold changes.** Do the top genes change?
3. **Threshold sweep.** Number of genes called across p and log2FC cutoffs, to show how arbitrary any single cutoff is.
4. **Pathway analysis instead of single genes.** Aggregating over many genes is more stable than ranking individual genes when n is tiny (section 5.9).

Note that a permutation test of the sample labels is not informative here: with three samples and one late label there are only three possible labelings, so the smallest attainable p-value is no smaller than about 1/3.

### 4.7 Negative and positive controls

- **Positive control:** CAF marker genes (for example *ACTA2*, *FAP*, *PDGFRB*, *COL1A1*) should be highly expressed in all three samples; their absence would suggest the wrong cell type or a mislabeled sample.
- **Negative control:** the simulation, where the truth is known, and the null genes of the unit tests, where the fraction of p < 0.05 should be near 5%.

---

## 5. Methods: tool by tool

For each tool: what it does, why it was chosen, the settings used and why, what it assumes, and how I checked it. The shell steps are `scripts/01` to `04`; the analysis is `scripts/deseq2.R`.

### 5.1 SRA Toolkit (`prefetch`, `fasterq-dump`)
- **What.** `prefetch` downloads the compressed `.sra` archive; `fasterq-dump` converts it to FASTQ using multiple threads. `--split-files` writes the two mates of a paired run to `_1` and `_2`.
- **Why.** Fetching the archive once and converting locally is more robust than streaming, and `fasterq-dump` is the maintained successor to `fastq-dump`. Temporary files go to internal disk because the conversion makes many small writes that are slow on an external drive.
- **Check.** Both mate files must have the same number of reads. Output is gzipped to a temporary name and renamed on success so a half-written file is never mistaken for a finished one.
- **Alternative.** Download FASTQ directly from the European Nucleotide Archive.

### 5.2 fastp
- **What.** One pass of adapter removal, quality filtering and QC reporting.
- **Settings.** `--detect_adapter_for_pe` finds adapters from the overlap of the two mates. `-q 20` sets the Phred score at which a base counts as "qualified" (the default is 15); by default a read is discarded if more than 40% of its bases are unqualified. **This is a filter, not a trimming threshold:** low quality tails are not cut unless sliding window trimming (`--cut_right`) is enabled.
- **Why trim at all.** HISAT2 can soft clip adapter sequence, so trimming gives a modest rather than dramatic gain; its main value is the QC report and removing very poor reads.
- **Check.** Read the HTML report: adapter content, quality by position, duplication, and the fraction of reads passing. **[RESULT: reads before and after.]**

### 5.3 HISAT2
- **What.** A splice aware aligner. Reads from mRNA span exon-exon junctions, so a plain DNA aligner would fail on them.
- **How it works.** HISAT2 uses a hierarchical graph FM index: one global index of the whole genome plus tens of thousands of small local indexes covering it in overlapping windows of roughly 60 kb. A read is anchored in the global index, and when part of it does not match, the local index near the anchor is used to find the next exon quickly. This is what lets it align spliced reads fast and with a small memory footprint.
- **Why not STAR.** STAR is excellent and often more sensitive for novel junctions, but it needs tens of gigabytes of memory for the human genome; HISAT2 needs a few. Salmon (alignment free quantification) is another strong choice for gene level counts and would be faster; I kept alignments because BAM files can be inspected in a genome browser and the alignment rate is an interpretable QC number.
- **Settings.** The prebuilt hg38 index; `--rna-strandness` is set from `STRANDEDNESS` (RF for reversely stranded paired data, the common dUTP protocols). `--summary-file` stores the alignment summary used in the dashboard's Run QC tab.
- **Assumption.** The genome build of the index and the annotation used later must match.
- **Check.** Overall alignment rate; most good human RNA-seq libraries align at roughly 85 to 95% or higher. **[RESULT per sample.]** A low rate suggests contamination, degraded RNA or the wrong genome.

### 5.4 samtools sort and index
- **What.** Sorts alignments by genomic coordinate and builds a `.bai` index for random access.
- **Settings.** Output is piped directly from HISAT2, so the large uncompressed SAM file is never written. `-m` caps memory per sorting thread (total about threads times `-m`); `-T` places temporary files on internal disk.
- **Check.** Both sides of the pipe are checked. With `pipefail` plus a copy of `PIPESTATUS`, a failure in either HISAT2 or `samtools sort` marks the sample failed, instead of silently producing a truncated BAM.

### 5.5 featureCounts
- **What.** Counts reads (or fragments) per gene by overlapping alignments with annotated exons and summing per gene (`-t exon -g gene_id`).
- **Decisions that change the numbers.**
  - **Strandedness (`-s`).** 0 ignores strand, 1 counts reads on the same strand as the gene, 2 counts reads on the opposite strand (reversely stranded libraries). The wrong value on a stranded library loses most reads; on an unstranded library any stranded setting loses about half.
  - **Paired end (`-p --countReadPairs`).** Counts each fragment once. In Subread 2.0.2 and later, `-p` alone no longer implies fragment counting.
  - **Ambiguity.** By default, reads overlapping more than one gene (`-O` off) and multi-mapping reads (`-M` off) are not counted, which is conservative.
- **Why raw integer counts.** DESeq2's model is on the count scale. TPM or FPKM remove the information the variance model needs, and are not valid input.
- **Check.** The `.summary` file: the percentage "Assigned", and the size of "Unassigned_NoFeatures" (intergenic or intronic reads, often genomic DNA or unspliced RNA) and "Unassigned_Ambiguity". Step 4 prints the assigned percentage per sample. **[RESULT.]**

### 5.6 DESeq2
- **Model.** For gene *i* and sample *j*, the count is negative binomial with mean (size factor *s_j*) times (expression level), so its variance is mean plus dispersion times mean squared. Dispersion captures biological variability beyond Poisson noise.
- **Size factors (median of ratios).** Each gene's count is divided by its geometric mean across samples; a sample's size factor is the median of those ratios. This corrects for sequencing depth and library composition, assuming most genes are not differentially expressed. It would be unreliable if a large share of genes shifted in one direction.
- **Dispersion.** A gene-wise estimate is noisy with few samples, so DESeq2 fits a smooth trend of dispersion against mean expression and shrinks each gene's estimate toward it (an empirical Bayes step). This borrowing of strength is what makes small designs workable, and it is also why the assumption that all genes follow one dispersion trend matters.
- **Test.** A Wald test of the coefficient: the log2 fold change divided by its standard error, compared to a normal distribution.
- **NA values.** `pvalue` is NA for all-zero genes and for genes flagged by Cook's distance (an outlier check that is only meaningful with several replicates). `padj` is additionally NA for genes removed by **independent filtering**, which drops genes with very low mean count that have no power, to improve the number of discoveries after adjustment.
- **Adjustment.** Benjamini-Hochberg, controlling the expected false discovery proportion.
- **Design and contrast.** `~ condition`, with `early_stage` as the reference, so a positive log2 fold change means higher in late stage. The variable of interest goes last in the formula.
- **Pre-filtering.** Genes with fewer than 10 total counts are removed for speed; this does not change results materially because independent filtering would remove them anyway.

### 5.7 apeglm shrinkage
- **What.** A Bayesian shrinkage of log2 fold changes using a heavy tailed (Cauchy) prior whose scale is estimated from the data. Noisy estimates for low count genes are pulled toward zero while large, well supported changes are left almost untouched; the reported value is the posterior mode.
- **Why.** Unshrunken fold changes of low count genes are huge and unreliable, which would dominate rankings and volcano plots.
- **Important detail.** Shrinkage changes effect sizes, not p-values. The saved table therefore contains the shrunken log2 fold change together with the original Wald p-value and padj. The dashboard colors genes by raw p-value and shrunken log2FC, which is a common heuristic and not a formal test of |log2FC| above a threshold (that would use `lfcThreshold`).
- **Note.** `lfcShrink` with `type = "apeglm"` needs a coefficient name (`resultsNames`), not a `contrast`, which is why the reference level is set explicitly.

### 5.8 Variance stabilizing transformation, PCA and sample distances
- **VST.** Raw counts have variance that grows with the mean, so distances between samples are dominated by a few highly expressed genes. `vst` makes the variance roughly constant across the mean (about log2 scale for large counts). `blind = FALSE` lets the design inform the dispersion trend. VST values are used for visualization only, never for the test.
- **PCA.** On the top 500 most variable genes, as `DESeq2::plotPCA` does; the dashboard reproduces this with `prcomp`. With three samples there are at most two informative components, so PCA is a sanity check on sample relationships, not a rich picture.
- **Check.** Do the two early samples sit closer to each other than to the late sample? If not, stage is not the main source of variation. **[RESULT.]**

### 5.9 Gene annotation and pathway analysis (org.Hs.eg.db, msigdbr, fgsea)
- **Symbols.** Ensembl IDs are mapped to symbols with `AnnotationDbi::mapIds`; version suffixes are stripped, unmapped IDs keep their Ensembl ID, and duplicated symbols are disambiguated for display while the plain symbol is kept for pathway analysis.
- **Preranked GSEA.** All genes are ranked by the Wald statistic. The statistic combines effect size and its uncertainty, so noisy genes do not dominate, which makes it a better ranking than log2 fold change alone or p-value alone. For each gene set a running sum walks down the ranked list, increasing at set members and decreasing otherwise; the enrichment score is its maximum deviation. The normalized enrichment score (NES) rescales it by the average of permutation scores of the same sign so sets of different sizes are comparable. `fgsea` computes p-values with an adaptive multilevel permutation procedure, then Benjamini-Hochberg adjusts them.
- **Gene sets.** MSigDB Hallmark (50 curated sets) plus a small custom CAF marker file. **The custom sets are illustrative and must be checked against the literature before being cited.** Sets of 5 to 500 genes are tested.
- **Assumption and caveat.** `fgsea` permutes genes, not samples, so it treats genes as independent, which they are not. This makes enrichment p-values optimistic. With three samples a sample-label permutation (the stricter alternative) is not possible. Read pathway results as hypotheses.
- **Check.** Does the sign of the enrichment match the top genes in the leading edge? **[RESULT.]**

### 5.10 R Shiny (with plotly, DT, bslib)
- **Why Shiny.** The reader's question after a volcano plot is "what does this gene look like?". A static figure cannot answer it, and a notebook requires re-running code. Shiny serves R directly to a browser, so the statistics stay in R and the reader needs only a link.
- **Reactivity.** Shiny builds a graph of dependencies: inputs feed `reactive()` expressions, which feed outputs. Reactive expressions are lazy and cached, so `filtered_res` classifies every gene once per threshold change and is shared by the volcano, MA plot, heatmap and table. `observeEvent` is for side effects (updating the dropdown when a point is clicked) and `reactiveVal` holds the shared "current gene".
- **Click to boxplot.** Each point carries its gene ID in plotly `customdata`; `event_register("plotly_click")` and `event_data("plotly_click", source = ...)` return that ID to the server, which updates the selected gene; the boxplot re-renders from the normalized count matrix.
- **Performance and hosting.** About 20,000 points are drawn with `scattergl` (WebGL). Code outside `server()` runs once per R process and is shared across sessions, so the data object is loaded once; per-session state lives inside `server()`. A slim data object (plain data frames and matrices) keeps the hosted app within free-tier memory limits.
- **Honesty built in.** A Methods and Caveats tab and a sidebar note appear whenever a group has a single sample.
- **Testing.** The pipeline output is tested on simulated data, and `shinytest2` checks that the app starts and renders without output errors.

---

## 6. Results

Fill from your run. For each item, give the figure, a one sentence reading, and what would change your mind.

| # | Result | Source in the repo | Reading prompt |
|---|---|---|---|
| 6.1 | Read and mapping QC table (reads, alignment rate, assigned %, detected genes) | Run QC tab, `results/deseq2/run_qc.csv` | Any sample more than 10 points below the others? Assigned below 50%? |
| 6.2 | PCA and sample distance heatmap | Sample QC tab | Do early samples cluster together? Is PC1 plausibly stage? |
| 6.3 | Positive control genes (CAF markers) | Gene Inspector | Expressed in all samples? |
| 6.4 | Differential expression overview: counts of up, down and total at stated cutoffs, with expected false positives | Volcano tab, p-value histogram | How does the p-value histogram look (uniform with a spike near 0, or skewed)? |
| 6.5 | MA plot | MA tab | Are large fold changes confined to low-count genes? |
| 6.6 | Top genes table (symbol, log2FC, p, padj, baseMean) | Results table | Are the top genes well expressed and consistent with the biology? |
| 6.7 | Leave one early sample out concordance | analysis (to be run) | Spearman correlation of the two fold change sets |
| 6.8 | Pathways (Hallmark and custom CAF sets) | Pathways tab | Direction, NES, leading edge genes, nominal vs adjusted |
| 6.9 | Simulation of designs | `analysis/simulate_design_power.R` | Which designs control false positives? |

---

## 7. Interpretation and limitations

Write four short paragraphs.

1. **What I would say with confidence** (technical quality, the consistency of the pipeline, the QC findings).
2. **What I would say as a hypothesis** (the top pathways and genes), each with the sentence "this needs replication in independent donors".
3. **What I would not say** (anything about mechanism, anything causal, any claim about stage specific biology generally).
4. **What would resolve it:** at least three independent donors per stage, preferably matched for culture passage and processed in one batch, and orthogonal validation (qPCR or protein) of a few candidate genes.

**Why not add normal tissue.** Bulk normal cervix is a mixture of epithelial, immune and stromal cells, so a CAF versus normal comparison would mostly measure cell composition. Merging a different study also confounds condition with batch (lab, protocol, sequencer), which cannot be corrected when every normal sample comes from the other study. A better anchor is normal cervical or other fibroblasts processed similarly, compared at the gene set level.

---

## 8. Engineering and reproducibility

- **Numbered, resumable scripts.** Each step skips finished work, writes outputs under a temporary name and renames on success, and reports failures per sample.
- **One settings file** for the shell steps (`config.sh`), one sample sheet, and one settings block in the R script.
- **Tests.** A stubbed test of the shell steps (paths with spaces, skipping, failure detection, strandedness and single-end plumbing); an end to end test of the DESeq2 step on simulated counts that checks direction, power, null calibration and the structure of the dashboard object; a dashboard smoke test; continuous integration in GitHub Actions.
- **Hardware workarounds** are documented in the README (scratch on internal disk, no intermediate SAM, capped sort memory, `caffeinate`).

**Problems I hit and how I fixed them** (a short, honest list is more convincing than a perfect story):

1. **A masked function.** `select()` failed in the dashboard because loading DESeq2 (through AnnotationDbi) masked dplyr's `select()`. Fix: call `dplyr::select` explicitly.
2. **A silent failure in a pipe.** A loop around `hisat2 | samtools sort` silently stopped after the first sample on macOS because referencing `PIPESTATUS` under `set -u` threw "unbound variable" in bash 3.2. Fix: copy `PIPESTATUS` into an array with defaults immediately and check both sides.
3. **Spaces in the project path.** A conda environment created inside a folder with a space failed in a Bioconductor data package's post-install script. Fix: create environments in the normal conda location, and quote every path in the scripts.
4. **A mis-named saved object.** The first slim dashboard file saved the sample table under the wrong variable name, so the app could not load it. The synthetic data test now checks the exact object names, so this cannot return.

---

## 9. Next steps

- Replace the illustrative CAF gene sets with published signatures and cite them.
- Add at least one independent dataset with replication, for example CAF versus normal fibroblasts, and compare at the gene set level.
- Add a donor covariate to the design if more samples become available (`~ donor + stage`).
- Add a deconvolution or marker based check that the three samples are fibroblast-like.

---

## 10. Defensibility checklist (questions I should be able to answer)

1. **Why raw counts and not TPM?** The negative binomial model needs counts; TPM discards the count scale that defines the variance.
2. **What does median of ratios assume?** That most genes are not differentially expressed; it fails when a large fraction shifts one way.
3. **What is dispersion and why shrink it?** Extra-Poisson variability, var = mean + alpha times mean squared. Gene-wise estimates are noisy at small n, so they are shrunk toward a trend.
4. **Why does DESeq2 run on 1 vs 2 and why do I still not trust it?** Dispersion is borrowed across genes, so the model is estimable, but the late group's biological variability is unobserved and donor is confounded with stage.
5. **Why shrunken log2FC for plots but Wald p-values for significance?** Shrinkage stabilizes effect sizes but does not test anything; p-values come from the unshrunken test. Cutoffs on both are a heuristic, not a test of |LFC| above a threshold.
6. **What does an NA padj mean?** Zero counts, an outlier flag, or removal by independent filtering for low mean count.
7. **Why are 1,000 genes at p < 0.05 not 1,000 discoveries?** About 5% of 20,000 null genes pass by chance; BH controls the expected proportion of false discoveries among called genes.
8. **What does BH control and what does it assume?** The expected false discovery proportion, under independence or positive dependence of tests.
9. **A sample has 40% assigned reads. What do I check first?** Strandedness (`-s`), then that the GTF matches the genome build, then intronic or genomic DNA content.
10. **Why rank by the Wald statistic for GSEA?** It combines effect size and uncertainty, so noisy genes do not dominate the ranking.
11. **What does fgsea assume?** Permuting genes treats them as independent, so p-values are optimistic; sample permutation is not possible with three samples.
12. **Why not merge in normal cervix data?** Cell composition and batch are confounded with condition.
13. **How does clicking a volcano point draw a boxplot?** `customdata` carries the gene ID; `event_data` returns it; a `reactiveVal` stores it; the boxplot renders from the normalized counts for that gene.
14. **Reactive expression versus observer?** A reactive expression is lazy, cached and returns a value; an observer runs eagerly for side effects.
15. **What would you do with 3 samples per group?** Test directly with `~ condition` (adding donor or batch covariates if available), use padj for calls, and run sample-label permutations for pathways.

---

## 11. References

Verify each against the publisher page before submitting.

- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35:2084-2092. (apeglm)
- Anders S, Huber W (2010). Differential expression analysis for sequence count data. *Genome Biology* 11:R106. (variance stabilization)
- Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology* 37:907-915.
- Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30:923-930.
- Chen S, Zhou Y, Chen Y, Gu J (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34:i884-i890.
- Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* 57:289-300.
- Subramanian A, et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545-15550.
- Korotkevich G, et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 060012. (fgsea)
- Liberzon A, et al. (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems* 1:417-425.
- Dobin A, et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29:15-21. (alternative aligner)
- Patro R, et al. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417-419. (alternative quantification)
- Chang W, et al. Shiny: web application framework for R. R package.
- Okabe M, Ito K (2008). Color Universal Design (CUD): a colorblind-safe palette.

---

## 12. Figure and table checklist

- [ ] Fig 1. Pipeline diagram (the Mermaid chart in the README, exported).
- [ ] Fig 2. fastp before and after quality, or the Run QC table.
- [ ] Fig 3. PCA and sample distance heatmap.
- [ ] Fig 4. Volcano plot with the cutoffs stated in the caption, plus the expected false positives.
- [ ] Fig 5. MA plot.
- [ ] Fig 6. p-value histogram.
- [ ] Fig 7. One example gene boxplot from the click-to-inspect panel, with the caption "the late stage group is a single sample".
- [ ] Fig 8. Pathway bar plot.
- [ ] Fig 9. Simulation of designs (`design_simulation.png`).
- [ ] Table 1. Samples and metadata (section 3).
- [ ] Table 2. Alignment and assignment QC.
- [ ] Table 3. Top differential genes.
- [ ] Dashboard screenshot and a short GIF of the click-to-boxplot interaction, saved in `docs/img/`.
