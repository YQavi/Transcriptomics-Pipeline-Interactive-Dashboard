# End-to-end tests of scripts/deseq2.R on simulated featureCounts output.
# Needs DESeq2, apeglm, EnhancedVolcano, tidyverse (pathway analysis is switched off for speed).
library(testthat)

find_repo_root <- function() {
  d <- normalizePath(getwd())
  for (i in 1:6) {
    if (file.exists(file.path(d, "scripts", "deseq2.R"))) return(d)
    d <- dirname(d)
  }
  stop("Could not locate the repository root from ", getwd())
}
repo <- find_repo_root()

# Simulated project: negative binomial counts, genes 1-20 up and 21-40 down in late_stage
simulate_project <- function(dir, n_early = 3, n_late = 3, n_genes = 600, seed = 42,
                             conditions = c("early_stage", "late_stage"), extra_sample = FALSE) {
  set.seed(seed)
  early <- paste0("SRR_E", seq_len(n_early)); late <- paste0("SRR_L", seq_len(n_late))
  samples <- c(early, late)
  cond <- c(rep(conditions[1], n_early), rep(conditions[2], n_late))
  mu <- exp(rnorm(n_genes, mean = 6, sd = 1.5))
  idx <- seq_len(n_genes)
  fold <- ifelse(idx <= 20, 4, ifelse(idx <= 40, 0.25, 1))
  cnt <- vapply(seq_along(samples), function(j) {
    m <- if (cond[j] == conditions[2]) mu * fold else mu
    as.numeric(rnbinom(n_genes, mu = m, size = 20))
  }, numeric(n_genes))

  bam <- file.path(dir, "aligned_bam", paste0(samples, ".sorted.bam"))
  tab <- data.frame(Geneid = sprintf("GENE%04d", idx), Chr = "chr1", Start = 1L, End = 100L,
                    Strand = "+", Length = 100L, cnt, check.names = FALSE)
  colnames(tab)[7:ncol(tab)] <- bam
  counts_file <- file.path(dir, "gene_counts.txt")
  writeLines("# Program:featureCounts v2.0.6; Command: simulated for tests", counts_file)
  suppressWarnings(write.table(tab, counts_file, sep = "\t", quote = FALSE,
                               row.names = FALSE, col.names = TRUE, append = TRUE))

  sheet_samples <- if (extra_sample) c(samples, "SRR_X9") else samples
  sheet_cond    <- if (extra_sample) c(cond, conditions[2]) else cond
  sheet <- data.frame(sample = sheet_samples, condition = sheet_cond, accession = sheet_samples)
  samples_file <- file.path(dir, "samples.tsv")
  write.table(sheet, samples_file, sep = "\t", quote = FALSE, row.names = FALSE)
  list(counts = counts_file, samples = samples_file, outdir = file.path(dir, "out"), dir = dir)
}

run_deseq2 <- function(p) {
  env <- c(COUNTS_FILE = p$counts, BAM_DIR = file.path(p$dir, "aligned_bam"),
           RNASEQ_SAMPLES = p$samples, RNASEQ_OUTDIR = p$outdir,
           RNASEQ_APP_COPY = "", RNASEQ_ORGDB = "", RNASEQ_MSIGDB_COLLECTION = "", RNASEQ_CUSTOM_SETS = "")
  old <- setwd(repo); on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "scripts/deseq2.R",
                                  env = paste0(names(env), "=", env), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  list(status = status, log = out)
}
tail_log <- function(r) paste(utils::tail(r$log, 25), collapse = "\n")

test_that("3 vs 3 run completes and writes every output", {
  p <- simulate_project(withr::local_tempdir())
  r <- run_deseq2(p)
  expect_equal(r$status, 0L, info = tail_log(r))
  for (f in c("app_data.RData", "deseq2_data.RData", "deseq2_all_results.csv", "normalized_counts.csv",
              "vst_counts.csv", "top_candidate_genes.csv", "gsea_results.csv", "qc_pca_volcano.pdf")) {
    expect_true(file.exists(file.path(p$outdir, f)), info = f)
  }
})

test_that("the dashboard object has exactly the names app.R loads, with consistent samples", {
  # Regression test: an earlier version saved the sample table under the wrong object name.
  p <- simulate_project(withr::local_tempdir())
  r <- run_deseq2(p); expect_equal(r$status, 0L, info = tail_log(r))
  e <- new.env(); load(file.path(p$outdir, "app_data.RData"), envir = e)
  expect_setequal(ls(e), c("res_df", "norm_counts", "vst_mat", "gsea_df", "qc_df", "contrast_info", "colData"))
  expect_setequal(rownames(e$colData), colnames(e$norm_counts))
  expect_setequal(rownames(e$colData), colnames(e$vst_mat))
  expect_true(all(c("gene_id", "gene_name", "symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "stat")
                  %in% names(e$res_df)))
  expect_equal(e$contrast_info$numerator, "late_stage")
  expect_false(is.unsorted(e$res_df$pvalue, na.rm = TRUE))      # results are ranked by p-value
})

test_that("the statistics behave: direction, power and null calibration", {
  p <- simulate_project(withr::local_tempdir())
  r <- run_deseq2(p); expect_equal(r$status, 0L, info = tail_log(r))
  e <- new.env(); load(file.path(p$outdir, "app_data.RData"), envir = e)
  res <- e$res_df
  up   <- res[res$gene_id %in% sprintf("GENE%04d", 1:20), ]
  down <- res[res$gene_id %in% sprintf("GENE%04d", 21:40), ]
  null <- res[!res$gene_id %in% sprintf("GENE%04d", 1:40), ]
  expect_gt(median(up$log2FoldChange), 1)                       # 4-fold up = 2 log2 units, shrunk a little
  expect_lt(median(down$log2FoldChange), -1)
  expect_gt(mean(up$padj < 0.05, na.rm = TRUE), 0.7)            # power on strong effects
  expect_lt(mean(null$pvalue < 0.05, na.rm = TRUE), 0.12)       # about 5% expected under the null
})

test_that("the real project layout (1 late vs 2 early) runs and warns that p-values are exploratory", {
  p <- simulate_project(withr::local_tempdir(), n_early = 2, n_late = 1)
  r <- run_deseq2(p)
  expect_equal(r$status, 0L, info = tail_log(r))
  expect_true(any(grepl("exploratory", r$log)))
  expect_true(file.exists(file.path(p$outdir, "app_data.RData")))
})

test_that("a sample missing from the counts table is a hard error", {
  p <- simulate_project(withr::local_tempdir(), extra_sample = TRUE)
  r <- run_deseq2(p)
  expect_false(r$status == 0L)
})

test_that("contrast levels that are not in the sample sheet are a hard error", {
  p <- simulate_project(withr::local_tempdir(), conditions = c("control", "treated"))
  r <- run_deseq2(p)
  expect_false(r$status == 0L)
  expect_true(any(grepl("Contrast levels", r$log)))
})
