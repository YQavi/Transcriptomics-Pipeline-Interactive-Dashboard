#!/usr/bin/env Rscript
# Step 5: DESeq2 differential expression + fgsea pathways + run QC + dashboard data.
# Run from the repository root (scripts/05_run_deseq2.sh does this for you), or directly:
#   COUNTS_FILE=/path/gene_counts.txt BAM_DIR=/path/aligned_bam Rscript scripts/deseq2.R
# Optional environment overrides (used by the automated tests): RNASEQ_SAMPLES, RNASEQ_OUTDIR,
# RNASEQ_APP_COPY, RNASEQ_ORGDB, RNASEQ_MSIGDB_COLLECTION, RNASEQ_CUSTOM_SETS.

suppressPackageStartupMessages({
  library(DESeq2); library(tidyverse)
})

# ============================== SETTINGS (edit here) ==============================
cfg <- list(
  dataset_title   = "CAFs, late vs early stage FIGO cervical cancer (GSE347724)",
  samples         = Sys.getenv("RNASEQ_SAMPLES", "config/samples.tsv"),          # columns: sample, condition, [covariates]
  counts          = Sys.getenv("COUNTS_FILE", "../gene_counts.txt"),   # raw featureCounts output
  hisat_dir       = Sys.getenv("BAM_DIR", "../aligned_bam"),           # holds <sample>_summary.txt
  outdir          = Sys.getenv("RNASEQ_OUTDIR", "results/deseq2"),
  app_data_copy   = Sys.getenv("RNASEQ_APP_COPY", "app/data/demo_deseq2_data.RData"),  # slim file the dashboard loads when deployed
  design          = "~ condition",                 # e.g. "~ batch + condition"; variable of interest LAST
  contrast        = c("condition", "late_stage", "early_stage"),   # variable, numerator, denominator
  min_total_count = 10,
  candidate_pvalue = 0.05,                         # only for the exported candidates table / static volcano
  candidate_lfc    = 1.5,
  orgdb           = Sys.getenv("RNASEQ_ORGDB", "org.Hs.eg.db"),                # "org.Mm.eg.db" for mouse; "" to skip symbol mapping
  gene_id_type    = "ENSEMBL",                     # ENSEMBL | ENTREZID | SYMBOL (SYMBOL = already symbols)
  msigdbr_species    = "Homo sapiens",
  msigdbr_collection = Sys.getenv("RNASEQ_MSIGDB_COLLECTION", "H"),                        # Hallmark; "" to skip MSigDB
  custom_gene_sets   = Sys.getenv("RNASEQ_CUSTOM_SETS", "config/gene_sets/caf_signatures.gmt")   # "" to skip
)
# ==================================================================================

cat("==========================================\nDESeq2 pipeline\n==========================================\n")
dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)

# Find the columns of a counts-like table that belong to each sample (BAM paths contain the sample ID).
# Returns column indices in sample order, or NULL if the names do not match uniquely.
match_cols <- function(ids, cols) {
  idx <- vapply(ids, function(x) { h <- grep(x, cols, fixed = TRUE); if (length(h) == 1) h else NA_integer_ }, integer(1))
  if (anyNA(idx)) NULL else idx
}

# ---- 1. Counts and metadata ----
if (!file.exists(cfg$counts)) stop("Counts file not found: ", cfg$counts, " (run step 04 or set COUNTS_FILE)")
sheet <- read.delim(cfg$samples, stringsAsFactors = FALSE, check.names = FALSE)
skip_n <- if (startsWith(readLines(cfg$counts, n = 1), "#")) 1L else 0L
raw <- read.table(cfg$counts, header = TRUE, skip = skip_n, row.names = 1,
                  check.names = FALSE, comment.char = "", quote = "")
counts <- raw[, 6:ncol(raw), drop = FALSE]    # columns 1-5 after the ID are Chr/Start/End/Strand/Length

idx <- match_cols(sheet$sample, colnames(counts))
if (!is.null(idx)) {
  counts <- counts[, idx, drop = FALSE]
} else {
  message("Could not match count columns to sample IDs by name; assuming samples.tsv order.")
  stopifnot("number of count columns != number of samples" = ncol(counts) == nrow(sheet))
}
colnames(counts) <- sheet$sample

colData <- sheet; rownames(colData) <- colData$sample
vars <- all.vars(as.formula(cfg$design))
miss <- setdiff(vars, colnames(colData))
if (length(miss)) stop("Design variables not in samples sheet: ", paste(miss, collapse = ", "))
colData[vars] <- lapply(colData[vars], factor)

grp <- cfg$contrast[1]; num <- cfg$contrast[2]; den <- cfg$contrast[3]
if (!all(c(num, den) %in% levels(colData[[grp]]))) stop("Contrast levels not found in '", grp, "'")
colData[[grp]] <- relevel(colData[[grp]], ref = den)    # coefficient becomes <grp>_<num>_vs_<den>

n_per <- table(colData[[grp]])
if (any(n_per < 2)) warning("A group has <2 replicates (", paste(names(n_per), n_per, collapse = ", "),
                            "): p-values are exploratory only.")

# ---- 2. DESeq2 ----
dds <- DESeqDataSetFromMatrix(counts[, rownames(colData)], colData, design = as.formula(cfg$design))
keep <- rowSums(counts(dds)) >= cfg$min_total_count
dds <- dds[keep, ]
cat("Genes passing low-count filter:", sum(keep), "\n")
dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    dds <- estimateSizeFactors(dds, type = "poscounts")
    DESeq(dds, fitType = "mean")
  }
)

res  <- results(dds, contrast = c(grp, num, den))
coef <- paste0(grp, "_", num, "_vs_", den)
stopifnot(coef %in% resultsNames(dds))
shr  <- tryCatch({
  suppressPackageStartupMessages(library(apeglm))
  lfcShrink(dds, coef = coef, type = "apeglm")   # stable effect sizes for ranking and plotting
}, error = function(e) {
  lfcShrink(dds, coef = coef, type = "normal")
})

# Shrunken LFC + the ORIGINAL Wald p-values (apeglm does not return p-values)
res_df <- as.data.frame(shr)
res_df$pvalue <- res[rownames(res_df), "pvalue"]
res_df$padj   <- res[rownames(res_df), "padj"]
res_df$stat   <- res[rownames(res_df), "stat"]          # Wald statistic, used to rank genes for GSEA
res_df <- res_df %>% tibble::rownames_to_column("gene_id") %>% dplyr::arrange(pvalue)

# ---- 3. Gene symbols ----
symbol <- res_df$gene_id
if (nzchar(cfg$orgdb) && toupper(cfg$gene_id_type) != "SYMBOL") {
  symbol <- tryCatch({
    suppressPackageStartupMessages(library(cfg$orgdb, character.only = TRUE))
    keys <- if (toupper(cfg$gene_id_type) == "ENSEMBL") sub("\\..*$", "", res_df$gene_id) else res_df$gene_id
    sym <- unname(suppressMessages(AnnotationDbi::mapIds(get(cfg$orgdb), keys = keys,
            keytype = toupper(cfg$gene_id_type), column = "SYMBOL", multiVals = "first")))
    ifelse(is.na(sym) | sym == "", res_df$gene_id, sym)
  }, error = function(e) {
    message("Symbol mapping skipped: ", conditionMessage(e))
    res_df$gene_id
  })
}
gene_name <- symbol
dup <- duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE)
gene_name[dup] <- paste0(gene_name[dup], " (", res_df$gene_id[dup], ")")
res_df$symbol    <- symbol        # plain symbol, used for pathway analysis
res_df$gene_name <- gene_name     # unique display label, used in the dashboard
res_df <- res_df[, c("gene_id", "gene_name", "symbol",
                     setdiff(names(res_df), c("gene_id", "gene_name", "symbol")))]

# ---- 4. Pathway analysis (fgsea). A failure here never stops the run. ----
run_gsea <- function() {
  suppressPackageStartupMessages(library(fgsea))
  ranks_df <- res_df %>%
    dplyr::filter(!is.na(stat), !is.na(symbol), symbol != "") %>%
    dplyr::group_by(symbol) %>%
    dplyr::slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  ranks <- setNames(ranks_df$stat, ranks_df$symbol)

  sets <- list(); custom <- list()
  if (nzchar(cfg$msigdbr_collection)) {
    suppressPackageStartupMessages(library(msigdbr))
    m <- tryCatch(msigdbr(species = cfg$msigdbr_species, collection = cfg$msigdbr_collection),
                  error = function(e) msigdbr(species = cfg$msigdbr_species, category = cfg$msigdbr_collection))
    sets <- split(m$gene_symbol, m$gs_name)
  }
  if (nzchar(cfg$custom_gene_sets) && file.exists(cfg$custom_gene_sets)) {
    custom <- gmtPathways(cfg$custom_gene_sets)
    sets <- c(sets, custom)
  }
  if (!length(sets)) stop("no gene sets configured")

  set.seed(1)
  out <- as.data.frame(fgsea(pathways = sets, stats = ranks, minSize = 5, maxSize = 500))
  out$leadingEdge <- vapply(out$leadingEdge, paste, character(1), collapse = ", ")
  out$source <- ifelse(out$pathway %in% names(custom), "Custom", "MSigDB")
  out[order(out$pval), ]
}
gsea_df <- tryCatch(run_gsea(), error = function(e) {
  message("GSEA skipped: ", conditionMessage(e)); data.frame()
})

# ---- 5. Run QC (HISAT2 alignment rate, featureCounts assignment rate); optional inputs ----
parse_hisat <- function(f) {
  if (!file.exists(f)) return(data.frame(total_reads = NA_real_, alignment_rate = NA_real_))
  l <- readLines(f)
  data.frame(
    total_reads    = suppressWarnings(as.numeric(sub(" .*", "", trimws(l[grep("reads; of these", l)[1]])))),
    alignment_rate = suppressWarnings(as.numeric(sub("%.*", "", grep("overall alignment rate", l, value = TRUE)[1]))))
}
qc_df <- tryCatch({
  hs <- dplyr::bind_rows(lapply(file.path(cfg$hisat_dir, paste0(sheet$sample, "_summary.txt")), parse_hisat))
  out <- data.frame(sample = sheet$sample, group = colData[sheet$sample, grp], hs,
                    detected_genes = colSums(counts[, sheet$sample, drop = FALSE] > 0), row.names = NULL)
  fc_file <- paste0(cfg$counts, ".summary")
  if (file.exists(fc_file)) {
    fcs <- read.delim(fc_file, row.names = 1, check.names = FALSE)
    ci  <- match_cols(sheet$sample, colnames(fcs))
    if (!is.null(ci)) {
      fcs <- fcs[, ci, drop = FALSE]
      out$assigned_reads <- as.numeric(fcs["Assigned", ])
      out$assigned_pct   <- round(100 * as.numeric(fcs["Assigned", ]) / colSums(fcs), 1)
    }
  }
  if (all(is.na(out$alignment_rate)) && is.null(out$assigned_pct)) NULL else out
}, error = function(e) { message("QC table skipped: ", conditionMessage(e)); NULL })
if (is.null(qc_df)) message("No run QC inputs found; the dashboard Run QC tab will be empty.")

# ---- 6. Outputs ----
norm_counts <- counts(dds, normalized = TRUE)
vsd     <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)
)
vst_mat <- assay(vsd)
out <- function(f) file.path(cfg$outdir, f)

write.csv(res_df, out("deseq2_all_results.csv"), row.names = FALSE)
write.csv(norm_counts, out("normalized_counts.csv"))
write.csv(vst_mat, out("vst_counts.csv"))
write.csv(gsea_df, out("gsea_results.csv"), row.names = FALSE)
if (!is.null(qc_df)) write.csv(qc_df, out("run_qc.csv"), row.names = FALSE)

top <- res_df %>% dplyr::filter(!is.na(pvalue), pvalue < cfg$candidate_pvalue,
                                abs(log2FoldChange) >= cfg$candidate_lfc)
write.csv(top, out("top_candidate_genes.csv"), row.names = FALSE)

contrast_info <- list(variable = grp, numerator = num, denominator = den,
                      design = cfg$design, title = cfg$dataset_title)

# Full object for your own downstream work
save(dds, res_df, norm_counts, vsd, colData, contrast_info, file = out("deseq2_data.RData"))

# Slim object for the dashboard: plain data frames and matrices only (no DESeq2 classes, no local paths).
# Built in a scratch environment so the saved object names are exactly what app.R loads.
app_env <- new.env()
assign("res_df",        res_df,        envir = app_env)
assign("norm_counts",   norm_counts,   envir = app_env)
assign("vst_mat",       vst_mat,       envir = app_env)
assign("gsea_df",       gsea_df,       envir = app_env)
assign("qc_df",         qc_df,         envir = app_env)
assign("contrast_info", contrast_info, envir = app_env)
assign("colData",       colData[, intersect(c("sample", vars), colnames(colData)), drop = FALSE], envir = app_env)
app_objs <- c("res_df", "norm_counts", "vst_mat", "gsea_df", "qc_df", "contrast_info", "colData")
save(list = app_objs, envir = app_env, file = out("app_data.RData"))
if (nzchar(cfg$app_data_copy)) {
  dir.create(dirname(cfg$app_data_copy), showWarnings = FALSE, recursive = TRUE)
  file.copy(out("app_data.RData"), cfg$app_data_copy, overwrite = TRUE)
}

# ---- 7. Static QC plots ----
tryCatch({
  pdf(out("qc_pca_volcano.pdf"), width = 10, height = 5)
  on.exit(dev.off(), add = TRUE)
  pca <- plotPCA(vsd, intgroup = grp, returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))
  print(ggplot(pca, aes(PC1, PC2, color = .data[[grp]], label = name)) +
    geom_point(size = 4) + geom_text(vjust = -1, size = 3) +
    labs(x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance"),
         title = cfg$dataset_title) + theme_bw())
  if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    print(EnhancedVolcano::EnhancedVolcano(res_df, lab = res_df$gene_name, x = "log2FoldChange", y = "pvalue",
      pCutoff = cfg$candidate_pvalue, FCcutoff = cfg$candidate_lfc, pointSize = 2, labSize = 3,
      title = paste(num, "vs", den),
      subtitle = paste0("p < ", cfg$candidate_pvalue, " & |log2FC| >= ", cfg$candidate_lfc)))
  }
}, error = function(e) message("Static QC plots skipped: ", conditionMessage(e)))

cat("==========================================\n",
    "Done. Outputs in '", cfg$outdir, "/'. Dashboard data: ", out("app_data.RData"), "\n",
    "==========================================\n", sep = "")