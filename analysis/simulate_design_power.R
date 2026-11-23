#!/usr/bin/env Rscript
# Why replication matters: simulate RNA-seq counts with KNOWN differentially expressed genes, then compare
#   (a) a per-gene Student t-test on log2(normalized counts + 1)   (the blueprint's "quick" test)
#   (b) DESeq2 (negative binomial model with shared dispersion trend)
# across designs from 1 vs 2 (this project) up to 10 vs 10, and measure what actually matters:
#   false positive rate at p < 0.05 among truly null genes,
#   power (fraction of true DE genes found at BH-adjusted p < 0.05),
#   false discovery proportion among genes called at BH-adjusted p < 0.05.
# Run from the repo root:  Rscript analysis/simulate_design_power.R
# Environment: SIM_REPS (default 10 simulations per design), SIM_OUTDIR (default results/design_simulation)

suppressPackageStartupMessages({ library(DESeq2); library(tidyverse) })
reps   <- as.integer(Sys.getenv("SIM_REPS", "10"))
outdir <- Sys.getenv("SIM_OUTDIR", "results/design_simulation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

n_genes <- 2000; de_frac <- 0.10; true_lfc <- 1            # 10% of genes, 2-fold change
designs <- data.frame(n_early = c(2, 3, 5, 10), n_late = c(1, 3, 5, 10))

simulate_counts <- function(n_early, n_late) {
  n <- n_early + n_late
  group <- factor(rep(c("early", "late"), c(n_early, n_late)), levels = c("early", "late"))
  mu    <- exp(rnorm(n_genes, mean = 5, sd = 1.5))             # gene means
  alpha <- 0.1 * exp(rnorm(n_genes, 0, 0.5))                   # gene dispersions (CV about 0.3 plus Poisson)
  sf    <- exp(rnorm(n, 0, 0.2))                               # library size differences
  is_de <- seq_len(n_genes) %in% sample(n_genes, round(n_genes * de_frac))
  fold  <- ifelse(is_de, 2^(sample(c(-1, 1), n_genes, TRUE) * true_lfc), 1)
  counts <- vapply(seq_len(n), function(j) {
    m <- mu * sf[j] * (if (group[j] == "late") fold else 1)
    as.numeric(rnbinom(n_genes, mu = m, size = 1 / alpha))
  }, numeric(n_genes))
  rownames(counts) <- sprintf("g%04d", seq_len(n_genes)); colnames(counts) <- paste0("s", seq_len(n))
  list(counts = counts, group = group, is_de = is_de)
}

ttest_p <- function(counts, group) {
  sf <- DESeq2::estimateSizeFactorsForMatrix(counts)
  x  <- log2(sweep(counts, 2, sf, "/") + 1)
  a <- group == "early"; b <- group == "late"
  apply(x, 1, function(v) tryCatch(t.test(v[b], v[a], var.equal = TRUE)$p.value, error = function(e) NA_real_))
}
deseq_p <- function(counts, group) {
  dds <- DESeqDataSetFromMatrix(counts, data.frame(group = group, row.names = colnames(counts)), ~ group)
  dds <- suppressMessages(suppressWarnings(DESeq(dds, quiet = TRUE)))
  results(dds)$pvalue
}
score <- function(p, is_de, method, design, rep) {
  padj <- p.adjust(p, "BH")
  called <- !is.na(padj) & padj < 0.05
  tibble(method = method, design = design, rep = rep,
         fpr   = mean(p[!is_de] < 0.05, na.rm = TRUE),
         power = mean(called[is_de]),
         fdp   = if (sum(called) > 0) mean(!is_de[called]) else 0)
}

out <- list()
for (d in seq_len(nrow(designs))) {
  label <- sprintf("%d early vs %d late", designs$n_early[d], designs$n_late[d])
  for (r in seq_len(reps)) {
    set.seed(1000 * d + r)
    s <- simulate_counts(designs$n_early[d], designs$n_late[d])
    out[[length(out) + 1]] <- score(ttest_p(s$counts, s$group), s$is_de, "t-test", label, r)
    out[[length(out) + 1]] <- score(deseq_p(s$counts, s$group), s$is_de, "DESeq2", label, r)
  }
  message("finished design: ", label)
}
res <- bind_rows(out)
res$design <- factor(res$design, levels = unique(res$design))
summary_tbl <- res %>% group_by(design, method) %>%
  summarise(across(c(fpr, power, fdp), mean, .names = "mean_{.col}"), .groups = "drop")
write.csv(res, file.path(outdir, "simulation_replicates.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(outdir, "simulation_summary.csv"), row.names = FALSE)
print(summary_tbl, n = Inf)

long <- res %>% pivot_longer(c(fpr, power, fdp), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, fpr = "False positive rate (p < 0.05, null genes)",
                         power = "Power (BH padj < 0.05, true DE genes)",
                         fdp = "False discovery proportion (BH padj < 0.05)"))
p <- ggplot(long, aes(design, value, fill = method)) +
  geom_boxplot(outlier.size = 0.8) + facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("DESeq2" = "#0072B2", "t-test" = "#E69F00")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Replication, not the choice of test, drives what you can conclude",
       subtitle = paste0(reps, " simulations per design; ", n_genes, " genes, ", 100 * de_frac, "% truly DE (2-fold)")) +
  theme_minimal(base_size = 12) + theme(legend.position = "top", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(outdir, "design_simulation.png"), p, width = 7.5, height = 9, dpi = 150)
cat("Saved results to", outdir, "\n")
