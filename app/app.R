# Interactive bulk RNA-seq dashboard.
# Reads the slim object written by the pipeline (app_data.RData): no DESeq2 needed at runtime,
# so it deploys to shinyapps.io / Posit Connect with a light package footprint.
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
})

# ---- Step 1: Ingest + validate -------------------------------------------
# Order of lookup: DESEQ2_DATA env var, pipeline output, bundled demo object
data_path <- Sys.getenv("DESEQ2_DATA", unset = "")
if (!nzchar(data_path)) {
  candidates <- c("../results/deseq2/app_data.RData", "results/deseq2/app_data.RData",
                  "data/demo_deseq2_data.RData", "app/data/demo_deseq2_data.RData")
  data_path <- candidates[file.exists(candidates)][1]
}
if (is.na(data_path) || !file.exists(data_path)) {
  stop("Cannot find app_data.RData. Run the pipeline first (or set DESEQ2_DATA).")
}
load(data_path)  # res_df, norm_counts, vst_mat, colData, contrast_info, gsea_df, qc_df
if (!exists("vst_mat")) stop("This file lacks 'vst_mat'. Re-run the pipeline to create app_data.RData.")
if (!exists("gsea_df")) gsea_df <- data.frame()
if (!exists("qc_df"))   qc_df   <- NULL
if (!exists("contrast_info")) {
  contrast_info <- list(variable = "condition", numerator = "late_stage",
                        denominator = "early_stage", title = "Bulk RNA-seq")
}

# Sample IDs in the matrices must exactly match the metadata rows
if (!setequal(colnames(norm_counts), rownames(colData)) ||
    !setequal(colnames(vst_mat), rownames(colData))) {
  stop("Sample IDs in the expression matrices do not match the metadata sheet.")
}
colData  <- colData[colnames(norm_counts), , drop = FALSE]
vst_mat  <- vst_mat[, rownames(colData), drop = FALSE]
res_df   <- as.data.frame(res_df)
if (!"gene_name" %in% names(res_df)) res_df$gene_name <- res_df$gene_id

grp_var    <- contrast_info$variable
grp_levels <- levels(factor(colData[[grp_var]]))
n_per      <- table(colData[[grp_var]])
cmp_label  <- paste(contrast_info$numerator, "vs", contrast_info$denominator)

# Dropdown shows gene name, value is the gene_id used to index the matrices
gene_choices <- setNames(res_df$gene_id, res_df$gene_name)
gene_label   <- function(id) res_df$gene_name[match(id, res_df$gene_id)]

# Okabe-Ito: colorblind-safe
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
cond_colors <- setNames(okabe[seq_along(grp_levels)], grp_levels)
sig_colors  <- c("Up-regulated" = "#D55E00", "Down-regulated" = "#0072B2", "Not Significant" = "grey75")

# Boxplot for one gene across all samples
make_gene_plot <- function(gene) {
  df <- data.frame(
    Sample    = rownames(colData),
    Value     = as.numeric(norm_counts[gene, rownames(colData)]),
    Condition = colData[[grp_var]]
  )
  ggplot(df, aes(x = Condition, y = Value, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    geom_jitter(width = 0.08, height = 0, size = 4, shape = 21, color = "black") +
    scale_fill_manual(values = cond_colors) +
    labs(title = paste("Gene Expression:", gene_label(gene)),
         y = "Normalized Count", x = NULL) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
}

# ---- UI -------------------------------------------------------------------
ui <- page_navbar(
  title = contrast_info$title,
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    title = "Threshold Controls",
    sliderInput("lfc_thresh", "Log2 Fold Change Cutoff:",
                min = 0, max = 5, value = 2, step = 0.1),
    # p-value cutoff moves in increments of 0.001
    sliderInput("pval_thresh", "p-value Cutoff:",
                min = 0.001, max = 0.1, value = 0.05, step = 0.001, sep = ""),
    sliderInput("top_labels", "Label top genes on plots:",
                min = 0, max = 30, value = 10, step = 1),
    selectizeInput("highlight", "Highlight genes:", choices = NULL, multiple = TRUE,
                   options = list(placeholder = "Type gene names...")),
    hr(),
    helpText(paste0("Comparison: ", cmp_label, " (n = ",
                    paste(names(n_per), n_per, sep = ": ", collapse = ", "), ")")),
    uiOutput("fpNote"),
    if (min(n_per) < 2) helpText("Note: a group has only one sample, so p-values are exploratory.")
  ),

  nav_panel("Sample QC & PCA",
    layout_columns(
      card(card_header("Principal Component Analysis (PCA)"),
           plotlyOutput("pcaPlot", height = "450px")),
      card(card_header("Sample-to-Sample Distance Matrix"),
           plotOutput("sampleDistPlot", height = "450px"))
    )
  ),

  nav_panel("Run QC",
    uiOutput("qcNote"),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Per-sample metrics"), DTOutput("qcTable")),
      card(card_header("Alignment and assignment rates"),
           plotlyOutput("qcPlot", height = "380px"))
    )
  ),

  nav_panel("Volcano Plot",
    layout_columns(
      col_widths = c(8, 4),
      card(card_header("Differential Gene Expression (click a point)"),
           plotlyOutput("volcanoPlot", height = "560px")),
      card(card_header("Selected gene by group"),
           plotOutput("volcanoBox", height = "520px"))
    )
  ),

  nav_panel("MA Plot",
    layout_columns(
      col_widths = c(8, 4),
      card(card_header("Mean expression vs fold change (click a point)"),
           plotlyOutput("maPlot", height = "560px")),
      card(card_header("Selected gene by group"),
           plotOutput("maBox", height = "520px"))
    )
  ),

  nav_panel("P-value Histogram",
    card(card_header("Distribution of raw p-values across all tested genes"),
         plotOutput("pvalHist", height = "420px"),
         helpText("How to read it: if nothing were truly different the bars would be flat. A flat floor with a ",
                  "spike near 0 means real signal on top of the null. A hill in the middle or a pile near 1 ",
                  "suggests the test is conservative or the model does not fit. The dashed line is the current cutoff."))
  ),

  nav_panel("Heatmap",
    card(card_header("Top differential genes (row z-scored VST expression)"),
         sliderInput("heat_n", "Number of genes:", min = 5, max = 100, value = 30, step = 5),
         textOutput("heatNote"),
         plotOutput("heatmapPlot", height = "650px"))
  ),

  nav_panel("Gene Inspector",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Gene Selection"),
           selectizeInput("selectedGene", "Search Gene Name:", choices = NULL,
                          options = list(placeholder = "Type to search gene...")),
           tableOutput("geneStatsTable"),
           downloadButton("dl_gene_png", "Download plot (PNG)")),
      card(card_header("Normalized Read Counts Across Samples"),
           plotOutput("geneExprPlot", height = "450px"))
    )
  ),

  nav_panel("Pathways (GSEA)",
    card(card_header("Gene set enrichment on the Wald statistic (fgsea)"),
         uiOutput("gseaNote"),
         layout_columns(
           col_widths = c(4, 4, 4),
           selectInput("gsea_src", "Gene set source:", c("All", "MSigDB", "Custom")),
           sliderInput("gsea_top", "Sets per direction in plot:", min = 3, max = 20, value = 8),
           div(style = "margin-top: 28px;", downloadButton("dl_gsea", "Download CSV"))
         ),
         plotOutput("gseaPlot", height = "420px"),
         DTOutput("gseaTable"))
  ),

  nav_panel("Results Table",
    card(card_header("Differential Expression Results"),
         checkboxGroupInput("sig_filter", "Show genes:",
                            choices  = c("Up-regulated", "Down-regulated", "Not Significant"),
                            selected = c("Up-regulated", "Down-regulated", "Not Significant"),
                            inline = TRUE),
         downloadButton("dl_table", "Download filtered table (CSV)"),
         DTOutput("resultsTable"))
  ),

  nav_panel("Methods & Caveats",
    card(card_header("How to read this dashboard"),
      tags$h5("Pipeline"),
      p("Reads were trimmed with fastp, aligned with HISAT2, counted with featureCounts and tested with ",
        "DESeq2 (negative binomial model). Log2 fold changes are shrunken with apeglm; p-values are the ",
        "original Wald p-values. PCA, distances and heatmaps use variance stabilized (VST) values."),
      tags$h5("Sample size"),
      p(paste0("Design: ", paste(names(n_per), n_per, sep = " = ", collapse = ", "), " samples. ")),
      if (min(n_per) < 2)
        p(strong("A group has a single sample."),
          "Within-group variability cannot be estimated for it, so p-values are optimistic and should be ",
          "used to rank genes, not to claim significance. Group is probably also confounded with the ",
          "individual donor."),
      tags$h5("Thresholds"),
      p("Up and down calls use the raw p-value and log2 fold change cutoffs in the sidebar. ",
        "Adjusted p-values are shown in the tables but not used for coloring."),
      tags$h5("Pathways"),
      p("fgsea ranks genes by the Wald statistic. Custom CAF sets are illustrative marker lists and ",
        "should be checked against primary literature before being cited."),
      tags$h5("Not included"),
      p("No normal tissue comparison. Merging a separate study would confound condition with batch.")
    )
  )
)

# ---- Server ---------------------------------------------------------------
server <- function(input, output, session) {

  updateSelectizeInput(session, "selectedGene", choices = gene_choices,
                       server = TRUE, selected = unname(gene_choices[1]))
  updateSelectizeInput(session, "highlight", choices = gene_choices, server = TRUE)

  # The "current gene" (a gene_id) is shared by click hooks and the dropdown
  current_gene <- reactiveVal(unname(gene_choices[1]))
  observeEvent(input$selectedGene, {
    if (nzchar(input$selectedGene)) current_gene(input$selectedGene)
  })

  # ---- classify genes ----
  filtered_res <- reactive({
    res_df %>%
      dplyr::mutate(
        Significance = dplyr::case_when(
          !is.na(pvalue) & pvalue < input$pval_thresh & log2FoldChange >=  input$lfc_thresh ~ "Up-regulated",
          !is.na(pvalue) & pvalue < input$pval_thresh & log2FoldChange <= -input$lfc_thresh ~ "Down-regulated",
          TRUE ~ "Not Significant"
        ),
        Significance = factor(Significance, levels = names(sig_colors))
      )
  })

  # Shared plot helpers
  label_annotations <- function(df, xcol, ycol) {
    lab <- df %>%
      dplyr::filter(Significance != "Not Significant") %>%
      dplyr::arrange(pvalue) %>%
      head(input$top_labels)
    lapply(seq_len(nrow(lab)), function(i) {
      list(x = lab[[xcol]][i], y = lab[[ycol]][i], text = lab$gene_name[i],
           showarrow = TRUE, arrowhead = 0, arrowsize = 0.5, ax = 20, ay = -25,
           font = list(size = 10))
    })
  }
  add_highlight <- function(p, df, xcol, ycol) {
    hl <- df[df$gene_id %in% input$highlight, , drop = FALSE]
    if (nrow(hl) == 0) return(p)
    add_trace(p, data = hl, x = hl[[xcol]], y = hl[[ycol]], type = "scatter", mode = "markers",
              inherit = FALSE, name = "Highlighted", customdata = hl$gene_id,
              marker = list(size = 13, color = "rgba(0,0,0,0)", line = list(color = "black", width = 2)),
              text = paste0("Gene: ", hl$gene_name), hoverinfo = "text")
  }

  # PCA (computed from the VST matrix, top 500 variable genes as in DESeq2::plotPCA)
  output$pcaPlot <- renderPlotly({
    v   <- apply(vst_mat, 1, var)
    top <- order(v, decreasing = TRUE)[seq_len(min(500, sum(v > 0)))]
    pc  <- prcomp(t(vst_mat[top, , drop = FALSE]), center = TRUE)
    pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2))
    df  <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], name = rownames(pc$x),
                      Group = colData[rownames(pc$x), grp_var])
    p <- ggplot(df, aes(PC1, PC2, color = Group, text = paste("Sample:", name))) +
      geom_point(size = 5, alpha = 0.9) +
      scale_color_manual(values = cond_colors) +
      labs(x = paste0("PC1: ", pct[1], "% variance"),
           y = paste0("PC2: ", pct[2], "% variance")) +
      theme_minimal()
    ggplotly(p, tooltip = c("text", "color"))
  })

  output$sampleDistPlot <- renderPlot({
    d <- dist(t(vst_mat))
    pheatmap(as.matrix(d), clustering_distance_rows = d, clustering_distance_cols = d,
             col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
  })

  # ---- Run QC ----
  output$qcNote <- renderUI({
    if (is.null(qc_df)) return(div(class = "alert alert-secondary", "No run QC table found in this data object."))
    if (any(qc_df$assigned_pct < 50, na.rm = TRUE)) {
      div(class = "alert alert-warning",
          "Low featureCounts assignment rate in at least one sample. Check strandedness (-s), ",
          "that the GTF matches the genome build, and paired-end settings.")
    }
  })
  output$qcTable <- renderDT({
    req(qc_df)
    datatable(qc_df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  output$qcPlot <- renderPlotly({
    req(qc_df)
    plot_ly(qc_df, x = ~sample) %>%
      add_bars(y = ~alignment_rate, name = "HISAT2 alignment %", marker = list(color = "#0072B2")) %>%
      add_bars(y = ~assigned_pct, name = "featureCounts assigned %", marker = list(color = "#E69F00")) %>%
      layout(barmode = "group", yaxis = list(title = "Percent", range = c(0, 100)),
             xaxis = list(title = ""), legend = list(orientation = "h"))
  })

  # ---- Volcano ----
  output$volcanoPlot <- renderPlotly({
    df <- filtered_res() %>%
      dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
      dplyr::mutate(neglog10p = -log10(pmax(pvalue, 1e-300)))

    plot_ly(
      data = df, x = ~log2FoldChange, y = ~neglog10p,
      color = ~Significance, colors = sig_colors,
      customdata = ~gene_id, source = "volcano",
      type = "scattergl", mode = "markers",
      marker = list(size = 6, opacity = 0.65),
      text = ~paste0("Gene: ", gene_name,
                     "<br>log2FC: ", round(log2FoldChange, 3),
                     "<br>p-value: ", signif(pvalue, 3)),
      hoverinfo = "text"
    ) %>%
      add_highlight(df, "log2FoldChange", "neglog10p") %>%
      layout(
        xaxis = list(title = paste0("Log2 Fold Change (", cmp_label, ")")),
        yaxis = list(title = "-Log10 p-value"),
        annotations = label_annotations(df, "log2FoldChange", "neglog10p"),
        shapes = list(
          list(type = "line", x0 = -input$lfc_thresh, x1 = -input$lfc_thresh,
               y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dash", color = "black", width = 1)),
          list(type = "line", x0 =  input$lfc_thresh, x1 =  input$lfc_thresh,
               y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dash", color = "black", width = 1)),
          list(type = "line", x0 = 0, x1 = 1, xref = "paper",
               y0 = -log10(input$pval_thresh), y1 = -log10(input$pval_thresh),
               line = list(dash = "dash", color = "black", width = 1))
        )
      ) %>%
      event_register("plotly_click")
  })

  # ---- MA plot ----
  output$maPlot <- renderPlotly({
    df <- filtered_res() %>%
      dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange), !is.na(baseMean), baseMean > 0) %>%
      dplyr::mutate(log10mean = log10(baseMean))

    plot_ly(
      data = df, x = ~log10mean, y = ~log2FoldChange,
      color = ~Significance, colors = sig_colors,
      customdata = ~gene_id, source = "ma",
      type = "scattergl", mode = "markers",
      marker = list(size = 6, opacity = 0.65),
      text = ~paste0("Gene: ", gene_name,
                     "<br>mean count: ", round(baseMean, 1),
                     "<br>log2FC: ", round(log2FoldChange, 3)),
      hoverinfo = "text"
    ) %>%
      add_highlight(df, "log10mean", "log2FoldChange") %>%
      layout(
        xaxis = list(title = "Log10 mean normalized count"),
        yaxis = list(title = paste0("Log2 Fold Change (", cmp_label, ")")),
        annotations = label_annotations(df, "log10mean", "log2FoldChange"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0,
               line = list(color = "black", width = 1)))
      ) %>%
      event_register("plotly_click")
  })

  # ---- Click-event hooks (volcano and MA both set the current gene) ----
  handle_click <- function(src) {
    observeEvent(event_data("plotly_click", source = src), {
      click <- event_data("plotly_click", source = src)
      gene  <- as.character(click$customdata[[1]])
      if (length(gene) == 1 && gene %in% unname(gene_choices)) {
        current_gene(gene)
        updateSelectizeInput(session, "selectedGene", choices = gene_choices,
                             selected = gene, server = TRUE)
      }
    })
  }
  handle_click("volcano")
  handle_click("ma")

  output$volcanoBox   <- renderPlot({ req(current_gene()); make_gene_plot(current_gene()) })
  output$maBox        <- renderPlot({ req(current_gene()); make_gene_plot(current_gene()) })
  output$geneExprPlot <- renderPlot({ req(current_gene()); make_gene_plot(current_gene()) })

  # ---- Multiple testing context ----
  output$fpNote <- renderUI({
    df <- filtered_res()
    n_tested <- sum(!is.na(df$pvalue))
    n_called <- sum(df$Significance != "Not Significant")
    expected <- round(n_tested * input$pval_thresh)
    div(class = "small text-muted",
        paste0(format(n_called, big.mark = ","), " genes pass both cutoffs. At p < ", input$pval_thresh,
               ", about ", format(expected, big.mark = ","), " of ", format(n_tested, big.mark = ","),
               " tested genes would pass by chance alone even if nothing were truly different."))
  })

  output$pvalHist <- renderPlot({
    d <- res_df[!is.na(res_df$pvalue), , drop = FALSE]
    ggplot(d, aes(pvalue)) +
      geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#0072B2", color = "white") +
      geom_vline(xintercept = input$pval_thresh, linetype = "dashed") +
      labs(x = "Raw p-value", y = "Number of genes") +
      theme_minimal(base_size = 14)
  })

  # ---- Heatmap of top differential genes ----
  heat_genes <- reactive({
    sig <- filtered_res() %>%
      dplyr::filter(Significance != "Not Significant") %>%
      dplyr::arrange(pvalue) %>%
      head(input$heat_n)
    fallback <- nrow(sig) < 2
    if (fallback) {
      sig <- filtered_res() %>% dplyr::filter(!is.na(pvalue)) %>%
        dplyr::arrange(pvalue) %>% head(input$heat_n)
    }
    list(genes = sig, fallback = fallback)
  })
  output$heatNote <- renderText({
    h <- heat_genes()
    if (h$fallback) paste0("Fewer than 2 genes pass the current cutoffs; showing the top ",
                           nrow(h$genes), " genes by p-value instead.")
    else paste0("Showing ", nrow(h$genes), " genes passing the current cutoffs, ranked by p-value.")
  })
  output$heatmapPlot <- renderPlot({
    g   <- heat_genes()$genes
    ids <- g$gene_id[g$gene_id %in% rownames(vst_mat)]
    m   <- vst_mat[ids, , drop = FALSE]
    m   <- m[apply(m, 1, sd) > 0, , drop = FALSE]
    req(nrow(m) >= 1)
    z <- t(scale(t(m)))
    rownames(z) <- gene_label(rownames(z))
    ann <- data.frame(Group = colData[[grp_var]], row.names = rownames(colData))
    pheatmap(z, annotation_col = ann, annotation_colors = list(Group = cond_colors),
             cluster_rows = nrow(z) >= 2, cluster_cols = TRUE,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             fontsize_row = 9, border_color = NA)
  })

  # ---- Gene inspector ----
  output$geneStatsTable <- renderTable({
    req(current_gene())
    res_df %>%
      dplyr::filter(gene_id == current_gene()) %>%
      dplyr::select(gene_name, gene_id, log2FoldChange, pvalue, padj)
  }, digits = 4)

  output$dl_gene_png <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_.-]", "_", gene_label(current_gene())), "_expression.png"),
    content  = function(file) ggsave(file, make_gene_plot(current_gene()), width = 6, height = 5, dpi = 150)
  )

  # ---- Pathways ----
  gsea_view <- reactive({
    validate(need(nrow(gsea_df) > 0,
                  "No pathway results in this data object (GSEA was skipped or found no sets)."))
    d <- gsea_df
    if (input$gsea_src != "All") d <- d[d$source == input$gsea_src, , drop = FALSE]
    validate(need(nrow(d) > 0, "No gene sets for this source."))
    d[order(d$pval), ]
  })
  output$gseaNote <- renderUI({
    if (min(n_per) < 2) div(class = "alert alert-warning",
      "Exploratory: with a single-sample group the ranking statistic is noisy. ",
      "Read pathways as hypotheses.")
  })
  output$gseaPlot <- renderPlot({
    d   <- gsea_view()
    n   <- input$gsea_top
    pos <- d[!is.na(d$NES) & d$NES > 0, , drop = FALSE]
    neg <- d[!is.na(d$NES) & d$NES < 0, , drop = FALSE]
    pos <- head(pos[order(-pos$NES), ], n)
    neg <- head(neg[order(neg$NES), ], n)
    top <- rbind(pos, neg)
    req(nrow(top) > 0)
    top$label <- gsub("^HALLMARK_", "", top$pathway)
    top$label <- factor(top$label, levels = top$label[order(top$NES)])
    top$Nominal <- ifelse(top$pval < 0.05, "p < 0.05", "p >= 0.05")
    ggplot(top, aes(NES, label, fill = Nominal)) +
      geom_col() +
      scale_fill_manual(values = c("p < 0.05" = "#D55E00", "p >= 0.05" = "grey70")) +
      labs(x = paste0("Normalized enrichment score (+ = higher in ", contrast_info$numerator, ")"),
           y = NULL, fill = NULL) +
      theme_minimal(base_size = 13)
  })
  output$gseaTable <- renderDT({
    d <- gsea_view()
    d <- d[, c("pathway", "source", "NES", "pval", "padj", "size", "leadingEdge")]
    datatable(d, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatRound(c("NES", "pval", "padj"), digits = 4)
  })
  output$dl_gsea <- downloadHandler(
    filename = function() "gsea_results.csv",
    content  = function(file) write.csv(gsea_view(), file, row.names = FALSE)
  )

  # ---- Results table ----
  table_data <- reactive({
    filtered_res() %>%
      dplyr::filter(Significance %in% input$sig_filter) %>%
      dplyr::select(gene_name, gene_id, log2FoldChange, pvalue, padj, Significance)
  })
  output$resultsTable <- renderDT({
    datatable(table_data(), options = list(pageLength = 15, autoWidth = TRUE),
              rownames = FALSE, filter = "top") %>%
      formatRound(columns = c("log2FoldChange", "pvalue", "padj"), digits = 4)
  })
  output$dl_table <- downloadHandler(
    filename = function() "de_results_filtered.csv",
    content  = function(file) write.csv(table_data(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)
