# Clustering module: Hierarchical + K-means over mvapp:: backends.
#  * Hierarchical: coloured dendrogram, trait heatmap, cluster distances, and an
#    ANOVA + Tukey letter-display comparison of clusters for any trait.
#  * K-means: cluster plot, optimal-k (elbow / silhouette / gap), and the same
#    ANOVA comparison.
# Reuses .corr_num_phenos() from mod_correlations.R.

# --- shared cluster-vs-trait comparison helpers -----------------------------

# map cluster labels back onto the original data rows and pull one trait
.cluster_compare_df <- function(ds, clusters, trait) {
  idcols <- c(ds$roles$id, ds$roles$group, ds$roles$time)
  key <- make.unique(do.call(paste, c(ds$data[idcols], sep = "_")))
  idx <- match(clusters$id, key)
  data.frame(cluster = factor(clusters$cluster),
             value = ds$data[[trait]][idx])
}

# screen every numeric trait for significant differences across the clusters
.cluster_screen <- function(ds, clusters) {
  idcols <- c(ds$roles$id, ds$roles$group, ds$roles$time)
  key <- make.unique(do.call(paste, c(ds$data[idcols], sep = "_")))
  idx <- match(clusters$id, key)
  df <- ds$data[idx, , drop = FALSE]
  df$.cluster <- factor(clusters$cluster)
  traits <- ds$roles$phenotypes[
    vapply(ds$data[ds$roles$phenotypes], is.numeric, logical(1))]
  sc <- anova_screen(df, traits, ".cluster")
  sc <- sc[!is.na(sc$anova_p) & sc$anova_p < 0.05, , drop = FALSE]
  sc$anova_p <- signif(sc$anova_p, 3)
  sc
}

# choices for the grouped-distribution plot type (shared)
.compare_geoms <- c("Box plot" = "box", "Box + points" = "box_jitter",
  "Violin" = "violin", "Violin + points" = "violin_jitter",
  "Bar +/- SE" = "bar_se", "Bar +/- SD" = "bar_sd",
  "Mean +/- SE + points" = "points_se")

# plot a trait across groups (clusters) with Tukey compact-letter display, in
# the chosen geom + palette. Reused for hierarchical and k-means comparisons.
.compare_plot <- function(cmp, ylab, res, geom = "box", palette = "ggsci: D3",
                          xlab = "cluster", test_label = "One-way ANOVA") {
  cmp$cluster <- factor(cmp$cluster)
  s <- res$summary; s$cluster <- factor(s$group)
  rng <- diff(range(cmp$value, na.rm = TRUE)); off <- 0.06 * (if (rng > 0) rng else 1)
  s$lab_y <- switch(geom,
    bar_se = s$mean + s$se, bar_sd = s$mean + s$sd, s$ymax) + off
  nlev <- nlevels(cmp$cluster)

  p <- switch(geom,
    box = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_boxplot(alpha = 0.85, colour = "grey30"),
    box_jitter = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_boxplot(alpha = 0.6, colour = "grey30", outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1),
    violin = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_violin(alpha = 0.85, colour = "grey30"),
    violin_jitter = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_violin(alpha = 0.6, colour = "grey30") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1),
    bar_se = ggplot2::ggplot(s, ggplot2::aes(cluster, mean, fill = cluster)) +
      ggplot2::geom_col(alpha = 0.85, colour = "grey30") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.2),
    bar_sd = ggplot2::ggplot(s, ggplot2::aes(cluster, mean, fill = cluster)) +
      ggplot2::geom_col(alpha = 0.85, colour = "grey30") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), width = 0.2),
    points_se = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, colour = cluster)) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
      ggplot2::stat_summary(fun = mean, geom = "point", colour = "black", size = 2,
        fun.min = function(z) mean(z) - stats::sd(z) / sqrt(length(z)),
        fun.max = function(z) mean(z) + stats::sd(z) / sqrt(length(z))) +
      ggplot2::stat_summary(fun.data = function(z) {
        m <- mean(z); e <- stats::sd(z) / sqrt(length(z))
        data.frame(y = m, ymin = m - e, ymax = m + e)
      }, geom = "errorbar", colour = "black", width = 0.15))

  sc <- cat_scale(palette, nlev, if (geom == "points_se") "colour" else "fill")
  p <- p + ggplot2::geom_text(data = s, inherit.aes = FALSE, fontface = "bold",
      ggplot2::aes(x = cluster, y = lab_y, label = letter)) +
    ggplot2::labs(x = xlab, y = ylab, subtitle = sprintf(
      "%s p = %.3g; groups sharing a letter are not significantly different.",
      test_label, res$anova_p)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "none")
  if (!is.null(sc)) p <- p + sc
  p
}

# Faceted version: a separate one-way ANOVA + Tukey letter display is run within
# each panel (level of `facet`). `cmp` has columns value, cluster, facet.
.compare_plot_faceted <- function(cmp, ylab, geom = "box", palette = "ggsci: D3",
                                  xlab = "group", facet_label = "facet",
                                  letter_fun = anova_tukey_letters,
                                  test_label = "one-way ANOVA + Tukey HSD") {
  cmp$cluster <- factor(cmp$cluster); cmp$facet <- factor(cmp$facet)
  parts <- lapply(split(cmp, cmp$facet), function(sub) {
    sub$cluster <- droplevels(factor(sub$cluster))
    r <- tryCatch(letter_fun(sub, "value", "cluster"),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    s <- r$summary; s$facet <- sub$facet[1]
    rng <- diff(range(sub$value, na.rm = TRUE)); off <- 0.06 * (if (rng > 0) rng else 1)
    s$lab_y <- switch(geom, bar_se = s$mean + s$se, bar_sd = s$mean + s$sd, s$ymax) + off
    s
  })
  s <- do.call(rbind, parts)
  if (is.null(s)) stop("No panel had two comparable groups.")
  s$cluster <- factor(s$group)
  nlev <- nlevels(cmp$cluster)

  base <- switch(geom,
    box = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_boxplot(alpha = 0.85, colour = "grey30"),
    box_jitter = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_boxplot(alpha = 0.6, colour = "grey30", outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 0.8),
    violin = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_violin(alpha = 0.85, colour = "grey30"),
    violin_jitter = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, fill = cluster)) +
      ggplot2::geom_violin(alpha = 0.6, colour = "grey30") +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 0.8),
    bar_se = ggplot2::ggplot(s, ggplot2::aes(cluster, mean, fill = cluster)) +
      ggplot2::geom_col(alpha = 0.85, colour = "grey30") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.2),
    bar_sd = ggplot2::ggplot(s, ggplot2::aes(cluster, mean, fill = cluster)) +
      ggplot2::geom_col(alpha = 0.85, colour = "grey30") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), width = 0.2),
    points_se = ggplot2::ggplot(cmp, ggplot2::aes(cluster, value, colour = cluster)) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
      ggplot2::stat_summary(fun = mean, geom = "point", colour = "black", size = 1.8) +
      ggplot2::stat_summary(fun.data = function(z) {
        m <- mean(z); e <- stats::sd(z) / sqrt(length(z))
        data.frame(y = m, ymin = m - e, ymax = m + e)
      }, geom = "errorbar", colour = "black", width = 0.15))

  sc <- cat_scale(palette, nlev, if (geom == "points_se") "colour" else "fill")
  p <- base + ggplot2::facet_wrap(~facet, scales = "free_y") +
    ggplot2::geom_text(data = s, inherit.aes = FALSE, fontface = "bold", size = 3,
      ggplot2::aes(x = cluster, y = lab_y, label = letter)) +
    ggplot2::labs(x = xlab, y = ylab, subtitle = sprintf(paste0(
      "A separate %s is run within each %s panel; letters are comparable only ",
      "within a panel."), test_label, facet_label)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "none")
  if (!is.null(sc)) p <- p + sc
  p
}

# a reusable "compare clusters" tabPanel body (trait picker + plot + code/legend)
.compare_ui <- function(ns, prefix) {
  tabPanel("Compare clusters (ANOVA)",
    checkboxInput(ns(paste0(prefix, "_screen_chk")),
      "List traits where clusters differ significantly", FALSE),
    conditionalPanel(sprintf("input.%s_screen_chk == true", prefix), ns = ns,
      h5("Traits with significant cluster differences (one-way ANOVA, p < 0.05)"),
      DT::DTOutput(ns(paste0(prefix, "_screen"))),
      dl_btn(ns(paste0(prefix, "_screen_dl"))), tags$hr()),
    uiOutput(ns(paste0(prefix, "_trait_ui"))),
    fluidRow(
      column(6, selectInput(ns(paste0(prefix, "_geom")), "Plot type",
        .compare_geoms)),
      column(6, selectInput(ns(paste0(prefix, "_pal")), "Colour palette",
        cat_palette_choices, selected = "ggsci: D3"))),
    checkboxInput(ns(paste0(prefix, "_show_legend")), "Show figure legend"),
    conditionalPanel(sprintf("input.%s_show_legend == true", prefix), ns = ns,
      textOutput(ns(paste0(prefix, "_legend")))),
    plotOutput(ns(paste0(prefix, "_plot")), height = "440px"),
    checkboxInput(ns(paste0(prefix, "_show_code")), "Show me the R code"),
    conditionalPanel(sprintf("input.%s_show_code == true", prefix), ns = ns,
      verbatimTextOutput(ns(paste0(prefix, "_code")))),
    tags$hr(),
    DT::DTOutput(ns(paste0(prefix, "_tab"))),
    dl_btn(ns(paste0(prefix, "_tab_dl"))))
}

mod_clustering_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    # ======================================================= Hierarchical =====
    tabPanel("Hierarchical",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("hc_pheno_ui")),
          checkboxInput(ns("hc_scale"), "Scale traits to unit variance", TRUE),
          selectInput(ns("hc_dist"), "Distance metric",
                      c("euclidean", "manhattan", "maximum", "canberra")),
          selectInput(ns("hc_link"), "Linkage",
                      c("complete", "average", "ward.D2", "single", "mcquitty")),
          numericInput(ns("hc_k"), "Number of clusters", 3, 1, 20),
          uiOutput(ns("hc_annot_ui"))
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Dendrogram", plotOutput(ns("dendro"), height = "520px")),
            tabPanel("Heatmap",
              fluidRow(
                column(6, selectInput(ns("hc_hmpal"), "Colour scale",
                  .corr_palette_choices, selected = "RdBu")),
                column(6, checkboxInput(ns("hc_hm_legend"), "Show figure legend"))),
              conditionalPanel("input.hc_hm_legend == true", ns = ns,
                textOutput(ns("hc_hm_legend_txt"))),
              plotOutput(ns("heatmap"), height = "560px")),
            tabPanel("Cluster distances",
              p("How far apart the clusters are, in the chosen distance metric on ",
                "the (scaled) traits. ",
                tags$b("Off-diagonal"), " cells are the mean distance between the ",
                "members of two different clusters - larger values mean the ",
                "clusters are more distinct and better separated. ",
                tags$b("Diagonal"), " cells are the mean distance among members of ",
                "the same cluster - smaller values mean a tighter, more ",
                "homogeneous cluster. Comparing the two tells you whether clusters ",
                "are genuinely separated (between-distances clearly larger than ",
                "within-distances)."),
              DT::DTOutput(ns("hc_dist_tab")), dl_btn(ns("hc_dist_dl"))),
            .compare_ui(ns, "hc_cmp")
          ),
          tags$hr(),
          DT::DTOutput(ns("hc_table")), dl_btn(ns("hc_table_dl")),
          checkboxInput(ns("hc_show_code"), "Show me the R code (clustering)"),
          conditionalPanel("input.hc_show_code == true", ns = ns,
            verbatimTextOutput(ns("hc_code")))
        )
      )
    ),
    # ============================================================ K-means =====
    tabPanel("K-means",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("km_pheno_ui")),
          checkboxInput(ns("km_scale"), "Scale traits to unit variance", TRUE),
          numericInput(ns("km_k"), "Number of clusters (k)", 3, 1, 20),
          numericInput(ns("km_maxk"), "Max k to evaluate", 10, 3, 20),
          actionButton(ns("km_go"), "Run k-means", class = "btn-primary"),
          helpText("The 'Optimal k' tab evaluates elbow, silhouette and gap",
                   "statistic to help you choose k.")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Clusters",
              checkboxInput(ns("km_show_legend"), "Show figure legend"),
              conditionalPanel("input.km_show_legend == true", ns = ns,
                textOutput(ns("km_plot_legend"))),
              plotly::plotlyOutput(ns("km_plot"), height = "440px")),
            tabPanel("Optimal k",
              radioButtons(ns("km_method"), "Method",
                c("Elbow (within-SS)" = "wss", "Silhouette" = "silhouette",
                  "Gap statistic" = "gap"), inline = TRUE),
              checkboxInput(ns("km_opt_show_legend"), "Show figure legend"),
              conditionalPanel("input.km_opt_show_legend == true", ns = ns,
                textOutput(ns("km_opt_legend"))),
              plotOutput(ns("km_opt"), height = "420px"),
              textOutput(ns("km_opt_txt"))),
            .compare_ui(ns, "km_cmp")
          ),
          tags$hr(),
          DT::DTOutput(ns("km_table")), dl_btn(ns("km_table_dl")),
          checkboxInput(ns("km_show_code"), "Show me the R code (k-means)"),
          conditionalPanel("input.km_show_code == true", ns = ns,
            verbatimTextOutput(ns("km_code")))
        )
      )
    )
  )
}

mod_clustering_server <- function(id, dataset, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    logmsg <- function(m) if (is.function(log_add)) log_add(m)

    # =================================================== Hierarchical ========
    output$hc_pheno_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("hc_phenos"), "Phenotypes", num, selected = num, multiple = TRUE)
    })

    hc <- reactive({
      ds <- dataset(); req(ds, input$hc_phenos)
      validate(need(length(input$hc_phenos) >= 2, "Select at least two phenotypes."))
      run_hclust(ds, phenotypes = input$hc_phenos, scale = input$hc_scale,
                 dist_method = input$hc_dist, link_method = input$hc_link,
                 k = input$hc_k)
    })

    output$dendro <- renderPlot({
      res <- hc()
      dend <- stats::as.dendrogram(res$hclust)
      dend <- dendextend::color_branches(dend, k = res$k)
      dend <- dendextend::color_labels(dend, k = res$k)
      graphics::par(mar = c(3, 1, 2, 7))
      plot(dend, horiz = TRUE,
           main = sprintf("Hierarchical clustering (%s linkage, %s distance) - %d clusters",
                          res$link_method, res$dist_method, res$k))
    })

    output$hc_annot_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("hc_annot"), "Annotate heatmap with (side colours)",
                  .corr_factors(ds), multiple = TRUE)
    })

    output$heatmap <- renderPlot({
      res <- hc(); ds <- dataset()
      ann <- data.frame(cluster = factor(res$clusters$cluster))
      rownames(ann) <- res$clusters$id
      if (length(input$hc_annot)) {                      # extra categorical stripes
        idcols <- c(ds$roles$id, ds$roles$group, ds$roles$time)
        key <- make.unique(do.call(paste, c(ds$data[idcols], sep = "_")))
        idx <- match(res$clusters$id, key)
        for (f in input$hc_annot) ann[[f]] <- factor(ds$data[[f]][idx])
      }
      pal <- if (is.null(input$hc_hmpal)) "RdBu" else input$hc_hmpal
      pheatmap::pheatmap(res$x, cluster_rows = res$hclust,
        cutree_rows = res$k,                              # gaps between clusters
        annotation_row = ann, scale = if (res$scaled) "none" else "column",
        color = .corr_palette(pal, 11),
        show_rownames = nrow(res$x) <= 50, fontsize_row = 7,
        main = "Scaled traits (columns) across samples (rows)")
    })

    output$hc_hm_legend_txt <- renderText({
      res <- hc()
      sprintf(paste0("Heatmap of %d scaled traits (columns) across %d samples ",
        "(rows). Samples are ordered by the hierarchical clustering tree and ",
        "split into %d clusters (white gaps); the coloured stripe(s) on the left ",
        "mark cluster membership%s. Cell colour shows the trait value (%s): low ",
        "values appear at one end of the colour key and high values at the other, ",
        "with the mid-tone marking values near the trait's average across samples. ",
        "Colour scale: %s."),
        ncol(res$x), nrow(res$x), res$k,
        if (length(input$hc_annot)) paste0(" and ",
          paste(input$hc_annot, collapse = ", ")) else "",
        if (res$scaled) "z-score across samples" else "column-scaled",
        if (is.null(input$hc_hmpal)) "RdBu" else input$hc_hmpal)
    })

    hc_dist_df <- reactive({
      M <- cluster_distances(hc())
      data.frame(cluster = rownames(M), round(M, 3), check.names = FALSE)
    })
    output$hc_dist_tab <- DT::renderDT(
      DT::datatable(hc_dist_df(), rownames = FALSE, options = list(dom = "t")))
    output$hc_dist_dl <- csv_dl(hc_dist_df, "cluster_distances.csv")

    output$hc_table <- DT::renderDT(
      DT::datatable(hc()$clusters, rownames = FALSE, options = list(scrollX = TRUE)))
    output$hc_table_dl <- csv_dl(reactive(hc()$clusters), "hierarchical_clusters.csv")
    output$hc_code <- renderText(
      code_hclust(input$hc_phenos, scale = input$hc_scale,
                  dist_method = input$hc_dist, link_method = input$hc_link,
                  k = input$hc_k))
    observeEvent(list(input$hc_k, input$hc_phenos, input$hc_link, input$hc_dist), {
      req(input$hc_phenos)
      logmsg(sprintf("Hierarchical clustering: %s linkage, %s distance, cut into %d clusters, on %d traits (%s).",
                     input$hc_link, input$hc_dist, input$hc_k,
                     length(input$hc_phenos),
                     paste(input$hc_phenos, collapse = ", ")))
    }, ignoreInit = TRUE)

    # cluster comparison (hierarchical)
    output$hc_cmp_trait_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("hc_cmp_trait"), "Trait to compare across clusters",
                  .corr_num_phenos(ds))
    })
    hc_cmp <- reactive({
      ds <- dataset(); req(ds, input$hc_cmp_trait)
      cmp <- .cluster_compare_df(ds, hc()$clusters, input$hc_cmp_trait)
      list(cmp = cmp, res = anova_tukey_letters(cmp, "value", "cluster"))
    })
    output$hc_cmp_plot <- renderPlot({
      x <- hc_cmp()
      .compare_plot(x$cmp, input$hc_cmp_trait, x$res, input$hc_cmp_geom,
                    input$hc_cmp_pal)
    })
    output$hc_cmp_tab <- DT::renderDT({
      s <- hc_cmp()$res$summary
      s[c("mean", "sd", "se")] <- lapply(s[c("mean", "sd", "se")], round, 3)
      DT::datatable(s[c("group", "n", "mean", "sd", "se", "letter")],
                    rownames = FALSE, options = list(dom = "t"))
    })
    output$hc_cmp_tab_dl <- csv_dl(reactive(hc_cmp()$res$summary),
                                   "hierarchical_cluster_anova.csv")
    output$hc_cmp_legend <- renderText({
      x <- hc_cmp()
      legend_anova_tukey(input$hc_cmp_trait, "cluster", anova_p = x$res$anova_p,
                         n_groups = nrow(x$res$summary), geom = input$hc_cmp_geom)
    })
    output$hc_cmp_code <- renderText(
      code_anova_tukey(input$hc_cmp_trait, "cluster"))
    hc_cmp_screen_df <- reactive(.cluster_screen(dataset(), hc()$clusters))
    output$hc_cmp_screen <- DT::renderDT(DT::datatable(hc_cmp_screen_df(),
      rownames = FALSE, options = list(pageLength = 10)))
    output$hc_cmp_screen_dl <- csv_dl(hc_cmp_screen_df,
      "hierarchical_cluster_significant_traits.csv")

    # ======================================================= K-means =========
    output$km_pheno_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("km_phenos"), "Phenotypes", num, selected = num, multiple = TRUE)
    })

    km <- eventReactive(input$km_go, {
      ds <- dataset(); req(ds, input$km_phenos)
      validate(need(length(input$km_phenos) >= 2, "Select at least two phenotypes."))
      res <- run_kmeans(ds, phenotypes = input$km_phenos, scale = input$km_scale,
                        k = input$km_k)
      logmsg(sprintf("K-means: k = %d (between/total SS = %.3f), on %d traits (%s).",
                     input$km_k, res$kmeans$betweenss / res$kmeans$totss,
                     length(input$km_phenos),
                     paste(input$km_phenos, collapse = ", ")))
      res
    })

    # fast: elbow + silhouette (computed with every run)
    kdiag <- eventReactive(input$km_go, {
      ds <- dataset(); req(input$km_phenos)
      validate(need(length(input$km_phenos) >= 2, "Select at least two phenotypes."))
      cluster_diagnostics(ds, phenotypes = input$km_phenos, scale = input$km_scale,
                          max_k = input$km_maxk, gap = FALSE)
    })
    # slow: gap statistic, only computed when the user views the Gap method
    kgap <- eventReactive(input$km_go, {
      ds <- dataset(); req(input$km_phenos)
      withProgress(message = "Computing gap statistic (bootstrap)...", value = 0.5,
        cluster_diagnostics(ds, phenotypes = input$km_phenos, scale = input$km_scale,
                            max_k = input$km_maxk, gap = TRUE, B = 25))
    })

    km_plotdf <- reactive({
      res <- km(); ds <- dataset()
      pca <- run_pca(ds, phenotypes = input$km_phenos, scale = input$km_scale)
      d <- merge(pca$ind[, c("id", "Dim.1", "Dim.2")], res$clusters, by = "id")
      d$cluster <- factor(d$cluster)
      list(df = d, var1 = pca$eig$variance_percent[1],
           var2 = pca$eig$variance_percent[2])
    })

    output$km_plot <- plotly::renderPlotly({
      validate(need(input$km_go, "Click 'Run k-means' first to compute clusters."))
      x <- km_plotdf(); d <- x$df
      p <- ggplot2::ggplot(d, ggplot2::aes(x = Dim.1, y = Dim.2,
                                           colour = cluster, text = id)) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::labs(title = "K-means clusters (projected onto PC1 vs PC2)",
                      x = sprintf("PC1 (%.1f%%)", x$var1),
                      y = sprintf("PC2 (%.1f%%)", x$var2)) +
        ggplot2::theme_minimal(base_size = 13)
      if (nlevels(d$cluster) <= 10) p <- p + ggsci::scale_colour_d3()
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "colour"))
    })

    output$km_plot_legend <- renderText({
      req(input$km_go); x <- km_plotdf()
      sprintf(paste0("K-means clusters shown in the space of the first two ",
        "principal components of the same %d scaled traits used for clustering. ",
        "PC1 and PC2 capture %.1f%% and %.1f%% of the total variation. Note: ",
        "k-means is computed in the full trait space (all %d traits), not on the ",
        "PCs - the PCA is used only to project the clusters into 2-D so they can ",
        "be seen; well-separated clusters may still overlap here if their ",
        "separation lies along higher components."),
        length(input$km_phenos), x$var1, x$var2, length(input$km_phenos))
    })

    output$km_opt <- renderPlot({
      validate(need(input$km_go, "Click 'Run k-means' first to evaluate the optimal number of clusters."))
      m <- input$km_method
      dg <- if (m == "gap") kgap() else kdiag()
      tab <- dg$table
      ylab <- c(wss = "Total within-cluster SS", silhouette = "Average silhouette width",
                gap = "Gap statistic")[[m]]
      d <- data.frame(k = tab$k, y = tab[[m]])
      d <- d[!is.na(d$y), ]
      p <- ggplot2::ggplot(d, ggplot2::aes(k, y)) +
        ggplot2::geom_line(colour = "#4C78A8") +
        ggplot2::geom_point(size = 2, colour = "#4C78A8") +
        ggplot2::scale_x_continuous(breaks = tab$k) +
        ggplot2::labs(x = "Number of clusters (k)", y = ylab,
                      title = sprintf("Optimal k - %s method", m)) +
        ggplot2::theme_minimal(base_size = 13)
      sug <- if (m == "silhouette") dg$best$silhouette else
             if (m == "gap") dg$best$gap else input$km_k
      if (!is.null(sug) && !is.na(sug))
        p <- p + ggplot2::geom_vline(xintercept = sug, linetype = "dashed",
                                     colour = "#E45756")
      p
    })
    output$km_opt_legend <- renderText({
      req(input$km_go)
      best <- if (input$km_method == "silhouette") kdiag()$best$silhouette
              else if (input$km_method == "gap") kgap()$best$gap else NA
      legend_optimal_k(input$km_method, best = best, B = 25)
    })
    output$km_opt_txt <- renderText({
      req(input$km_go)
      b <- kdiag()$best
      gaptxt <- if (input$km_method == "gap")
        sprintf("; gap statistic (firstSEmax): %s", kgap()$best$gap)
      else "; select the 'Gap statistic' method to compute it (slower, bootstrapped)"
      sprintf("Suggested k - silhouette: %s%s. (For the elbow, look for the 'knee' where the curve flattens.)",
              b$silhouette, gaptxt)
    })

    output$km_table <- DT::renderDT(
      DT::datatable(km()$clusters, rownames = FALSE, options = list(scrollX = TRUE)))
    output$km_table_dl <- csv_dl(reactive(km()$clusters), "kmeans_clusters.csv")
    output$km_code <- renderText(
      code_kmeans(input$km_phenos, scale = input$km_scale, k = input$km_k))

    # cluster comparison (k-means)
    output$km_cmp_trait_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("km_cmp_trait"), "Trait to compare across clusters",
                  .corr_num_phenos(ds))
    })
    km_cmp <- reactive({
      ds <- dataset(); req(ds, input$km_cmp_trait)
      cmp <- .cluster_compare_df(ds, km()$clusters, input$km_cmp_trait)
      list(cmp = cmp, res = anova_tukey_letters(cmp, "value", "cluster"))
    })
    output$km_cmp_plot <- renderPlot({
      validate(need(input$km_go, "Click 'Run k-means' first to compute clusters."))
      x <- km_cmp()
      .compare_plot(x$cmp, input$km_cmp_trait, x$res, input$km_cmp_geom,
                    input$km_cmp_pal)
    })
    output$km_cmp_tab <- DT::renderDT({
      s <- km_cmp()$res$summary
      s[c("mean", "sd", "se")] <- lapply(s[c("mean", "sd", "se")], round, 3)
      DT::datatable(s[c("group", "n", "mean", "sd", "se", "letter")],
                    rownames = FALSE, options = list(dom = "t"))
    })
    output$km_cmp_tab_dl <- csv_dl(reactive(km_cmp()$res$summary),
                                   "kmeans_cluster_anova.csv")
    output$km_cmp_legend <- renderText({
      x <- km_cmp()
      legend_anova_tukey(input$km_cmp_trait, "cluster", anova_p = x$res$anova_p,
                         n_groups = nrow(x$res$summary), geom = input$km_cmp_geom)
    })
    output$km_cmp_code <- renderText(
      code_anova_tukey(input$km_cmp_trait, "cluster"))
    km_cmp_screen_df <- reactive({
      req(input$km_go); .cluster_screen(dataset(), km()$clusters)
    })
    output$km_cmp_screen <- DT::renderDT(DT::datatable(km_cmp_screen_df(),
      rownames = FALSE, options = list(pageLength = 10)))
    output$km_cmp_screen_dl <- csv_dl(km_cmp_screen_df,
      "kmeans_cluster_significant_traits.csv")
  })
}
