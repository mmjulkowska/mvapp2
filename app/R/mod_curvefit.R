# Curve-fitting module: fit a model per sample (individual plant / sample ID),
# curate by R^2, view fit-plots (single or a grid of all samples), and compare
# the extracted parameters (DELTA / INTERCEPT / R2) across a factor, optionally
# faceted. Reuses .corr_num_phenos / .corr_factors (mod_correlations.R),
# .compare_plot / .compare_plot_faceted / .compare_geoms (mod_clustering.R),
# cat_palette_choices / legend_ui / csv_dl / dl_btn (util.R).

.cf_models <- c("Linear" = "linear", "Quadratic" = "quadratic",
  "Exponential" = "exponential", "Square-root" = "square_root",
  "Logarithmic" = "logarithmic")
.cf_numcols <- function(ds) names(ds$data)[vapply(ds$data, is.numeric, logical(1))]

# per-facet ANOVA + Tukey summary table (one block per facet level)
.cf_faceted_summary <- function(cmp) {
  parts <- lapply(split(cmp, cmp$facet), function(sub) {
    sub$cluster <- droplevels(factor(sub$cluster))
    r <- tryCatch(anova_tukey_letters(sub, "value", "cluster"),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    s <- r$summary
    data.frame(facet = as.character(sub$facet[1]), s, stringsAsFactors = FALSE)
  })
  do.call(rbind, parts)
}

mod_curvefit_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    # ---------------------------------------------------------------- Fit
    tabPanel("Fit",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("resp_ui")),
          uiOutput(ns("time_ui")),
          uiOutput(ns("id_ui")),
          uiOutput(ns("group_ui")),
          selectInput(ns("model"), "Model", .cf_models),
          numericInput(ns("r2cut"), "R² cut-off (flag fits below)", 0.7, 0, 1, 0.05),
          actionButton(ns("go"), "Fit curves", class = "btn-primary"),
          checkboxInput(ns("estimate"),
            "Also estimate the best model per group (all plants pooled)", FALSE)),
        mainPanel(
          legend_ui(ns, "cf", "Show explanation of this table"),
          textOutput(ns("cf_msg")),
          checkboxInput(ns("cf_code_chk"), "Show me the R code"),
          conditionalPanel("input.cf_code_chk == true", ns = ns,
            verbatimTextOutput(ns("cf_code"))),
          tags$hr(),
          h5("Per-plant fit (worst R² first)"),
          DT::DTOutput(ns("fit_tab")),
          fluidRow(
            column(6, dl_btn(ns("fit_dl"))),
            column(6, downloadButton(ns("curated_dl"),
              "Download curated data (R² ≥ cut-off)", class = "btn-sm"))),
          conditionalPanel("input.estimate == true", ns = ns,
            tags$hr(),
            h5("Best model per group, pooling all plants (by adjusted R²)"),
            DT::DTOutput(ns("estimate_tab")), dl_btn(ns("estimate_dl")))))),
    # ------------------------------------------------------------ Fit-plot
    tabPanel("Fit-plot",
      sidebarLayout(
        sidebarPanel(
          radioButtons(ns("fp_view"), "View",
            c("Grid of samples" = "grid", "Single sample" = "single"),
            selected = "grid"),
          conditionalPanel("input.fp_view == 'single'", ns = ns,
            uiOutput(ns("sample_ui"))),
          conditionalPanel("input.fp_view == 'grid'", ns = ns,
            uiOutput(ns("fp_sort_ui")),
            radioButtons(ns("fp_sortdir"), "Order",
              c("Ascending" = "asc", "Descending" = "desc"), selected = "asc"),
            fluidRow(
              column(6, numericInput(ns("fp_ncol"), "Columns", 3, 1, 8, 1)),
              column(6, numericInput(ns("fp_nrow"), "Rows per page", 3, 1, 8, 1))),
            uiOutput(ns("fp_page_ui")),
            helpText("Each page shows Columns × Rows panels. Sort the panels,",
                     "then step through pages; the plot area scrolls if it runs",
                     "past the window."))),
        mainPanel(
          legend_ui(ns, "fp"),
          div(style = "overflow-y: auto; max-height: 78vh;",
            plotOutput(ns("fit_plot"), height = "auto"))))),
    # -------------------------------------------------- Examine differences
    tabPanel("Examine differences",
      sidebarLayout(
        sidebarPanel(
          selectInput(ns("param"), "Parameter to compare",
            c("DELTA (rate)" = "DELTA", "INTERCEPT" = "INTERCEPT", "R²" = "r2")),
          uiOutput(ns("byfac_ui")),
          uiOutput(ns("exfacet_ui")),
          selectInput(ns("ex_geom"), "Plot type", .compare_geoms),
          selectInput(ns("ex_pal"), "Colour palette", cat_palette_choices,
            selected = "ggsci: D3")),
        mainPanel(
          legend_ui(ns, "ex"),
          plotOutput(ns("ex_plot"), height = "460px"),
          h5("Group summary (Tukey letters)"),
          DT::DTOutput(ns("ex_tab")), dl_btn(ns("ex_dl")))))
  )
}

mod_curvefit_server <- function(id, dataset, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    logmsg <- function(m) if (is.function(log_add)) log_add(m)

    output$resp_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("resp"), "Response (phenotype to model)", .corr_num_phenos(ds))
    })
    output$time_ui <- renderUI({
      ds <- dataset(); req(ds); nc <- .cf_numcols(ds)
      sel <- if (length(ds$roles$time) && ds$roles$time %in% nc) ds$roles$time else nc[1]
      selectInput(ns("time"), "Time / dose (x-axis)", nc, selected = sel)
    })
    # one curve per unique combination of these columns -- the individual plant
    output$id_ui <- renderUI({
      ds <- dataset(); req(ds); f <- .corr_factors(ds)
      sel <- intersect(ds$roles$id, f); if (!length(sel)) sel <- f[1]
      selectInput(ns("samples"), "Fit one curve per (sample / plant ID)",
        f, selected = sel, multiple = TRUE)
    })
    # factors carried along for comparing parameters & estimating best model;
    # the column acting as time / dose is excluded from these choices.
    output$group_ui <- renderUI({
      ds <- dataset(); req(ds)
      avail <- setdiff(.corr_factors(ds), input$time)
      guess <- setdiff(intersect(ds$roles$group, avail), isolate(input$samples))
      cur <- intersect(isolate(input$groups), avail)
      sel <- if (length(cur)) cur else guess
      selectInput(ns("groups"), "Grouping factors (genotype, treatment, …)",
        avail, selected = sel, multiple = TRUE)
    })

    fit <- eventReactive(input$go, {
      ds <- dataset(); req(input$resp, input$time, input$samples)
      validate(need(!(input$time %in% c(input$resp, input$samples)),
        "Time must differ from the response and the sample ID."))
      res <- fit_curves(ds, input$time, input$resp, input$samples,
                        groups = input$groups, model = input$model)
      nb <- sum(res$table$r2 < input$r2cut, na.rm = TRUE)
      logmsg(sprintf(paste0("Curve fit: %s model of %s ~ %s, one curve per %s",
        "%s; %d/%d plants below R2 %.2f."),
        input$model, input$resp, input$time, paste(input$samples, collapse = "+"),
        if (length(input$groups)) sprintf(" (grouped by %s)",
          paste(input$groups, collapse = "+")) else "",
        nb, nrow(res$table), input$r2cut))
      res
    })

    output$cf_msg <- renderText({
      res <- fit(); nb <- sum(res$table$r2 < input$r2cut, na.rm = TRUE)
      sprintf("%d of %d plants have R² below the cut-off (%.2f) with the %s model (median R² = %.2f).",
        nb, nrow(res$table), input$r2cut, input$model,
        stats::median(res$table$r2, na.rm = TRUE))
    })

    fit_tbl <- reactive({
      t <- fit()$table
      t[c("INTERCEPT", "DELTA", "r2")] <-
        lapply(t[c("INTERCEPT", "DELTA", "r2")], round, 4)
      t$sample <- NULL; t
    })
    output$fit_tab <- DT::renderDT(DT::datatable(fit_tbl(), rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE)))
    output$fit_dl <- csv_dl(fit_tbl, "curve_fit_parameters.csv")

    # curated raw data: keep only rows belonging to well-fit plants
    output$curated_dl <- downloadHandler(
      filename = function() "curated_data_r2_above_cutoff.csv",
      content = function(f) {
        res <- fit(); ds <- dataset()
        good <- res$table$sample[!is.na(res$table$r2) & res$table$r2 >= input$r2cut]
        key <- do.call(paste, c(ds$data[res$samples], sep = " | "))
        utils::write.csv(ds$data[key %in% good, , drop = FALSE], f, row.names = FALSE)
        logmsg(sprintf("Downloaded curated data: kept %d well-fit plants (R2 >= %.2f).",
          length(good), input$r2cut))
      })

    output$legend_cf <- renderText({
      res <- fit()
      legend_curvefit(input$resp, input$time, input$model,
        n_samples = nrow(res$table),
        median_r2 = stats::median(res$table$r2, na.rm = TRUE),
        r2_cut = input$r2cut,
        n_below = sum(res$table$r2 < input$r2cut, na.rm = TRUE))
    })
    output$cf_code <- renderText(
      code_curvefit(input$time, input$resp, input$samples, input$model))

    estimate <- reactive({
      req(input$estimate); ds <- dataset(); req(input$resp, input$time)
      validate(need(length(input$groups) >= 1,
        "Pick at least one grouping factor above to estimate the best model per group."))
      em <- estimate_models(ds, input$time, input$resp, input$groups)
      num <- setdiff(names(em), c("sample", "n", "best_model"))
      em[num] <- lapply(em[num], round, 3); em$n <- NULL
      names(em)[names(em) == "sample"] <- "group"
      em
    })
    output$estimate_tab <- DT::renderDT(DT::datatable(estimate(), rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE)))
    output$estimate_dl <- csv_dl(estimate, "curve_model_estimation.csv")

    # ---- Fit-plot ----
    output$sample_ui <- renderUI({
      s <- fit()$table$sample
      selectInput(ns("sample"), "Sample", s, selected = s[1])
    })
    # sortable columns: fit parameters + the sample-ID and grouping columns
    output$fp_sort_ui <- renderUI({
      res <- fit()
      choices <- c("R² (fit quality)" = "r2", "DELTA (rate)" = "DELTA",
        "INTERCEPT" = "INTERCEPT", "Sample ID" = "sample")
      extra <- setdiff(c(res$samples, res$groups), "sample")
      choices <- c(choices, stats::setNames(extra, extra))
      selectInput(ns("fp_sort"), "Sort panels by", choices, selected = "r2")
    })
    # samples in the chosen order (worst-R² first by default)
    fp_ordered <- reactive({
      tab <- fit()$table
      sb <- input$fp_sort; if (is.null(sb) || !sb %in% names(tab)) sb <- "r2"
      ord <- order(tab[[sb]], decreasing = identical(input$fp_sortdir, "desc"))
      tab$sample[ord]
    })
    fp_per_page <- reactive({
      max(1, (input$fp_ncol %||% 3) * (input$fp_nrow %||% 3))
    })
    fp_npages <- reactive(max(1, ceiling(length(fp_ordered()) / fp_per_page())))
    output$fp_page_ui <- renderUI({
      np <- fp_npages()
      if (np <= 1) return(helpText(sprintf("%d sample(s) — one page.",
                                           length(fp_ordered()))))
      sliderInput(ns("fp_page"), "Page", min = 1, max = np, value = 1, step = 1)
    })

    output$fit_plot <- renderPlot({
      res <- fit(); ds <- dataset()
      key <- do.call(paste, c(ds$data[res$samples], sep = " | "))
      r2v <- setNames(res$table$r2, res$table$sample)

      if (identical(input$fp_view, "single")) {
        req(input$sample)
        sub <- ds$data[key == input$sample, , drop = FALSE]
        x <- sub[[res$time]]; y <- sub[[res$response]]
        ln <- curve_fit_line(x, y, res$model)
        p <- ggplot2::ggplot(data.frame(x = x, y = y), ggplot2::aes(x, y)) +
          ggplot2::geom_point(size = 2.5, alpha = 0.85, colour = "#08519C") +
          ggplot2::labs(x = res$time, y = res$response,
            title = sprintf("%s   (%s model, R² = %.3f)", input$sample,
                            res$model, r2v[[input$sample]])) +
          ggplot2::theme_minimal(base_size = 13)
        if (!is.null(ln))
          p <- p + ggplot2::geom_line(data = ln, ggplot2::aes(x, y),
                                      colour = "#E45756", linewidth = 1)
        return(p)
      }

      # grid: one page of panels, in the chosen sort order
      all_s <- fp_ordered(); per <- fp_per_page()
      pg <- min(max(1, input$fp_page %||% 1), fp_npages())
      idx <- seq((pg - 1) * per + 1, min(pg * per, length(all_s)))
      samps <- all_s[idx]
      keep <- key %in% samps
      pts <- data.frame(.sample = factor(key[keep], levels = samps),
                        x = ds$data[[res$time]][keep],
                        y = ds$data[[res$response]][keep])
      lns <- do.call(rbind, lapply(samps, function(s) {
        sub <- ds$data[key == s, , drop = FALSE]
        ln <- curve_fit_line(sub[[res$time]], sub[[res$response]], res$model)
        if (is.null(ln)) return(NULL)
        data.frame(.sample = factor(s, levels = samps), x = ln$x, y = ln$y)
      }))
      labs_df <- data.frame(.sample = factor(samps, levels = samps),
        lab = sprintf("R²=%.2f", r2v[samps]))
      p <- ggplot2::ggplot(pts, ggplot2::aes(x, y)) +
        ggplot2::geom_point(size = 1.4, alpha = 0.8, colour = "#08519C") +
        ggplot2::facet_wrap(~.sample, ncol = input$fp_ncol %||% 3,
                            scales = "free") +
        ggplot2::geom_text(data = labs_df, inherit.aes = FALSE, size = 3,
          colour = "grey30", hjust = 0, vjust = 1.3,
          ggplot2::aes(x = -Inf, y = Inf, label = lab)) +
        ggplot2::labs(x = res$time, y = res$response,
          title = sprintf("%s model — page %d of %d (sorted by %s, %s)",
            res$model, pg, fp_npages(),
            if (is.null(input$fp_sort)) "r2" else input$fp_sort,
            if (identical(input$fp_sortdir, "desc")) "high→low" else "low→high")) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
      if (!is.null(lns))
        p <- p + ggplot2::geom_line(data = lns, ggplot2::aes(x, y),
                                    colour = "#E45756", linewidth = 0.8)
      p
    }, height = function() {
      if (identical(input$fp_view, "single")) return(460)
      max(320, (input$fp_nrow %||% 3) * 230)
    })
    output$legend_fp <- renderText({
      res <- fit()
      sprintf(paste0("Observed values of %s (points) against %s, with the fitted ",
        "%s model (red line) for each plant. R² is the fraction of the variance ",
        "in %s that the fitted growth curve explains: the closer R² is to 1, the ",
        "more of the variation the curve captures, so the better the points ",
        "follow the line. An R² near 0 means the model explains almost none of ",
        "the variance."),
        res$response, res$time, res$model, res$response)
    })

    # ---- Examine differences ----
    output$byfac_ui <- renderUI({
      g <- input$groups
      validate(need(length(g) >= 1,
        "Pick grouping factors in the Fit tab to compare the fitted parameters."))
      selectInput(ns("byfac"), "Compare across", g)
    })
    output$exfacet_ui <- renderUI({
      g <- setdiff(input$groups, input$byfac)
      selectInput(ns("exfac"), "Facet plot by (one or more)", g, multiple = TRUE)
    })
    ex_res <- reactive({
      t <- fit()$table; req(input$byfac, input$param)
      validate(need(input$byfac %in% names(t), "Choose a grouping factor used in the fit."))
      fv <- setdiff(input$exfac, input$byfac)
      cmp <- data.frame(cluster = factor(t[[input$byfac]]), value = t[[input$param]])
      if (length(fv))
        cmp$facet <- factor(do.call(paste, c(t[fv], sep = " | ")))
      cmp <- cmp[is.finite(cmp$value), ]
      faceted <- length(fv) > 0
      list(cmp = cmp, faceted = faceted, fvars = fv,
           res = if (!faceted) anova_tukey_letters(cmp, "value", "cluster") else NULL)
    })
    output$ex_plot <- renderPlot({
      x <- ex_res(); ylab <- sprintf("%s (%s model)", input$param, fit()$model)
      if (x$faceted)
        .compare_plot_faceted(x$cmp, ylab, input$ex_geom, input$ex_pal,
          xlab = input$byfac, facet_label = paste(x$fvars, collapse = " | "))
      else
        .compare_plot(x$cmp, ylab, x$res, input$ex_geom, input$ex_pal,
          xlab = input$byfac)
    })
    ex_tbl <- reactive({
      x <- ex_res()
      if (x$faceted) {
        s <- .cf_faceted_summary(x$cmp)
        validate(need(!is.null(s), "No panel had two comparable groups."))
        s[c("mean", "sd", "se")] <- lapply(s[c("mean", "sd", "se")], round, 4)
        fname <- paste(x$fvars, collapse = " | ")
        names(s)[names(s) == "group"] <- input$byfac
        names(s)[names(s) == "facet"] <- fname
        s[c(fname, input$byfac, "n", "mean", "sd", "se", "letter")]
      } else {
        s <- x$res$summary
        s[c("mean", "sd", "se")] <- lapply(s[c("mean", "sd", "se")], round, 4)
        names(s)[names(s) == "group"] <- input$byfac
        s[c(input$byfac, "n", "mean", "sd", "se", "letter")]
      }
    })
    output$ex_tab <- DT::renderDT(DT::datatable(ex_tbl(), rownames = FALSE,
      options = list(dom = "tp", pageLength = 15, scrollX = TRUE)))
    output$ex_dl <- csv_dl(ex_tbl, "curve_parameter_anova.csv")
    output$legend_ex <- renderText({
      x <- ex_res(); ylab <- sprintf("%s (%s model)", input$param, fit()$model)
      if (x$faceted)
        sprintf(paste0("The fitted %s is compared across %s, with a separate ",
          "one-way ANOVA + Tukey HSD run within each %s panel. Letters are ",
          "comparable only within a panel: groups sharing a letter are not ",
          "significantly different (p >= 0.05)."),
          ylab, input$byfac, paste(x$fvars, collapse = " | "))
      else
        legend_anova_tukey(ylab, input$byfac, anova_p = x$res$anova_p,
          n_groups = nrow(x$res$summary), geom = input$ex_geom)
    })
  })
}
