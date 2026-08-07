# Data Exploration module: normality, homogeneity of variance, t-test, one-way
# ANOVA (Tukey letters), and two-way ANOVA. Reuses .corr_num_phenos /
# .corr_factors (mod_correlations.R), .compare_plot / .compare_geoms
# (mod_clustering.R), and cat_palette_choices / legend_ui / csv_dl (util.R).

# factors with exactly two levels (for the t-test)
.explore_two_level <- function(ds) {
  f <- .corr_factors(ds)
  f[vapply(f, function(x) length(unique(ds$data[[x]])) == 2, logical(1))]
}

mod_explore_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    # ------------------------------------------------------------- Normality
    tabPanel("Normality",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("nrm_trait_ui")),
          checkboxInput(ns("nrm_bygroup"), "Test per group?", FALSE),
          uiOutput(ns("nrm_group_ui")),
          uiOutput(ns("nrm_facet_ui"))),
        mainPanel(
          legend_ui(ns, "nrm"),
          fluidRow(column(6, plotOutput(ns("nrm_hist"), height = "320px")),
                   column(6, plotOutput(ns("nrm_qq"), height = "320px"))),
          h5("Shapiro-Wilk test"),
          DT::DTOutput(ns("nrm_tab")), dl_btn(ns("nrm_dl")),
          checkboxInput(ns("nrm_code_chk"), "Show me the R code"),
          conditionalPanel("input.nrm_code_chk == true", ns = ns,
            verbatimTextOutput(ns("nrm_code")))))),
    # -------------------------------------------------------- Equal variance
    tabPanel("Equal variance",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("var_trait_ui")),
          uiOutput(ns("var_group_ui")),
          selectInput(ns("var_method"), "Test",
            c("Levene (Brown-Forsythe)" = "levene", "Bartlett" = "bartlett")),
          uiOutput(ns("var_facet_ui"))),
        mainPanel(
          legend_ui(ns, "var"),
          plotOutput(ns("var_plot"), height = "340px"),
          h5("Test result(s)"),
          DT::DTOutput(ns("var_tab")), dl_btn(ns("var_dl")),
          checkboxInput(ns("var_code_chk"), "Show me the R code"),
          conditionalPanel("input.var_code_chk == true", ns = ns,
            verbatimTextOutput(ns("var_code")))))),
    # ---------------------------------------------------------------- t-test
    tabPanel("t-test",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("tt_trait_ui")),
          uiOutput(ns("tt_group_ui")),
          checkboxInput(ns("tt_var_equal"), "Assume equal variances (Student)?", FALSE),
          checkboxInput(ns("tt_paired"), "Paired?", FALSE),
          conditionalPanel("input.tt_paired == true", ns = ns,
            uiOutput(ns("tt_pairby_ui")),
            helpText("Values are averaged per pairing level x group, then paired ",
                     "across the pairing variable (e.g. one Control/Salt pair per ",
                     "genotype).")),
          uiOutput(ns("tt_facet_ui")),
          checkboxInput(ns("tt_screen_chk"),
            "List traits with significant differences", FALSE),
          helpText("The grouping factor must have exactly two levels.")),
        mainPanel(
          conditionalPanel("input.tt_screen_chk == true", ns = ns,
            h5("Traits differing between the two groups (t-test, p < 0.05)"),
            DT::DTOutput(ns("tt_screen_tab")), dl_btn(ns("tt_screen_dl")),
            tags$hr()),
          legend_ui(ns, "tt"),
          plotOutput(ns("tt_plot"), height = "360px"),
          h5("Test result(s)"),
          DT::DTOutput(ns("tt_tab")), dl_btn(ns("tt_dl")),
          checkboxInput(ns("tt_code_chk"), "Show me the R code"),
          conditionalPanel("input.tt_code_chk == true", ns = ns,
            verbatimTextOutput(ns("tt_code")))))),
    # ----------------------------------------------------- One-way ANOVA + Tukey
    tabPanel("ANOVA (one-way)",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("aov_trait_ui")),
          uiOutput(ns("aov_group_ui")),
          selectInput(ns("aov_geom"), "Plot type", .compare_geoms),
          selectInput(ns("aov_pal"), "Colour palette", cat_palette_choices,
            selected = "ggsci: D3"),
          uiOutput(ns("aov_facet_ui")),
          checkboxInput(ns("aov_screen_chk"),
            "List traits with significant differences", FALSE)),
        mainPanel(
          conditionalPanel("input.aov_screen_chk == true", ns = ns,
            h5("Traits with significant group differences (one-way ANOVA)"),
            DT::DTOutput(ns("aov_screen_tab")), dl_btn(ns("aov_screen_dl")),
            tags$hr()),
          legend_ui(ns, "aov"),
          plotOutput(ns("aov_plot"), height = "440px"),
          checkboxInput(ns("aov_code_chk"), "Show me the R code"),
          conditionalPanel("input.aov_code_chk == true", ns = ns,
            verbatimTextOutput(ns("aov_code"))),
          h5("ANOVA p-value(s)"),
          DT::DTOutput(ns("aov_pvals")), dl_btn(ns("aov_pvals_dl")),
          h5("Group summary (Tukey letters, per facet)"),
          DT::DTOutput(ns("aov_tab")), dl_btn(ns("aov_dl"))))),
    # ------------------------------------------------------- Non-parametric
    tabPanel("Non-parametric",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("np_trait_ui")),
          uiOutput(ns("np_group_ui")),
          selectInput(ns("np_geom"), "Plot type", .compare_geoms),
          selectInput(ns("np_pal"), "Colour palette", cat_palette_choices,
            selected = "ggsci: D3"),
          uiOutput(ns("np_facet_ui")),
          checkboxInput(ns("np_screen_chk"),
            "List traits with significant differences", FALSE),
          helpText("Use these rank-based tests when the data are not normally ",
                   "distributed or have unequal variances (check the Normality ",
                   "and Equal variance tabs first).")),
        mainPanel(
          conditionalPanel("input.np_screen_chk == true", ns = ns,
            h5("Traits with significant group differences (Kruskal-Wallis, p < 0.05)"),
            DT::DTOutput(ns("np_screen_tab")), dl_btn(ns("np_screen_dl")),
            tags$hr()),
          wellPanel(
            strong("When to use this: "),
            "Kruskal-Wallis is the non-parametric analogue of one-way ANOVA, and ",
            "pairwise Wilcoxon (Mann-Whitney) tests provide the post-hoc letters. ",
            "They rank the data instead of using the raw values, so they do ", em("not"),
            " assume normality or equal variances - a robust fallback when those ",
            "ANOVA assumptions fail. With two groups this is the Mann-Whitney U test."),
          legend_ui(ns, "np"),
          plotOutput(ns("np_plot"), height = "440px"),
          checkboxInput(ns("np_code_chk"), "Show me the R code"),
          conditionalPanel("input.np_code_chk == true", ns = ns,
            verbatimTextOutput(ns("np_code"))),
          h5("Kruskal-Wallis p-value(s)"),
          DT::DTOutput(ns("np_pvals")), dl_btn(ns("np_pvals_dl")),
          h5("Group summary (Wilcoxon letters, per facet)"),
          DT::DTOutput(ns("np_tab")), dl_btn(ns("np_dl"))))),
    # --------------------------------------------------------- Two-way ANOVA
    tabPanel("ANOVA (two-way)",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("a2_trait_ui")),
          uiOutput(ns("a2_fa_ui")),
          uiOutput(ns("a2_fb_ui")),
          checkboxInput(ns("a2_screen_chk"),
            "List traits with significant effects", FALSE)),
        mainPanel(
          conditionalPanel("input.a2_screen_chk == true", ns = ns,
            h5("Traits with significant main effects or interaction (two-way ANOVA)"),
            DT::DTOutput(ns("a2_screen_tab")), dl_btn(ns("a2_screen_dl")),
            tags$hr()),
          legend_ui(ns, "a2"),
          plotOutput(ns("a2_plot"), height = "380px"),
          h5("ANOVA table"),
          DT::DTOutput(ns("a2_tab")), dl_btn(ns("a2_dl")),
          checkboxInput(ns("a2_code_chk"), "Show me the R code"),
          conditionalPanel("input.a2_code_chk == true", ns = ns,
            verbatimTextOutput(ns("a2_code"))))))
  )
}

mod_explore_server <- function(id, dataset, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    logmsg <- function(m) if (is.function(log_add)) log_add(m)
    numsel <- function(inputId, lab = "Trait") renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns(inputId), lab, .corr_num_phenos(ds))
    })
    facsel <- function(inputId, lab, choices_fun = .corr_factors, sel = 1) renderUI({
      ds <- dataset(); req(ds); ch <- choices_fun(ds); req(length(ch) > 0)
      selectInput(ns(inputId), lab, ch, selected = ch[min(sel, length(ch))])
    })
    facetsel <- function(inputId) renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns(inputId), "Facet plot by (one or more)", .corr_factors(ds),
                  multiple = TRUE)
    })
    facet_vars <- function(v) if (is.null(v) || !length(v)) NULL else v
    build_facet <- function(data, vars) {
      if (is.null(vars) || !length(vars)) return(NULL)
      factor(do.call(paste, c(data[vars], sep = " | ")))
    }
    facet_lab <- function(vars) paste(vars, collapse = " | ")

    # ---- Normality ----
    output$nrm_trait_ui <- numsel("nrm_trait")
    output$nrm_group_ui <- renderUI({
      req(input$nrm_bygroup); ds <- dataset(); req(ds)
      selectInput(ns("nrm_group"), "Group", .corr_factors(ds))
    })
    output$nrm_facet_ui <- facetsel("nrm_facet")
    # variables that define each independently-tested subset
    nrm_split <- reactive(c(if (isTRUE(input$nrm_bygroup)) input$nrm_group,
                            facet_vars(input$nrm_facet)))
    nrm_data <- reactive({
      ds <- dataset(); req(ds, input$nrm_trait)
      sv <- nrm_split()
      d <- data.frame(value = ds$data[[input$nrm_trait]])
      if (length(sv)) d$facet <- build_facet(ds$data, sv)
      d
    })
    output$nrm_hist <- renderPlot({
      d <- nrm_data()
      p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
          bins = 30, fill = "#9ECAE1", colour = "white") +
        ggplot2::geom_density(colour = "#08519C", linewidth = 1) +
        ggplot2::labs(x = input$nrm_trait, title = "Distribution") +
        ggplot2::theme_minimal(base_size = 13)
      if ("facet" %in% names(d)) p <- p + ggplot2::facet_wrap(~facet, scales = "free")
      p
    })
    output$nrm_qq <- renderPlot({
      d <- nrm_data()
      p <- ggplot2::ggplot(d, ggplot2::aes(sample = value)) +
        ggplot2::stat_qq(alpha = 0.6) + ggplot2::stat_qq_line(colour = "red") +
        ggplot2::labs(title = "Normal Q-Q plot") +
        ggplot2::theme_minimal(base_size = 13)
      if ("facet" %in% names(d)) p <- p + ggplot2::facet_wrap(~facet, scales = "free")
      p
    })
    nrm_tbl <- reactive({
      ds <- dataset(); req(input$nrm_trait); sv <- nrm_split()
      st <- if (length(sv)) {
        dd <- ds$data; dd$.subset <- as.character(build_facet(dd, sv))
        s <- shapiro_by_group(dd, input$nrm_trait, ".subset")
        names(s)[1] <- paste(sv, collapse = " | "); s
      } else shapiro_by_group(ds, input$nrm_trait)
      st$W <- round(st$W, 4); st$p <- signif(st$p, 3); st
    })
    output$nrm_tab <- DT::renderDT(DT::datatable(nrm_tbl(), rownames = FALSE,
      options = list(pageLength = 10)))
    output$nrm_dl <- csv_dl(nrm_tbl, "normality_shapiro.csv")
    output$legend_nrm <- renderText({
      sv <- nrm_split()
      legend_normality(input$nrm_trait,
        if (length(sv)) paste(sv, collapse = " & ") else NULL)
    })
    output$nrm_code <- renderText(code_normality(input$nrm_trait,
      if (length(nrm_split())) nrm_split()[1] else NULL))

    # ---- Equal variance ----
    output$var_trait_ui <- numsel("var_trait")
    output$var_group_ui <- facsel("var_group", "Group")
    output$var_facet_ui <- facetsel("var_facet")
    output$var_plot <- renderPlot({
      ds <- dataset(); req(input$var_trait, input$var_group)
      fv <- facet_vars(input$var_facet)
      d <- ds$data[, c(input$var_trait, input$var_group, fv), drop = FALSE]
      names(d)[1:2] <- c("value", "grp"); d$grp <- factor(d$grp)
      fac <- build_facet(d, fv)
      # normalise each point by the mean of its group (within a facet cell), so
      # the plot shows *relative spread* (variance) on a common scale around 1
      cell <- if (is.null(fac)) d$grp else interaction(d$grp, fac, drop = TRUE)
      d$rel <- d$value / ave(d$value, cell, FUN = function(z) mean(z, na.rm = TRUE))
      if (!is.null(fac)) d$facet <- fac
      p <- ggplot2::ggplot(d, ggplot2::aes(grp, rel, fill = grp)) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA, colour = "grey30") +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
        ggsci::scale_fill_d3() +
        ggplot2::labs(x = input$var_group,
          y = sprintf("%s / group mean", input$var_trait),
          title = "Relative spread within each group (points ÷ group mean)") +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "none")
      if (!is.null(fac)) p <- p + ggplot2::facet_wrap(~facet)
      p
    })
    var_res_tbl <- reactive({
      ds <- dataset(); req(input$var_trait, input$var_group)
      fv <- facet_vars(input$var_facet)
      run1 <- function(dsub, label = NA) {
        r <- variance_test(dsub, input$var_trait, input$var_group,
                           method = input$var_method)
        data.frame(subset = label, method = r$method,
                   statistic = signif(r$statistic, 4), df = as.character(r$df),
                   p = signif(r$p, 4),
                   verdict = ifelse(r$p < 0.05, "variances differ (p<0.05)",
                                    "no difference"))
      }
      if (is.null(fv)) {
        out <- run1(ds$data); out$subset <- NULL; out
      } else {
        dd <- ds$data; dd$.f <- as.character(build_facet(dd, fv))
        out <- do.call(rbind, lapply(split(dd, dd$.f), function(sub)
          tryCatch(run1(sub, sub$.f[1]), error = function(e) NULL)))
        names(out)[names(out) == "subset"] <- paste(fv, collapse = " | "); out
      }
    })
    output$var_tab <- DT::renderDT(DT::datatable(var_res_tbl(), rownames = FALSE,
      options = list(dom = "t")))
    output$var_dl <- csv_dl(var_res_tbl, "variance_test.csv")
    output$legend_var <- renderText({
      r <- var_res_tbl()
      legend_variance(input$var_trait, input$var_group,
        method = sub(" .*", "", r$method[1]),
        p = if (nrow(r) == 1) r$p[1] else NA)
    })
    output$var_code <- renderText(
      code_variance(input$var_trait, input$var_group, input$var_method))

    # ---- t-test ----
    output$tt_trait_ui <- numsel("tt_trait")
    output$tt_group_ui <- facsel("tt_group", "Group (2 levels)", .explore_two_level)
    output$tt_pairby_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("tt_pairby"), "Pair by", .corr_factors(ds))
    })
    output$tt_facet_ui <- facetsel("tt_facet")
    # test result(s): one row overall, or one row per facet subset
    tt_res_tbl <- reactive({
      ds <- dataset(); req(input$tt_trait, input$tt_group)
      fv <- facet_vars(input$tt_facet)
      paired <- isTRUE(input$tt_paired)
      run1 <- function(dsub, label = NA) {
        r <- t_test_two(dsub, input$tt_trait, input$tt_group,
                        var_equal = isTRUE(input$tt_var_equal), paired = paired,
                        pair_by = if (paired) input$tt_pairby else NULL)
        data.frame(subset = label,
                   group1 = r$groups[1], mean1 = signif(r$mean1, 4),
                   group2 = r$groups[2], mean2 = signif(r$mean2, 4),
                   difference = signif(r$diff, 4),
                   CI_low = signif(r$conf_int[1], 4),
                   CI_high = signif(r$conf_int[2], 4),
                   t = signif(r$statistic, 4), df = signif(r$df, 4),
                   p = signif(r$p, 4), test = r$method)
      }
      if (is.null(fv)) {
        out <- run1(ds$data); out$subset <- NULL; out
      } else {
        dd <- ds$data; dd$.f <- as.character(build_facet(dd, fv))
        out <- do.call(rbind, lapply(split(dd, dd$.f), function(sub)
          tryCatch(run1(sub, sub$.f[1]), error = function(e) NULL)))
        names(out)[names(out) == "subset"] <- paste(fv, collapse = " | ")
        out
      }
    })
    output$tt_tab <- DT::renderDT(DT::datatable(tt_res_tbl(), rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)))
    output$tt_dl <- csv_dl(tt_res_tbl, "t_test_results.csv")
    output$tt_plot <- renderPlot({
      ds <- dataset(); req(input$tt_trait, input$tt_group)
      fv <- facet_vars(input$tt_facet)
      paired <- isTRUE(input$tt_paired) && !is.null(input$tt_pairby)
      cols <- c(input$tt_trait, input$tt_group,
                if (paired) input$tt_pairby, fv)
      d <- ds$data[stats::complete.cases(ds$data[cols]), cols, drop = FALSE]
      names(d)[1:2] <- c("value", "grp"); d$grp <- factor(d$grp)
      if (paired) names(d)[3] <- "pairid"
      fac <- build_facet(d, fv); if (!is.null(fac)) d$facet <- fac

      if (paired) {
        # one mean per pairing level x group (x facet), with connecting lines
        by_cols <- c("pairid", "grp", if (!is.null(fac)) "facet")
        agg <- stats::aggregate(value ~ ., data = d[c("value", by_cols)], FUN = mean)
        p <- ggplot2::ggplot(agg, ggplot2::aes(grp, value)) +
          ggplot2::geom_boxplot(ggplot2::aes(fill = grp), alpha = 0.5,
            outlier.shape = NA, colour = "grey30") +
          ggplot2::geom_line(ggplot2::aes(group = pairid), colour = "grey55",
            alpha = 0.6) +
          ggplot2::geom_point(ggplot2::aes(colour = grp), size = 2) +
          ggsci::scale_fill_d3() + ggsci::scale_colour_d3()
      } else {
        p <- ggplot2::ggplot(d, ggplot2::aes(grp, value, fill = grp)) +
          ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA, colour = "grey30") +
          ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
          ggsci::scale_fill_d3()
      }
      p <- p + ggplot2::labs(x = input$tt_group, y = input$tt_trait) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "none")
      if (!is.null(fac)) p <- p + ggplot2::facet_wrap(~facet)
      p
    })
    output$legend_tt <- renderText({
      tt <- tt_res_tbl()
      legend_ttest(input$tt_trait, input$tt_group,
        p = if (nrow(tt) == 1) tt$p[1] else NA,
        method = trimws(as.character(tt$test[1])))
    })
    output$tt_code <- renderText(code_ttest(input$tt_trait, input$tt_group,
      var_equal = isTRUE(input$tt_var_equal), paired = isTRUE(input$tt_paired),
      pair_by = if (isTRUE(input$tt_paired)) input$tt_pairby else NULL))
    tt_screen_df <- reactive({
      ds <- dataset(); req(input$tt_group)
      traits <- .corr_num_phenos(ds)
      fv <- facet_vars(input$tt_facet)
      if (!is.null(fv)) {
        dd <- ds$data; dd$.facet <- as.character(build_facet(dd, fv))
        sc <- ttest_screen_faceted(dd, traits, input$tt_group, ".facet")
        pcols <- setdiff(names(sc), "trait")
        sc <- sc[apply(sc[pcols], 1, function(p) any(p < 0.05, na.rm = TRUE)), ,
                 drop = FALSE]
        sc[pcols] <- lapply(sc[pcols], signif, 3)
      } else {
        sc <- ttest_screen(ds, traits, input$tt_group)
        sc <- sc[!is.na(sc$p) & sc$p < 0.05, , drop = FALSE]
        sc$p <- signif(sc$p, 3)
      }
      sc
    })
    output$tt_screen_tab <- DT::renderDT(DT::datatable(tt_screen_df(),
      rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE),
      caption = "One p-value column per facet subset (p < 0.05 in at least one)"))
    output$tt_screen_dl <- csv_dl(tt_screen_df, "ttest_significant_traits.csv")

    # ---- One-way ANOVA + Tukey ----
    output$aov_trait_ui <- numsel("aov_trait")
    output$aov_group_ui <- facsel("aov_group", "Group")
    output$aov_facet_ui <- facetsel("aov_facet")
    aov_res <- reactive({
      ds <- dataset(); req(ds, input$aov_trait, input$aov_group)
      cmp <- data.frame(cluster = factor(ds$data[[input$aov_group]]),
                        value = ds$data[[input$aov_trait]])
      list(cmp = cmp, res = anova_tukey_letters(cmp, "value", "cluster"))
    })
    output$aov_plot <- renderPlot({
      ds <- dataset(); req(input$aov_trait, input$aov_group)
      fv <- facet_vars(input$aov_facet)
      if (!is.null(fv)) {
        validate(need(!(input$aov_group %in% fv),
                      "Facet by variables other than the grouping factor."))
        cmp <- data.frame(cluster = factor(ds$data[[input$aov_group]]),
                          value = ds$data[[input$aov_trait]],
                          facet = build_facet(ds$data, fv))
        .compare_plot_faceted(cmp, input$aov_trait, input$aov_geom, input$aov_pal,
                              xlab = input$aov_group, facet_label = facet_lab(fv))
      } else {
        x <- aov_res()
        .compare_plot(x$cmp, input$aov_trait, x$res, input$aov_geom, input$aov_pal,
                      xlab = input$aov_group)
      }
    })
    # facet-aware summary: Tukey letters + ANOVA p computed WITHIN each facet
    aov_summary <- reactive({
      ds <- dataset(); req(input$aov_trait, input$aov_group)
      fv <- facet_vars(input$aov_facet)
      mk <- function(sub, fl = NA) {
        sub$cluster <- droplevels(factor(sub$cluster))
        r <- anova_tukey_letters(sub, "value", "cluster")
        s <- r$summary; s$anova_p <- r$anova_p
        if (!is.na(fl)) s$facet <- fl
        s
      }
      if (is.null(fv) || (input$aov_group %in% fv)) {
        cmp <- data.frame(cluster = factor(ds$data[[input$aov_group]]),
                          value = ds$data[[input$aov_trait]])
        list(faceted = FALSE, summary = mk(cmp))
      } else {
        cmp <- data.frame(cluster = factor(ds$data[[input$aov_group]]),
                          value = ds$data[[input$aov_trait]],
                          facet = build_facet(ds$data, fv))
        parts <- lapply(split(cmp, cmp$facet), function(sub)
          tryCatch(mk(sub, as.character(sub$facet[1])), error = function(e) NULL))
        list(faceted = TRUE, summary = do.call(rbind, parts),
             facet_name = facet_lab(fv))
      }
    })

    aov_pval_tbl <- reactive({
      a <- aov_summary(); s <- a$summary
      if (a$faceted) {
        agg <- do.call(rbind, lapply(split(s, s$facet), function(g)
          data.frame(facet = g$facet[1], n_groups = nrow(g),
                     anova_p = signif(g$anova_p[1], 3))))
        rownames(agg) <- NULL; names(agg)[1] <- a$facet_name; agg
      } else {
        data.frame(comparison = input$aov_group, n_groups = nrow(s),
                   anova_p = signif(s$anova_p[1], 3))
      }
    })
    output$aov_pvals <- DT::renderDT(DT::datatable(aov_pval_tbl(),
      rownames = FALSE, options = list(pageLength = 10)))
    output$aov_pvals_dl <- csv_dl(aov_pval_tbl, "anova_pvalues_per_facet.csv")

    aov_tbl <- reactive({
      a <- aov_summary(); s <- a$summary
      s[c("mean", "sd", "se")] <- lapply(s[c("mean", "sd", "se")], round, 3)
      cols <- c(if (a$faceted) "facet", "group", "n", "mean", "sd", "se", "letter")
      out <- s[cols]
      names(out)[names(out) == "group"] <- input$aov_group
      if (a$faceted) names(out)[names(out) == "facet"] <- a$facet_name
      out
    })
    output$aov_tab <- DT::renderDT(DT::datatable(aov_tbl(), rownames = FALSE,
      options = list(pageLength = 15)))
    output$aov_dl <- csv_dl(aov_tbl, "anova_tukey_summary.csv")

    aov_screen_df <- reactive({
      ds <- dataset(); req(input$aov_group)
      traits <- .corr_num_phenos(ds)
      fv <- facet_vars(input$aov_facet)
      if (!is.null(fv) && !(input$aov_group %in% fv)) {
        dd <- ds$data; dd$.facet <- as.character(build_facet(dd, fv))
        sc <- anova_screen_faceted(dd, traits, input$aov_group, ".facet")
        pcols <- setdiff(names(sc), "trait")
        sc <- sc[apply(sc[pcols], 1, function(p) any(p < 0.05, na.rm = TRUE)), ,
                 drop = FALSE]
        sc[pcols] <- lapply(sc[pcols], signif, 3)
      } else {
        sc <- anova_screen(ds, traits, input$aov_group)
        sc <- sc[!is.na(sc$anova_p) & sc$anova_p < 0.05, , drop = FALSE]
        sc$anova_p <- signif(sc$anova_p, 3)
      }
      sc
    })
    output$aov_screen_tab <- DT::renderDT(DT::datatable(aov_screen_df(),
      rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE),
      caption = "One p-value column per facet subset (p < 0.05 in at least one)"))
    output$aov_screen_dl <- csv_dl(aov_screen_df, "anova_significant_traits.csv")
    output$legend_aov <- renderText({
      x <- aov_res()
      legend_anova_tukey(input$aov_trait, input$aov_group, anova_p = x$res$anova_p,
        n_groups = nrow(x$res$summary), geom = input$aov_geom)
    })
    output$aov_code <- renderText(code_anova_tukey(input$aov_trait, input$aov_group))
    observeEvent(list(input$aov_trait, input$aov_group), {
      req(input$aov_trait, input$aov_group)
      logmsg(sprintf("One-way ANOVA: %s ~ %s (p = %.3g).", input$aov_trait,
        input$aov_group, aov_res()$res$anova_p))
    }, ignoreInit = TRUE)

    # ---- Non-parametric (Kruskal-Wallis + Wilcoxon letters) ----
    output$np_trait_ui <- numsel("np_trait")
    output$np_group_ui <- facsel("np_group", "Group")
    output$np_facet_ui <- facetsel("np_facet")

    # facet-aware summary: Kruskal-Wallis + Wilcoxon letters WITHIN each facet
    np_summary <- reactive({
      ds <- dataset(); req(input$np_trait, input$np_group)
      fv <- facet_vars(input$np_facet)
      mk <- function(sub, fl = NA) {
        sub$cluster <- droplevels(factor(sub$cluster))
        r <- kruskal_letters(sub, "value", "cluster")
        s <- r$summary; s$kruskal_p <- r$anova_p
        if (!is.na(fl)) s$facet <- fl
        s
      }
      if (is.null(fv) || (input$np_group %in% fv)) {
        cmp <- data.frame(cluster = factor(ds$data[[input$np_group]]),
                          value = ds$data[[input$np_trait]])
        list(faceted = FALSE, summary = mk(cmp))
      } else {
        cmp <- data.frame(cluster = factor(ds$data[[input$np_group]]),
                          value = ds$data[[input$np_trait]],
                          facet = build_facet(ds$data, fv))
        parts <- lapply(split(cmp, cmp$facet), function(sub)
          tryCatch(mk(sub, as.character(sub$facet[1])), error = function(e) NULL))
        list(faceted = TRUE, summary = do.call(rbind, parts),
             facet_name = facet_lab(fv))
      }
    })

    output$np_plot <- renderPlot({
      ds <- dataset(); req(input$np_trait, input$np_group)
      fv <- facet_vars(input$np_facet)
      if (!is.null(fv) && !(input$np_group %in% fv)) {
        cmp <- data.frame(cluster = factor(ds$data[[input$np_group]]),
                          value = ds$data[[input$np_trait]],
                          facet = build_facet(ds$data, fv))
        .compare_plot_faceted(cmp, input$np_trait, input$np_geom, input$np_pal,
          xlab = input$np_group, facet_label = facet_lab(fv),
          letter_fun = kruskal_letters,
          test_label = "Kruskal-Wallis + Wilcoxon test")
      } else {
        cmp <- data.frame(cluster = factor(ds$data[[input$np_group]]),
                          value = ds$data[[input$np_trait]])
        .compare_plot(cmp, input$np_trait,
          kruskal_letters(ds, input$np_trait, input$np_group),
          input$np_geom, input$np_pal, xlab = input$np_group,
          test_label = "Kruskal-Wallis")
      }
    })

    np_pval_tbl <- reactive({
      a <- np_summary(); s <- a$summary
      if (a$faceted) {
        agg <- do.call(rbind, lapply(split(s, s$facet), function(g)
          data.frame(facet = g$facet[1], n_groups = nrow(g),
                     kruskal_p = signif(g$kruskal_p[1], 3))))
        rownames(agg) <- NULL; names(agg)[1] <- a$facet_name; agg
      } else {
        data.frame(comparison = input$np_group, n_groups = nrow(s),
                   kruskal_p = signif(s$kruskal_p[1], 3))
      }
    })
    output$np_pvals <- DT::renderDT(DT::datatable(np_pval_tbl(),
      rownames = FALSE, options = list(pageLength = 10)))
    output$np_pvals_dl <- csv_dl(np_pval_tbl, "kruskal_pvalues_per_facet.csv")

    np_tbl <- reactive({
      a <- np_summary(); s <- a$summary
      s[c("mean", "sd", "se", "median")] <-
        lapply(s[c("mean", "sd", "se", "median")], round, 3)
      cols <- c(if (a$faceted) "facet", "group", "n", "median", "mean", "sd",
                "se", "letter")
      out <- s[cols]
      names(out)[names(out) == "group"] <- input$np_group
      if (a$faceted) names(out)[names(out) == "facet"] <- a$facet_name
      out
    })
    output$np_tab <- DT::renderDT(DT::datatable(np_tbl(), rownames = FALSE,
      options = list(pageLength = 15)))
    output$np_dl <- csv_dl(np_tbl, "kruskal_wilcoxon_summary.csv")
    output$legend_np <- renderText({
      a <- np_summary(); s <- a$summary
      legend_nonparam(input$np_trait, input$np_group,
        kruskal_p = if (a$faceted) NA else s$kruskal_p[1],
        n_groups = if (a$faceted) NA else nrow(s), geom = input$np_geom)
    })
    output$np_code <- renderText(code_kruskal(input$np_trait, input$np_group))
    np_screen_df <- reactive({
      ds <- dataset(); req(input$np_group)
      traits <- .corr_num_phenos(ds)
      fv <- facet_vars(input$np_facet)
      if (!is.null(fv) && !(input$np_group %in% fv)) {
        dd <- ds$data; dd$.facet <- as.character(build_facet(dd, fv))
        sc <- kruskal_screen_faceted(dd, traits, input$np_group, ".facet")
        pcols <- setdiff(names(sc), "trait")
        sc <- sc[apply(sc[pcols], 1, function(p) any(p < 0.05, na.rm = TRUE)), ,
                 drop = FALSE]
        sc[pcols] <- lapply(sc[pcols], signif, 3)
      } else {
        sc <- kruskal_screen(ds, traits, input$np_group)
        sc <- sc[!is.na(sc$kruskal_p) & sc$kruskal_p < 0.05, , drop = FALSE]
        sc$kruskal_p <- signif(sc$kruskal_p, 3)
      }
      sc
    })
    output$np_screen_tab <- DT::renderDT(DT::datatable(np_screen_df(),
      rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE),
      caption = "One p-value column per facet subset (p < 0.05 in at least one)"))
    output$np_screen_dl <- csv_dl(np_screen_df, "kruskal_significant_traits.csv")

    # ---- Two-way ANOVA ----
    output$a2_trait_ui <- numsel("a2_trait")
    output$a2_fa_ui <- facsel("a2_fa", "Factor A", .corr_factors, 1)
    output$a2_fb_ui <- facsel("a2_fb", "Factor B", .corr_factors, 2)
    a2_res <- reactive({
      ds <- dataset(); req(ds, input$a2_trait, input$a2_fa, input$a2_fb)
      validate(need(input$a2_fa != input$a2_fb, "Choose two different factors."))
      anova_two_way(ds, input$a2_trait, input$a2_fa, input$a2_fb)
    })
    output$a2_plot <- renderPlot({
      r <- a2_res(); m <- r$means
      names(m) <- c("A", "B", "mean"); m$A <- factor(m$A); m$B <- factor(m$B)
      ggplot2::ggplot(m, ggplot2::aes(A, mean, colour = B, group = B)) +
        ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2.5) +
        ggsci::scale_colour_d3() +
        ggplot2::labs(x = r$factor_a, colour = r$factor_b,
          y = sprintf("mean %s", input$a2_trait), title = "Interaction plot") +
        ggplot2::theme_minimal(base_size = 13)
    })
    a2_tbl <- reactive({
      t <- a2_res()$table
      t$Sum_Sq <- round(t$Sum_Sq, 2); t$Mean_Sq <- round(t$Mean_Sq, 2)
      t$F_value <- round(t$F_value, 3); t$p <- signif(t$p, 3); t
    })
    output$a2_tab <- DT::renderDT(DT::datatable(a2_tbl(), rownames = FALSE,
      options = list(dom = "t")))
    output$a2_dl <- csv_dl(a2_tbl, "two_way_anova.csv")
    output$legend_a2 <- renderText({
      t <- a2_res()$table
      gp <- function(term) { p <- t$p[t$term == term]; if (length(p)) p else NA }
      legend_anova2way(input$a2_trait, input$a2_fa, input$a2_fb,
        p_a = gp(input$a2_fa), p_b = gp(input$a2_fb),
        p_int = gp(paste0(input$a2_fa, ":", input$a2_fb)))
    })
    output$a2_code <- renderText(
      code_anova2way(input$a2_trait, input$a2_fa, input$a2_fb))
    a2_screen_df <- reactive({
      ds <- dataset(); req(input$a2_fa, input$a2_fb)
      validate(need(input$a2_fa != input$a2_fb, "Choose two different factors."))
      sc <- anova2_screen(ds, .corr_num_phenos(ds), input$a2_fa, input$a2_fb)
      sig <- sc[apply(sc[c("p_A", "p_B", "p_interaction")], 1,
        function(p) any(p < 0.05, na.rm = TRUE)), , drop = FALSE]
      sig[c("p_A", "p_B", "p_interaction")] <-
        lapply(sig[c("p_A", "p_B", "p_interaction")], signif, 3)
      names(sig) <- c("trait", paste0("p_", input$a2_fa),
                      paste0("p_", input$a2_fb), "p_interaction")
      sig
    })
    output$a2_screen_tab <- DT::renderDT(DT::datatable(a2_screen_df(),
      rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)))
    output$a2_screen_dl <- csv_dl(a2_screen_df, "two_way_significant_traits.csv")
  })
}
