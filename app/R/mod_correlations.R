# Correlations module: full port of the legacy Correlations tab.
#   * Correlation Plot sub-tab: method, plot style, type, ordering, palette
#     (RColorBrewer + ggsci) + levels, non-significant cross-out, split-by-factor
#     upper/lower-triangle overlay, correlation table, PDF/CSV download.
#   * Scatterplots sub-tab: x/y, colour (ggsci D3), shape, regression line + CI, rug.
# All statistics come from the mvapp backend; the module only wires + renders.

.corr_num_phenos <- function(ds) {
  ds$roles$phenotypes[vapply(ds$data[ds$roles$phenotypes], is.numeric, logical(1))]
}
.corr_factors <- function(ds) {
  cand <- unique(c(ds$roles$group, ds$roles$id, ds$roles$time))
  cand[cand %in% names(ds$data)]
}

# Palette choices offered in the UI (RColorBrewer diverging + ggsci ramps).
.corr_palette_choices <- c(
  "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG",
  "ggsci: GSEA", "ggsci: D3", "ggsci: NPG", "ggsci: Lancet")

# Resolve a palette name to a smooth 200-colour ramp for corrplot's `col`.
# ggsci qualitative palettes are turned into diverging ramps (pole-white-pole);
# ggsci GSEA is already diverging.
.corr_palette <- function(name, n = 10) {
  n <- max(3L, min(as.integer(n), 11L))
  cols <- if (startsWith(name, "ggsci:")) {
    key <- trimws(sub("ggsci:", "", name))
    switch(key,
      "GSEA"   = ggsci::pal_gsea("default")(12),
      "D3"     = { p <- ggsci::pal_d3("category10")(10); c(p[1], "white", p[2]) },
      "NPG"    = { p <- ggsci::pal_npg("nrc")(10);       c(p[2], "white", p[1]) },
      "Lancet" = { p <- ggsci::pal_lancet("lanonc")(9);  c(p[2], "white", p[1]) },
      ggsci::pal_gsea("default")(12))
  } else {
    RColorBrewer::brewer.pal(n, name)
  }
  grDevices::colorRampPalette(cols)(200)
}

# Categorical colour palettes for the scatterplot points (ggsci journals +
# RColorBrewer qualitative). Each maps to a discrete ggplot colour scale.
.scatter_palette_choices <- c(
  "ggsci: D3", "ggsci: NPG", "ggsci: Lancet", "ggsci: JAMA", "ggsci: NEJM",
  "ggsci: AAAS", "ggsci: JCO", "ggsci: IGV", "ggsci: Simpsons",
  "Brewer: Set1", "Brewer: Dark2", "Brewer: Paired", "Brewer: Set2",
  "Default (ggplot)")

# Max number of distinct colours each palette can supply.
.scatter_palette_capacity <- c(
  "ggsci: D3" = 10, "ggsci: NPG" = 10, "ggsci: Lancet" = 9, "ggsci: JAMA" = 7,
  "ggsci: NEJM" = 8, "ggsci: AAAS" = 10, "ggsci: JCO" = 10, "ggsci: IGV" = 51,
  "ggsci: Simpsons" = 16, "Brewer: Set1" = 9, "Brewer: Dark2" = 8,
  "Brewer: Paired" = 12, "Brewer: Set2" = 8)

# Resolve a palette name + number of levels to a ggplot colour scale, or NULL
# (ggplot's default hue) when the palette can't supply enough distinct colours.
.scatter_colour_scale <- function(name, n) {
  cap <- .scatter_palette_capacity[name]
  if (!is.na(cap) && n > cap) return(NULL)   # graceful fallback to default hue
  if (startsWith(name, "ggsci:")) {
    key <- trimws(sub("ggsci:", "", name))
    return(switch(key,
      "D3"       = ggsci::scale_colour_d3(),
      "NPG"      = ggsci::scale_colour_npg(),
      "Lancet"   = ggsci::scale_colour_lancet(),
      "JAMA"     = ggsci::scale_colour_jama(),
      "NEJM"     = ggsci::scale_colour_nejm(),
      "AAAS"     = ggsci::scale_colour_aaas(),
      "JCO"      = ggsci::scale_colour_jco(),
      "IGV"      = ggsci::scale_colour_igv(),
      "Simpsons" = ggsci::scale_colour_simpsons(),
      ggsci::scale_colour_d3()))
  }
  if (startsWith(name, "Brewer:")) {
    return(ggplot2::scale_colour_brewer(palette = trimws(sub("Brewer:", "", name))))
  }
  NULL
}

mod_correlations_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    # ------------------------------------------------------------ Correlation
    tabPanel("Correlation Plot",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("phenos_ui")),
          checkboxInput(ns("subset"), "Subset the data?", FALSE),
          uiOutput(ns("subset_ui")),
          tags$hr(),
          selectInput(ns("method"), "Correlation method",
                      c("pearson", "spearman")),
          selectInput(ns("plot_method"), "Plot style",
                      c("circle", "square", "ellipse", "number", "shade", "pie")),
          selectInput(ns("type"), "Plot type", c("full", "lower", "upper")),
          selectInput(ns("order"), "Order labels by",
                      c("Original order" = "original",
                        "Angular order of eigenvectors" = "AOE",
                        "First principal component" = "FPC",
                        "Hierarchical clustering" = "hclust",
                        "Alphabetical" = "alphabet")),
          numericInput(ns("levels"), "Number of colour levels", 10, 3, 11),
          selectInput(ns("palette"), "Palette", .corr_palette_choices),
          checkboxInput(ns("sig"), "Cross out non-significant correlations", FALSE),
          conditionalPanel("input.sig == true", ns = ns,
            numericInput(ns("sig_level"), "Significance threshold",
                         0.05, 0.001, 0.5, step = 0.01)),
          tags$hr(),
          checkboxInput(ns("split"),
                        "Split into upper/lower triangle by a factor", FALSE),
          conditionalPanel("input.split == true", ns = ns,
            uiOutput(ns("split_ui")),
            helpText("In split mode the two triangles show your two levels.",
                     "Ordering and the significance cross-out both apply (a",
                     "shared trait order keeps the triangles aligned, and each",
                     "triangle is crossed out by its own level's significance).",
                     "'Plot type' is set by the split and is ignored."))
        ),
        mainPanel(
          textOutput(ns("nobs")),
          checkboxInput(ns("show_legend"), "Show figure legend"),
          conditionalPanel("input.show_legend == true", ns = ns,
            textOutput(ns("legend"))),
          plotOutput(ns("corrplot"), height = "560px"),
          downloadButton(ns("dl_plot"), "Download plot (PDF)"),
          checkboxInput(ns("show_code"), "Show me the R code"),
          conditionalPanel("input.show_code == true", ns = ns,
            verbatimTextOutput(ns("code"))),
          tags$hr(),
          h5("Correlation table"),
          DT::DTOutput(ns("cortable")),
          downloadButton(ns("dl_table"), "Download table (CSV)")
        )
      )
    ),
    # ------------------------------------------------------------ Scatterplots
    tabPanel("Scatterplots",
      sidebarLayout(
        sidebarPanel(
          uiOutput(ns("sx_ui")),
          uiOutput(ns("sy_ui")),
          uiOutput(ns("colour_ui")),
          selectInput(ns("scatter_pal"), "Colour palette",
                      .scatter_palette_choices, selected = "ggsci: D3"),
          checkboxInput(ns("use_shape"), "Map shape to a factor?", FALSE),
          uiOutput(ns("shape_ui")),
          tags$hr(),
          checkboxInput(ns("lm"), "Add regression line", TRUE),
          conditionalPanel("input.lm == true", ns = ns,
            checkboxInput(ns("lm_ci"), "Show confidence interval", TRUE)),
          checkboxInput(ns("rug"), "Add marginal rug", FALSE)
        ),
        mainPanel(
          checkboxInput(ns("sc_show_legend"), "Show figure legend"),
          conditionalPanel("input.sc_show_legend == true", ns = ns,
            textOutput(ns("scatter_legend"))),
          plotly::plotlyOutput(ns("scatter"), height = "500px"),
          checkboxInput(ns("sc_show_code"), "Show me the R code"),
          conditionalPanel("input.sc_show_code == true", ns = ns,
            verbatimTextOutput(ns("scatter_code"))),
          tags$hr(),
          h5("Regression fit"),
          DT::DTOutput(ns("scatter_stats")),
          dl_btn(ns("scatter_stats_dl")))
      )
    )
  )
}

mod_correlations_server <- function(id, dataset, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    logmsg <- function(m) if (is.function(log_add)) log_add(m)

    # ---- shared inputs ----
    output$phenos_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("phenos"), "Phenotypes", choices = num,
                  selected = num, multiple = TRUE)
    })
    output$subset_ui <- renderUI({
      ds <- dataset(); req(ds, input$subset)
      facs <- .corr_factors(ds); req(length(facs) > 0)
      tagList(
        selectInput(ns("sub_factor"), "Subset factor", facs),
        uiOutput(ns("sub_value_ui")))
    })
    output$sub_value_ui <- renderUI({
      ds <- dataset(); req(ds, input$sub_factor)
      selectInput(ns("sub_value"), "Keep level",
                  sort(unique(as.character(ds$data[[input$sub_factor]]))))
    })

    # split controls (factors with >= 2 levels)
    output$split_ui <- renderUI({
      ds <- dataset(); req(ds)
      facs <- .corr_factors(ds)
      facs <- facs[vapply(facs,
        function(f) length(unique(ds$data[[f]])) >= 2, logical(1))]
      req(length(facs) > 0)
      sel <- if (length(ds$roles$group)) ds$roles$group[1] else facs[1]
      tagList(
        selectInput(ns("split_by"), "Split factor", facs, selected = sel),
        uiOutput(ns("split_levels_ui")))
    })
    output$split_levels_ui <- renderUI({
      ds <- dataset(); req(ds, input$split_by)
      levs <- sort(unique(as.character(ds$data[[input$split_by]])))
      req(length(levs) >= 2)
      tagList(
        selectInput(ns("upper_level"), "Upper triangle level", levs,
                    selected = levs[1]),
        selectInput(ns("lower_level"), "Lower triangle level", levs,
                    selected = levs[2]))
    })

    # data after optional subsetting
    sub_data <- reactive({
      ds <- dataset(); req(ds)
      d <- ds$data
      if (isTRUE(input$subset) && !is.null(input$sub_factor) &&
          !is.null(input$sub_value)) {
        d <- d[as.character(d[[input$sub_factor]]) == input$sub_value, ,
               drop = FALSE]
      }
      d
    })

    # ---- Correlation Plot ----
    corM <- reactive({
      req(input$phenos)
      validate(need(length(input$phenos) >= 2, "Select at least two phenotypes."))
      cor_matrix(sub_data(), input$phenos, method = input$method)
    })

    draw_corr <- function() {
      req(input$phenos)
      validate(need(length(input$phenos) >= 2, "Select at least two phenotypes."))
      col <- .corr_palette(input$palette, input$levels)

      if (isTRUE(input$split)) {
        req(input$split_by, input$upper_level, input$lower_level)
        validate(need(input$upper_level != input$lower_level,
                      "Choose two different levels for the two triangles."))
        sp <- cor_split(sub_data(), input$phenos, by = input$split_by,
                        upper_level = input$upper_level,
                        lower_level = input$lower_level, method = input$method)

        # One shared trait order for BOTH triangles (from the mean of the two
        # matrices) so the upper/lower overlay stays aligned.
        Mref <- (sp$upper + sp$lower) / 2
        idx <- if (input$order == "original") {
          seq_len(ncol(Mref))
        } else {
          corrplot::corrMatOrder(Mref, order = input$order,
                                 hclust.method = "complete")
        }
        Mu <- sp$upper[idx, idx]; Ml <- sp$lower[idx, idx]

        upper_args <- list(Mu, type = "upper", method = input$plot_method,
          order = "original", tl.col = "black", col = col, tl.pos = "lt",
          mar = c(0, 0, 2, 0),
          title = sprintf("upper: %s  (n=%d)   |   lower: %s  (n=%d)",
                          sp$upper_level, sp$n_upper, sp$lower_level, sp$n_lower))
        lower_args <- list(Ml, type = "lower", method = input$plot_method,
          order = "original", tl.col = "black", col = col, add = TRUE,
          tl.pos = "n", cl.pos = "n")

        if (isTRUE(input$sig)) {
          Pu <- sp$p_upper[idx, idx]; Pl <- sp$p_lower[idx, idx]
          crossout <- list(sig.level = input$sig_level, insig = "pch",
                           pch.col = "grey40")
          upper_args <- c(upper_args, list(p.mat = Pu), crossout)
          lower_args <- c(lower_args, list(p.mat = Pl), crossout)
        }
        do.call(corrplot::corrplot, upper_args)
        do.call(corrplot::corrplot, lower_args)
        return(invisible())
      }

      args <- list(corr = corM(), method = input$plot_method, type = input$type,
                   order = input$order, tl.col = "black", col = col)
      if (isTRUE(input$sig)) {
        P <- cor_pmatrix(sub_data(), input$phenos, method = input$method)
        args <- c(args, list(p.mat = P, sig.level = input$sig_level,
                             insig = "pch", pch.col = "grey40"))
      }
      do.call(corrplot::corrplot, args)
    }

    output$nobs <- renderText({
      req(input$phenos)
      paste0("n = ", cor_n(sub_data(), input$phenos),
             " complete observations across the selected traits.")
    })
    output$corrplot <- renderPlot(draw_corr())

    output$cortable <- DT::renderDT({
      DT::datatable(as.data.frame(round(corM(), 3)),
                    options = list(scrollX = TRUE))
    })

    output$dl_plot <- downloadHandler(
      filename = function() "MVApp_correlation_plot.pdf",
      content = function(f) {
        grDevices::pdf(f); draw_corr(); grDevices::dev.off()
        logmsg(sprintf("Downloaded correlation plot (%s, %s%s).", input$method,
                       input$plot_method,
                       if (isTRUE(input$split))
                         sprintf(", split by %s", input$split_by) else ""))
      })
    output$dl_table <- downloadHandler(
      filename = function() "MVApp_correlation_table.csv",
      content = function(f) {
        utils::write.csv(round(corM(), 4), f)
        logmsg("Downloaded correlation table (CSV).")
      })

    # --- Show me the code / figure legend (Correlation Plot) ---
    output$legend <- renderText({
      req(input$phenos)
      legend_correlation(input$phenos, method = input$method,
        n = cor_n(sub_data(), input$phenos), split = isTRUE(input$split),
        upper_level = input$upper_level, lower_level = input$lower_level,
        sig = isTRUE(input$sig), sig_level = input$sig_level)
    })
    output$code <- renderText({
      req(input$phenos)
      code_correlation(input$phenos, method = input$method,
        plot_style = input$plot_method, type = input$type, order = input$order,
        palette = input$palette, levels = input$levels, sig = isTRUE(input$sig),
        sig_level = input$sig_level, split = isTRUE(input$split),
        by = input$split_by, upper_level = input$upper_level,
        lower_level = input$lower_level)
    })

    # ---- Scatterplots ----
    output$sx_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("sx"), "X trait", num, selected = num[1])
    })
    output$sy_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("sy"), "Y trait", num,
                  selected = if (length(num) > 1) num[2] else num[1])
    })
    output$colour_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("colour"), "Colour by", c("<none>", .corr_factors(ds)),
                  selected = if (length(ds$roles$group)) ds$roles$group[1] else "<none>")
    })
    output$shape_ui <- renderUI({
      ds <- dataset(); req(ds, input$use_shape)
      selectInput(ns("shape"), "Shape by", .corr_factors(ds))
    })

    # sample-name column for the plotly tooltip (first id role, else first factor)
    id_col <- reactive({
      ds <- dataset(); req(ds)
      if (length(ds$roles$id)) ds$roles$id[1]
      else { f <- .corr_factors(ds); if (length(f)) f[1] else NULL }
    })

    output$scatter <- plotly::renderPlotly({
      d <- sub_data(); req(input$sx, input$sy)
      colour_on <- !is.null(input$colour) && input$colour != "<none>"
      shape_on <- isTRUE(input$use_shape) && !is.null(input$shape)
      idc <- id_col()
      # colour/shape are categorical: coerce to factor so numeric factors
      # (e.g. DAY) map to a discrete scale instead of a continuous one
      if (colour_on) d[[input$colour]] <- factor(d[[input$colour]])
      if (shape_on)  d[[input$shape]]  <- factor(d[[input$shape]])

      aes_args <- list(x = as.name(input$sx), y = as.name(input$sy))
      if (!is.null(idc) && idc %in% names(d))
        aes_args$text <- as.name(idc)                 # -> hover shows sample name
      if (colour_on) aes_args$colour <- as.name(input$colour)
      if (shape_on)  aes_args$shape <- as.name(input$shape)

      p <- ggplot2::ggplot(d, do.call(ggplot2::aes, aes_args)) +
        ggplot2::geom_point(size = 2, alpha = 0.85) +
        ggplot2::theme_minimal(base_size = 13)
      if (isTRUE(input$lm))
        p <- p + ggplot2::geom_smooth(ggplot2::aes(text = NULL), method = "lm",
                                      se = isTRUE(input$lm_ci), formula = y ~ x)
      if (isTRUE(input$rug)) p <- p + ggplot2::geom_rug(alpha = 0.4)
      if (colour_on) {                                 # chosen categorical palette
        sc <- .scatter_colour_scale(input$scatter_pal,
                                    length(unique(d[[input$colour]])))
        if (!is.null(sc)) p <- p + sc
      }

      plotly::ggplotly(p, tooltip = c("text", "x", "y", "colour", "shape"))
    })

    # per-group (or overall) regression stats, shared by table + legend
    scatter_stat_tbl <- reactive({
      req(input$sx, input$sy)
      scatter_stats(sub_data(), input$sx, input$sy, group = scatter_colour())
    })
    output$scatter_stats <- DT::renderDT({
      st <- scatter_stat_tbl()
      st$r  <- round(st$r, 3); st$R2 <- round(st$R2, 3); st$p <- signif(st$p, 3)
      names(st)[names(st) == "R2"] <- "R²"
      DT::datatable(st, rownames = FALSE, options = list(dom = "t"),
                    caption = paste0(input$sy, " ~ ", input$sx))
    })
    output$scatter_stats_dl <- csv_dl(scatter_stat_tbl,
                                      "scatter_regression_fit.csv")

    # --- Show me the code / figure legend (Scatterplots) ---
    scatter_colour <- reactive(
      if (!is.null(input$colour) && input$colour != "<none>") input$colour else NULL)

    output$scatter_legend <- renderText({
      st <- scatter_stat_tbl(); col <- scatter_colour()
      grp_r2 <- if (!is.null(col)) stats::setNames(st$R2, st$group) else NULL
      grp_p  <- if (!is.null(col)) stats::setNames(st$p, st$group) else NULL
      legend_scatter(input$sx, input$sy, colour = col, lm = isTRUE(input$lm),
                     ci = isTRUE(input$lm_ci), r2 = st$R2[1], p = st$p[1],
                     group_r2 = grp_r2, group_p = grp_p)
    })
    output$scatter_code <- renderText({
      req(input$sx, input$sy)
      code_scatter(input$sx, input$sy, colour = scatter_colour(),
                   lm = isTRUE(input$lm), ci = isTRUE(input$lm_ci),
                   rug = isTRUE(input$rug))
    })
  })
}
