# PCA module: thin UI over mvapp::run_pca().
#  * pre-PCA subsetting by a categorical variable
#  * individuals scatter (PCx vs PCy) coloured by a categorical variable
#  * scree, variable contributions, eigenvalues
# Reuses .corr_num_phenos / .corr_factors from mod_correlations.R.

mod_pca_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      uiOutput(ns("pheno_ui")),
      checkboxInput(ns("scale"), "Scale traits to unit variance", TRUE),
      checkboxInput(ns("subset"), "Subset the data before PCA?", FALSE),
      uiOutput(ns("subset_ui")),
      tags$hr(),
      uiOutput(ns("colour_ui")),
      fluidRow(
        column(6, numericInput(ns("pcx"), "PC on x-axis", 1, min = 1, step = 1)),
        column(6, numericInput(ns("pcy"), "PC on y-axis", 2, min = 1, step = 1))
      ),
      actionButton(ns("go"), "Run PCA", class = "btn-primary")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Individuals", legend_ui(ns, "ind"),
          plotly::plotlyOutput(ns("indplot"), height = "480px")),
        tabPanel("Scree", legend_ui(ns, "scree"),
          plotOutput(ns("scree"), height = "400px")),
        tabPanel("Contributions", legend_ui(ns, "contrib"),
          plotOutput(ns("varplot"), height = "480px")),
        tabPanel("Eigenvalues", legend_ui(ns, "eig"),
          DT::DTOutput(ns("eig")), dl_btn(ns("eig_dl")))
      ),
      checkboxInput(ns("show_code"), "Show me the R code"),
      conditionalPanel("input.show_code == true", ns = ns,
        verbatimTextOutput(ns("code"))),
      tags$hr(),
      h5("Sample coordinates"),
      DT::DTOutput(ns("ind_tab")), dl_btn(ns("ind_dl"))
    )
  )
}

mod_pca_server <- function(id, dataset, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    logmsg <- function(m) if (is.function(log_add)) log_add(m)

    output$pheno_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("phenos"), "Phenotypes", num, selected = num, multiple = TRUE)
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
    output$colour_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("colour"), "Colour individuals by",
        c("<none>", .corr_factors(ds)),
        selected = if (length(ds$roles$group)) ds$roles$group[1] else "<none>")
    })

    # dataset after optional subsetting (roles preserved)
    sub_ds <- reactive({
      ds <- dataset(); req(ds)
      d <- ds$data
      if (isTRUE(input$subset) && !is.null(input$sub_factor) &&
          !is.null(input$sub_value)) {
        d <- d[as.character(d[[input$sub_factor]]) == input$sub_value, , drop = FALSE]
      }
      mvapp_dataset(d, id = ds$roles$id, group = ds$roles$group,
                    time = ds$roles$time, phenotypes = ds$roles$phenotypes)
    })

    pca <- eventReactive(input$go, {
      ds <- sub_ds(); req(input$phenos)
      validate(need(length(input$phenos) >= 2, "Select at least two phenotypes."))
      res <- run_pca(ds, phenotypes = input$phenos, scale = input$scale)
      logmsg(sprintf("PCA on %d traits (%s), scale = %s%s; PC1/PC2 explain %.1f%%/%.1f%%.",
                     length(input$phenos), paste(input$phenos, collapse = ", "),
                     input$scale,
                     if (isTRUE(input$subset))
                       sprintf(", subset %s=%s", input$sub_factor, input$sub_value) else "",
                     res$eig$variance_percent[1], res$eig$variance_percent[2]))
      res
    })

    # sample coordinates with categorical columns attached (for colouring)
    ind_df <- reactive({
      res <- pca(); ds <- sub_ds()
      d <- res$ind
      idcols <- c(ds$roles$id, ds$roles$group, ds$roles$time)
      key <- make.unique(do.call(paste, c(ds$data[idcols], sep = "_")))
      idx <- match(d$id, key)
      for (f in .corr_factors(ds)) d[[f]] <- ds$data[[f]][idx]
      d
    })

    output$indplot <- plotly::renderPlotly({
      res <- pca(); d <- ind_df()
      cx <- paste0("Dim.", input$pcx); cy <- paste0("Dim.", input$pcy)
      validate(need(all(c(cx, cy) %in% names(d)), "Chosen PCs are out of range."))
      col_on <- !is.null(input$colour) && input$colour != "<none>"
      if (col_on) d[[input$colour]] <- factor(d[[input$colour]])  # discrete colour
      aes_args <- list(x = as.name(cx), y = as.name(cy), text = as.name("id"))
      if (col_on) aes_args$colour <- as.name(input$colour)
      vx <- res$eig$variance_percent[input$pcx]
      vy <- res$eig$variance_percent[input$pcy]
      p <- ggplot2::ggplot(d, do.call(ggplot2::aes, aes_args)) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
        ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
        ggplot2::labs(title = "PCA of individuals",
                      x = sprintf("PC%d (%.1f%%)", input$pcx, vx),
                      y = sprintf("PC%d (%.1f%%)", input$pcy, vy)) +
        ggplot2::theme_minimal(base_size = 13)
      if (col_on && length(unique(d[[input$colour]])) <= 10)
        p <- p + ggsci::scale_colour_d3()
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "colour"))
    })

    output$scree <- renderPlot({
      e <- pca()$eig
      ggplot2::ggplot(e, ggplot2::aes(factor(component), variance_percent)) +
        ggplot2::geom_col(fill = "#4C78A8") +
        ggplot2::geom_line(ggplot2::aes(group = 1), colour = "#E45756") +
        ggplot2::geom_point(colour = "#E45756") +
        ggplot2::labs(x = "Principal component", y = "% variance explained",
                      title = "Scree plot") +
        ggplot2::theme_minimal(base_size = 13)
    })

    output$varplot <- renderPlot({
      res <- pca()
      cx <- paste0("Dim.", input$pcx); cy <- paste0("Dim.", input$pcy)
      vc <- res$var_coord; ct <- res$contrib
      validate(need(all(c(cx, cy) %in% names(vc)), "Chosen PCs are out of range."))
      df <- data.frame(trait = vc$trait, x = vc[[cx]], y = vc[[cy]],
                       contrib = ct[[cx]] + ct[[cy]])
      vx <- res$eig$variance_percent[input$pcx]
      vy <- res$eig$variance_percent[input$pcy]
      ggplot2::ggplot(df) +
        ggplot2::annotate("path", colour = "grey80",
          x = cos(seq(0, 2 * pi, length.out = 100)),
          y = sin(seq(0, 2 * pi, length.out = 100))) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
        ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
        ggplot2::geom_segment(ggplot2::aes(0, 0, xend = x, yend = y, colour = contrib),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
        ggplot2::geom_text(ggplot2::aes(x, y, label = trait), vjust = -0.5, size = 4) +
        ggplot2::scale_colour_gradient(low = "#9ECAE1", high = "#08306B",
                                       name = "contribution") +
        ggplot2::coord_equal() +
        ggplot2::labs(title = "Variable contributions",
          x = sprintf("PC%d (%.1f%%)", input$pcx, vx),
          y = sprintf("PC%d (%.1f%%)", input$pcy, vy)) +
        ggplot2::theme_minimal(base_size = 13)
    })

    output$eig <- DT::renderDT({
      e <- pca()$eig
      e[] <- lapply(e, function(z) if (is.numeric(z)) round(z, 3) else z)
      DT::datatable(e, rownames = FALSE, options = list(dom = "t"))
    })
    output$eig_dl <- csv_dl(reactive(pca()$eig), "PCA_eigenvalues.csv")

    output$ind_tab <- DT::renderDT(
      DT::datatable(ind_df(), rownames = FALSE, options = list(scrollX = TRUE)))
    output$ind_dl <- csv_dl(ind_df, "PCA_sample_coordinates.csv")

    # per-sub-tab legends (each describes its own plot; PC-aware where relevant)
    output$legend_ind <- renderText({
      res <- pca()
      legend_pca(input$phenos, n = res$n, scale = input$scale, what = "individuals",
        pcx = input$pcx, pcy = input$pcy,
        varx = res$eig$variance_percent[input$pcx],
        vary = res$eig$variance_percent[input$pcy],
        colour = if (!is.null(input$colour) && input$colour != "<none>")
          input$colour else NULL)
    })
    output$legend_scree <- renderText(
      legend_pca(input$phenos, n = pca()$n, scale = input$scale, what = "scree"))
    output$legend_contrib <- renderText({
      res <- pca()
      legend_pca(input$phenos, n = res$n, scale = input$scale,
        what = "contributions", pcx = input$pcx, pcy = input$pcy,
        varx = res$eig$variance_percent[input$pcx],
        vary = res$eig$variance_percent[input$pcy])
    })
    output$legend_eig <- renderText(
      legend_pca(input$phenos, n = pca()$n, scale = input$scale, what = "eigenvalues"))

    output$code <- renderText(code_pca(input$phenos, scale = input$scale))
  })
}
