# MDS module: thin UI over mvapp::run_mds(). Reuses the .corr_num_phenos /
# .corr_factors helpers defined in mod_correlations.R (all app R files share one
# environment in dev mode).

mod_mds_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      uiOutput(ns("pheno_ui")),
      checkboxInput(ns("scale"), "Scale traits to unit variance", TRUE),
      selectInput(ns("dist"), "Distance metric",
                  c("euclidean", "manhattan", "maximum", "canberra")),
      uiOutput(ns("colour_ui")),
      actionButton(ns("go"), "Run MDS", class = "btn-primary")
    ),
    mainPanel(
      textOutput(ns("gof")),
      checkboxInput(ns("show_legend"), "Show figure legend"),
      conditionalPanel("input.show_legend == true", ns = ns,
        textOutput(ns("legend"))),
      plotly::plotlyOutput(ns("plot"), height = "500px"),
      checkboxInput(ns("show_code"), "Show me the R code"),
      conditionalPanel("input.show_code == true", ns = ns,
        verbatimTextOutput(ns("code"))),
      tags$hr(),
      DT::DTOutput(ns("points")),
      dl_btn(ns("points_dl"))
    )
  )
}

mod_mds_server <- function(id, dataset, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    logmsg <- function(m) if (is.function(log_add)) log_add(m)

    output$pheno_ui <- renderUI({
      ds <- dataset(); req(ds); num <- .corr_num_phenos(ds)
      selectInput(ns("phenos"), "Phenotypes", num, selected = num, multiple = TRUE)
    })
    output$colour_ui <- renderUI({
      ds <- dataset(); req(ds)
      selectInput(ns("colour"), "Colour by", c("<none>", .corr_factors(ds)),
        selected = if (length(ds$roles$group)) ds$roles$group[1] else "<none>")
    })

    mds <- eventReactive(input$go, {
      ds <- dataset(); req(ds, input$phenos)
      validate(need(length(input$phenos) >= 2, "Select at least two phenotypes."))
      res <- run_mds(ds, phenotypes = input$phenos, scale = input$scale,
                     dist_method = input$dist)
      logmsg(sprintf("MDS run on %d traits (%s), %s distance; GOF = %.3f.",
                     length(input$phenos), paste(input$phenos, collapse = ", "),
                     input$dist, res$gof[1]))
      res
    })

    plot_df <- reactive({
      res <- mds(); ds <- dataset()
      d <- res$points
      # attach the colour factor by matching sample id back to the data
      if (!is.null(input$colour) && input$colour != "<none>") {
        idc <- c(ds$roles$id, ds$roles$group, ds$roles$time)
        key <- do.call(paste, c(ds$data[idc], sep = "_"))
        d[[input$colour]] <- ds$data[[input$colour]][match(d$id, make.unique(key))]
      }
      d
    })

    output$gof <- renderText({
      res <- mds()
      sprintf("Goodness-of-fit: %.3f (proportion of distance variation captured by the 2-D map).",
              res$gof[1])
    })

    output$plot <- plotly::renderPlotly({
      d <- plot_df()
      col_on <- !is.null(input$colour) && input$colour != "<none>"
      if (col_on) d[[input$colour]] <- factor(d[[input$colour]])  # discrete colour
      aes_args <- list(x = as.name("MDS1"), y = as.name("MDS2"),
                       text = as.name("id"))
      if (col_on)
        aes_args$colour <- as.name(input$colour)
      p <- ggplot2::ggplot(d, do.call(ggplot2::aes, aes_args)) +
        ggplot2::geom_point(size = 2.5, alpha = 0.85) +
        ggplot2::labs(title = "MDS configuration") +
        ggplot2::theme_minimal(base_size = 13)
      if (!is.null(input$colour) && input$colour != "<none>" &&
          length(unique(d[[input$colour]])) <= 10)
        p <- p + ggsci::scale_colour_d3()
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "colour"))
    })

    output$points <- DT::renderDT(
      DT::datatable(mds()$points, rownames = FALSE, options = list(scrollX = TRUE)))
    output$points_dl <- csv_dl(reactive(mds()$points), "MDS_coordinates.csv")

    output$legend <- renderText({
      res <- mds()
      legend_mds(input$phenos, n = res$n, scale = input$scale,
                 dist_method = input$dist, gof = res$gof[1])
    })
    output$code <- renderText(
      code_mds(input$phenos, scale = input$scale, dist_method = input$dist))
  })
}
