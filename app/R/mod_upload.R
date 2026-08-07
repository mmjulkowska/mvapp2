# Upload module: read a long-format CSV into an mvapp_dataset and preview it.

mod_upload_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("file"), "Long-format CSV", accept = ".csv"),
      helpText("No file? The bundled example loads automatically."),
      uiOutput(ns("role_ui"))
    ),
    mainPanel(DT::DTOutput(ns("preview")))
  )
}

mod_upload_server <- function(id, log_add = NULL) {
  moduleServer(id, function(input, output, session) {
    logmsg <- function(m) if (is.function(log_add)) log_add(m)

    raw <- reactive({
      path <- if (is.null(input$file)) {
        system.file("extdata", "example_long.csv", package = "mvapp")
      } else input$file$datapath
      if (identical(path, "")) {           # dev mode: pkg not installed
        path <- normalizePath(file.path("..", "mvapp", "inst", "extdata",
                                        "example_long.csv"))
      }
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    })

    observeEvent(input$file, {
      logmsg(sprintf("Loaded file '%s' (%d rows, %d columns).",
                     input$file$name, nrow(raw()), ncol(raw())))
    })

    output$role_ui <- renderUI({
      ns <- session$ns
      cols <- names(raw())
      tagList(
        selectInput(ns("id"), "ID column(s) — identify a single plant / sample",
                    cols, selected = cols[1], multiple = TRUE),
        uiOutput(ns("group_ui")),          # excludes whatever is chosen as ID
        selectInput(ns("time"),  "Time column",     c("<none>", cols),
                    selected = if ("DAY" %in% cols) "DAY" else "<none>"),
        helpText("Both ID and Grouping accept several columns. Pick every column",
                 "that together pins down one plant as the ID, and every factor",
                 "you want to compare across as the Grouping. Columns used for",
                 "the ID are removed from the Grouping choices.")
      )
    })

    # Grouping choices exclude the ID columns; re-renders when the ID changes but
    # keeps whatever grouping picks are still valid.
    output$group_ui <- renderUI({
      ns <- session$ns
      avail <- setdiff(names(raw()), input$id)
      guess <- intersect(c("ACCESSION", "GENOTYPE", "LINE", "TREATMENT",
                           "CONDITION", "GROUP"), avail)
      cur <- intersect(isolate(input$group), avail)
      sel <- if (length(cur)) cur else if (length(guess)) guess else avail[1]
      selectInput(ns("group"),
                  "Grouping column(s) — factors to compare (genotype, treatment, …)",
                  avail, selected = sel, multiple = TRUE)
    })

    dataset <- reactive({
      req(input$id, input$group)
      time <- if (identical(input$time, "<none>")) NULL else input$time
      phen <- setdiff(names(raw()), c(input$id, input$group, time))
      mvapp_dataset(raw(), id = input$id, group = input$group,
                    time = time, phenotypes = phen)
    })

    observeEvent(list(input$id, input$group, input$time), {
      req(input$id, input$group)
      time <- if (identical(input$time, "<none>")) "<none>" else input$time
      phen <- setdiff(names(raw()), c(input$id, input$group,
                                      if (time == "<none>") NULL else time))
      logmsg(sprintf("Dataset roles set: id=%s, group=%s, time=%s; %d traits (%s).",
                     paste(input$id, collapse = "+"),
                     paste(input$group, collapse = "+"), time, length(phen),
                     paste(phen, collapse = ", ")))
    }, ignoreInit = TRUE)

    output$preview <- DT::renderDT(DT::datatable(raw(), options = list(pageLength = 8)))
    dataset
  })
}
