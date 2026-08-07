# Session log-book: a running, timestamped record of what the user did, to
# complement the reproducible R script. Modules receive `log_add()` and call it
# on meaningful actions (data load, PCA run, downloads, and -- later -- outlier
# curation). Downloadable as a plain-text provenance record.

# Create a log-book handle: a reactiveVal of lines plus an add() function.
new_logbook <- function() {
  rv <- shiny::reactiveVal(character(0))
  list(
    entries = rv,
    add = function(msg) {
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      shiny::isolate(rv(c(rv(), paste0("[", ts, "] ", msg))))
    }
  )
}

mod_logbook_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    h4("Processing log-book"),
    p("A running record of what you did in this session — data loaded, ",
      "analyses run, outputs downloaded (and, later, outlier curation). ",
      "Download it alongside your figures so the processing history travels ",
      "with the results, complementing the R script."),
    downloadButton(ns("dl"), "Download log-book (.txt)"),
    tags$hr(),
    verbatimTextOutput(ns("log"))
  )
}

mod_logbook_server <- function(id, logbook) {
  moduleServer(id, function(input, output, session) {
    output$log <- renderText({
      e <- logbook$entries()
      if (!length(e)) "No actions logged yet." else paste(e, collapse = "\n")
    })
    output$dl <- downloadHandler(
      filename = function() "MVApp_logbook.txt",
      content = function(f) {
        writeLines(c("MVApp 2.0 processing log-book", "", logbook$entries()), f)
      })
  })
}
