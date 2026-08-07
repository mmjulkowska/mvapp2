# MVApp 2.0 -- thin Shiny UI over the mvapp:: backend.
#
# The app deliberately contains NO statistics: every computation calls a
# function in the mvapp package. During development we load the package from
# ../mvapp so the app runs before the package is installed; in the Docker image
# mvapp is installed and picked up via library().

suppressWarnings(suppressMessages({
  library(shiny)
  library(DT)
  library(corrplot)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  library(ggsci)
}))

# --- load the mvapp backend (installed pkg preferred, else source ../mvapp/R) -
if (requireNamespace("mvapp", quietly = TRUE)) {
  library(mvapp)
} else {
  pkg_R <- normalizePath(file.path("..", "mvapp", "R"), mustWork = FALSE)
  if (dir.exists(pkg_R)) {
    invisible(lapply(list.files(pkg_R, "[.]R$", full.names = TRUE), source))
    message("mvapp loaded in dev mode from ", pkg_R)
  } else {
    stop("Cannot find the mvapp backend. Install it or run from app/.")
  }
}

# Shiny modules live in app/R -- one per analysis tab (start of Phase 1).
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

ui <- navbarPage(
  title = "MVApp 2.0",
  theme = NULL,
  tabPanel("Upload",       icon = icon("seedling"),        mod_upload_ui("upload")),
  tabPanel("Data exploration", icon = icon("magnifying-glass-chart"),
           mod_explore_ui("explore")),
  tabPanel("Curve fitting", icon = icon("wrench"),         mod_curvefit_ui("curvefit")),
  tabPanel("Correlations", icon = icon("compress"),        mod_correlations_ui("corr")),
  tabPanel("PCA",          icon = icon("braille"),         mod_pca_ui("pca")),
  tabPanel("MDS",          icon = icon("diagram-project"), mod_mds_ui("mds")),
  tabPanel("Clustering",   icon = icon("sitemap"),         mod_clustering_ui("clust")),
  tabPanel("Log-book",     icon = icon("book"),            mod_logbook_ui("logbook")),
  tabPanel("About",        icon = icon("circle-info"),
    fluidPage(
      h3("MVApp 2.0 -- development build"),
      p("This is the modernized, package-backed MVApp. The production app on ",
        "the KAUST server is untouched; this repo is the sandbox."),
      p("Every analysis tab is a Shiny module backed by tested functions in ",
        "the mvapp package. Add features by adding a function + a module.")
    )
  )
)

server <- function(input, output, session) {
  logbook <- new_logbook()                   # session-wide processing log
  ds <- mod_upload_server("upload", log_add = logbook$add)
  mod_explore_server("explore", dataset = ds, log_add = logbook$add)
  mod_curvefit_server("curvefit", dataset = ds, log_add = logbook$add)
  mod_correlations_server("corr", dataset = ds, log_add = logbook$add)
  mod_pca_server("pca", dataset = ds, log_add = logbook$add)
  mod_mds_server("mds", dataset = ds, log_add = logbook$add)
  mod_clustering_server("clust", dataset = ds, log_add = logbook$add)
  mod_logbook_server("logbook", logbook)
}

shinyApp(ui, server)
