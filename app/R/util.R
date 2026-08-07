# Small shared UI/server helpers for the MVApp modules.

# A CSV download handler for a reactive data frame. Assign its result to an
# output slot whose id matches a downloadButton in the UI, e.g.
#   output$foo_dl <- csv_dl(reactive(mytable()), "mytable.csv")
csv_dl <- function(react, fname) {
  downloadHandler(
    filename = function() fname,
    content = function(f) utils::write.csv(react(), f, row.names = FALSE))
}

# A compact download button used under tables.
dl_btn <- function(id) downloadButton(id, "Download CSV", class = "btn-sm")

# "Show figure legend" checkbox + (hidden) text, placed above a plot. The server
# should render output `legend_<id>`; the checkbox input is `show_legend_<id>`.
# null-coalescing helper: a %||% b returns a unless it is NULL/empty, else b
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

legend_ui <- function(ns, id, label = "Show figure legend") {
  tagList(
    checkboxInput(ns(paste0("show_legend_", id)), label),
    conditionalPanel(sprintf("input.show_legend_%s == true", id), ns = ns,
      textOutput(ns(paste0("legend_", id)))))
}

# ---- categorical colour/fill palettes (ggsci journals + Brewer) -------------
cat_palette_choices <- c(
  "ggsci: D3", "ggsci: NPG", "ggsci: Lancet", "ggsci: JAMA", "ggsci: NEJM",
  "ggsci: AAAS", "ggsci: JCO", "ggsci: IGV", "ggsci: Simpsons",
  "Brewer: Set1", "Brewer: Dark2", "Brewer: Paired", "Brewer: Set2",
  "Default (ggplot)")

cat_palette_capacity <- c(
  "ggsci: D3" = 10, "ggsci: NPG" = 10, "ggsci: Lancet" = 9, "ggsci: JAMA" = 7,
  "ggsci: NEJM" = 8, "ggsci: AAAS" = 10, "ggsci: JCO" = 10, "ggsci: IGV" = 51,
  "ggsci: Simpsons" = 16, "Brewer: Set1" = 9, "Brewer: Dark2" = 8,
  "Brewer: Paired" = 12, "Brewer: Set2" = 8)

# Resolve a palette name + number of levels to a ggplot scale for the given
# aesthetic ("colour" or "fill"), or NULL (ggplot default) when it can't supply
# enough distinct colours.
cat_scale <- function(name, n, aes = c("colour", "fill")) {
  aes <- match.arg(aes)
  cap <- cat_palette_capacity[name]
  if (!is.na(cap) && n > cap) return(NULL)
  if (startsWith(name, "ggsci:")) {
    key <- tolower(trimws(sub("ggsci:", "", name)))
    fn <- sprintf("scale_%s_%s", if (aes == "fill") "fill" else "color", key)
    return(tryCatch(getExportedValue("ggsci", fn)(),
                    error = function(e) NULL))
  }
  if (startsWith(name, "Brewer:")) {
    pal <- trimws(sub("Brewer:", "", name))
    return(if (aes == "fill") ggplot2::scale_fill_brewer(palette = pal)
           else ggplot2::scale_colour_brewer(palette = pal))
  }
  NULL
}
