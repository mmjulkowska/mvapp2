#' Time Of Experiment (hours since start)
#'
#' Derives elapsed time from timestamps *relative to the earliest observation*
#' rather than a hard-coded epoch (a fragility called out in the PhenoSight
#' notebooks). Accepts POSIXct or any string parseable by [as.POSIXct()].
#'
#' @param datetime Vector of timestamps (POSIXct or character).
#' @param origin Optional reference time; defaults to the minimum of `datetime`.
#' @param units Passed to [difftime()]; default "hours".
#' @return Numeric vector of elapsed time in `units`.
#' @export
compute_toe <- function(datetime, origin = NULL, units = "hours") {
  dt <- if (inherits(datetime, "POSIXct")) datetime else as.POSIXct(datetime, tz = "UTC")
  if (all(is.na(dt))) stop("compute_toe(): could not parse any timestamps", call. = FALSE)
  if (is.null(origin)) origin <- min(dt, na.rm = TRUE)
  as.numeric(difftime(dt, origin, units = units))
}
