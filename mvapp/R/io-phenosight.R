# ---------------------------------------------------------------------------
# PhenoSight raw-export readers  ---  Roadmap Phase 2
#
# These are intentionally stubs: they fix the *interface* (arguments + return
# type = mvapp_dataset) so the import wizard and downstream code can be built
# against them, while the raw-parsing bodies are filled in during Phase 2.
# ---------------------------------------------------------------------------

#' Read a raw PhenoSight PS2 (chlorophyll-fluorescence) export  [Phase 2a]
#'
#' TODO (Phase 2a):
#'  * Apply the Fv/Fm split-row fix: PhenoVation writes Fv/Fm on one set of rows
#'    and the light-adapted traits on others (tagged by `nTmPam` or a negative
#'    sentinel). Split by sign / pulse count and re-merge on row-identity keys.
#'  * Parse `tray.ID` from the image filename (3rd `_`-delimited token).
#'  * Derive TOE from Date/Time via [compute_toe()].
#'  * Join a coding file per the chosen aggregation regime.
#'
#' @param path Path to a `.INF` / `.TXT` / `*_Processed.csv` export.
#' @param coding_file Optional tray coding/decode file.
#' @param regime Aggregation regime: "A" per-position, "B" per-tray, "C" per-pot.
#' @return An `mvapp_dataset` with `source = "phenosight_ps2"`.
#' @export
read_ps2 <- function(path, coding_file = NULL, regime = c("A", "B", "C")) {
  stop("read_ps2(): not yet implemented -- see Roadmap Phase 2a", call. = FALSE)
}

#' Read a raw PhenoSight RGB export (digital biomass / growth)  [Phase 2a]
#'
#' TODO (Phase 2a): digital biomass = sum of the 7 area measurements
#' (`area` + `Side01_00..05_area`); sum objects back per image; recover keys.
#'
#' @inheritParams read_ps2
#' @return An `mvapp_dataset` with `source = "phenosight_rgb"`.
#' @export
read_rgb <- function(path, coding_file = NULL, regime = c("A", "B", "C")) {
  stop("read_rgb(): not yet implemented -- see Roadmap Phase 2a", call. = FALSE)
}

#' Read a PhenoSight watering / evapotranspiration export  [Phase 2a]
#'
#' TODO (Phase 2a): support all three shapes -- Arduino Before/After, per-tray
#' H2OData, and the semicolon-delimited platform log (Water = Mass2 - Mass1).
#'
#' @inheritParams read_ps2
#' @return An `mvapp_dataset` with `source = "phenosight_watering"`.
#' @export
read_watering <- function(path, coding_file = NULL, regime = c("A", "B", "C")) {
  stop("read_watering(): not yet implemented -- see Roadmap Phase 2a", call. = FALSE)
}
