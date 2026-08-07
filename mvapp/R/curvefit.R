# Curve fitting -- ports the MVApp curve-fitting tab. For each "sample" (a unique
# combination of grouping columns) a model is fitted of a response vs a
# continuous time/dose variable, and the intercept, primary rate coefficient
# (DELTA), and goodness-of-fit (R^2) are extracted. Non-linear shapes are fitted
# as linear models on transformed variables; R^2 is always computed on the
# original scale (correlation between observed and fitted, squared) so it is
# comparable across model types.

# supported models (label = internal key)
.curve_models <- c("Linear" = "linear", "Quadratic" = "quadratic",
  "Exponential" = "exponential", "Square-root" = "square_root",
  "Logarithmic" = "logarithmic")

# fit one sample; returns intercept, delta (primary coefficient), r2
.fit_one <- function(x, y, model) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  na <- list(intercept = NA_real_, delta = NA_real_, r2 = NA_real_)
  if (length(x) < 3 || length(unique(x)) < 2) return(na)
  fit <- tryCatch(switch(model,
    linear      = { m <- stats::lm(y ~ x); list(p = stats::fitted(m), co = stats::coef(m)) },
    quadratic   = { m <- stats::lm(y ~ x + I(x^2)); list(p = stats::fitted(m), co = stats::coef(m)) },
    exponential = { if (any(y <= 0)) return(na)
                    m <- stats::lm(log(y) ~ x); list(p = exp(stats::fitted(m)), co = stats::coef(m)) },
    square_root = { if (any(x < 0)) return(na)
                    m <- stats::lm(y ~ sqrt(x)); list(p = stats::fitted(m), co = stats::coef(m)) },
    logarithmic = { if (any(x <= 0)) return(na)
                    m <- stats::lm(y ~ log(x)); list(p = stats::fitted(m), co = stats::coef(m)) },
    return(na)), error = function(e) NULL)
  if (is.null(fit)) return(na)
  r2 <- suppressWarnings(stats::cor(y, fit$p)^2)
  list(intercept = unname(fit$co[1]), delta = unname(fit$co[2]),
       r2 = if (is.na(r2)) NA_real_ else r2)
}

#' Fit a curve per sample
#'
#' @param data A data frame or [mvapp_dataset].
#' @param time Continuous predictor column (e.g. day, dose).
#' @param response Numeric response column (the phenotype to model).
#' @param samples Character vector of columns whose unique combinations define a
#'   single curve to fit -- typically the individual plant / sample ID (so one
#'   curve is fitted per plant, not per whole genotype x treatment group).
#' @param groups Optional character vector of factor columns (e.g. genotype,
#'   treatment) carried along with each sample as metadata, so the extracted
#'   parameters can later be compared across those factors. Assumed constant
#'   within a sample; the first value per sample is kept.
#' @param model One of "linear", "quadratic", "exponential", "square_root",
#'   "logarithmic".
#' @return An object of class `mvapp_curvefit`: a `table` with one row per
#'   sample (the sample-ID and group columns, n points, INTERCEPT, DELTA, r2),
#'   plus the model and column roles.
#' @export
fit_curves <- function(data, time, response, samples, groups = character(0),
                       model = "linear") {
  if (inherits(data, "mvapp_dataset")) data <- data$data
  stopifnot(is.data.frame(data))
  for (v in c(time, response, samples, groups)) {
    if (!v %in% names(data)) stop("fit_curves(): column not found: ", v, call. = FALSE)
  }
  if (length(samples) < 1) stop("fit_curves(): need at least one sample-ID column", call. = FALSE)
  # a grouping factor may also be part of the plant ID (e.g. accession is both);
  # keep each underlying column once, but still remember the full `groups` so the
  # extracted parameters can be compared across them.
  keep <- unique(c(samples, groups))
  d <- data[stats::complete.cases(data[c(time, response, samples)]),
            c(keep, time, response), drop = FALSE]
  d$.sample <- do.call(paste, c(d[samples], sep = " | "))
  parts <- split(d, d$.sample)
  rows <- lapply(parts, function(s) {
    f <- .fit_one(s[[time]], s[[response]], model)
    cbind(data.frame(sample = s$.sample[1], stringsAsFactors = FALSE),
          s[1, keep, drop = FALSE],
          data.frame(n = nrow(s), INTERCEPT = f$intercept, DELTA = f$delta,
                     r2 = f$r2))
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  tab <- tab[order(tab$r2), , drop = FALSE]     # worst fits first
  structure(list(table = tab, model = model, time = time,
                 response = response, samples = samples, groups = groups),
            class = "mvapp_curvefit")
}

#' Estimate which model fits best per group
#'
#' Pools ALL plants within each group (unique combination of `groups`, e.g.
#' genotype x treatment) and fits every model to the whole group at once, then
#' returns each model's R^2 plus the best-fitting model. This answers "which
#' shape best describes this group of plants overall", rather than fitting a
#' separate model to every individual plant. "Best" is chosen by *adjusted* R^2,
#' which penalises extra parameters, so a simpler model wins unless a more
#' complex one fits genuinely better (avoiding the bias where a quadratic always
#' beats a straight line).
#'
#' @inheritParams fit_curves
#' @param models Models to compare.
#' @return A data frame: sample (group) + one R^2 column per model + `best_model`.
#' @export
estimate_models <- function(data, time, response, groups,
                            models = unname(.curve_models)) {
  npred <- c(linear = 1, quadratic = 2, exponential = 1, square_root = 1,
             logarithmic = 1)
  base <- NULL; nvec <- NULL
  for (m in models) {
    fc <- fit_curves(data, time, response, samples = groups, model = m)$table
    if (is.null(nvec)) nvec <- fc[c("sample", "n")]
    col <- fc[c("sample", "r2")]; names(col)[2] <- m
    base <- if (is.null(base)) col else merge(base, col, by = "sample", sort = FALSE)
  }
  base <- merge(base, nvec, by = "sample", sort = FALSE)
  adj <- vapply(models, function(m) {
    p <- npred[[m]]
    a <- 1 - (1 - base[[m]]) * (base$n - 1) / (base$n - p - 1)
    ifelse(base$n > p + 1, a, base[[m]])
  }, numeric(nrow(base)))
  if (is.null(dim(adj))) adj <- matrix(adj, nrow = nrow(base))
  adj[is.na(adj)] <- -Inf
  base$best_model <- models[max.col(adj, ties.method = "first")]
  base
}

#' Fitted curve for one sample (for plotting)
#'
#' @param x,y Numeric predictor and response for a single sample.
#' @param model Model key (see [fit_curves()]).
#' @param n Number of points along the fitted curve.
#' @return A data frame (x, y) of the fitted curve, or NULL if it can't be fit.
#' @export
curve_fit_line <- function(x, y, model, n = 100) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2) return(NULL)
  xg <- seq(min(x), max(x), length.out = n)
  pred <- tryCatch(switch(model,
    linear      = stats::predict(stats::lm(y ~ x), newdata = data.frame(x = xg)),
    quadratic   = stats::predict(stats::lm(y ~ x + I(x^2)), newdata = data.frame(x = xg)),
    exponential = { if (any(y <= 0)) return(NULL)
                    exp(stats::predict(stats::lm(log(y) ~ x), newdata = data.frame(x = xg))) },
    square_root = { if (any(x < 0)) return(NULL)
                    stats::predict(stats::lm(y ~ sqrt(x)), newdata = data.frame(x = xg)) },
    logarithmic = { if (any(x <= 0)) return(NULL)
                    stats::predict(stats::lm(y ~ log(x)), newdata = data.frame(x = xg)) },
    return(NULL)), error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  data.frame(x = xg, y = as.numeric(pred))
}

#' @export
print.mvapp_curvefit <- function(x, ...) {
  cat("<mvapp_curvefit>", x$model, " model of", x$response, "vs", x$time,
      "\n  samples:", nrow(x$table),
      " median R2:", round(stats::median(x$table$r2, na.rm = TRUE), 3), "\n")
  invisible(x)
}
