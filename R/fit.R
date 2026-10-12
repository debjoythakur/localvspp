#' Fit global and Haar-localized penalized point-process models
#'
#' Runs the whole pipeline on a data frame of events: Berman-Turner
#' quadrature, covariate smoothing, the finest-level Haar-localized design
#' and the requested CBIC-tuned fits.
#'
#' @param data Data frame of events with longitude/latitude and covariate
#'   columns, e.g. [stl_crime].
#' @param methods Any of `"GLI"`, `"GSI"`, `"GALI"`, `"LLI"`, `"LSI"`.
#' @param ntile Quadrature grid size (x, y).
#' @param J Number of Haar resolution levels; the localized design uses the
#'   finest wavelet level `J - 1`.
#' @param sigma Covariate smoothing bandwidth.
#' @param quad_weights Passed to [bt_quadrature()] as `weights`.
#' @param covariates Optional data frame of covariates (one row per event);
#'   defaults to [crime_covariates()]`(data)`.
#' @param tol Threshold for treating a coefficient as non-zero.
#' @param scad_args Extra arguments for [fit_scad_cbic()], e.g. `list(mc.cores = 1)`.
#' @param verbose Print progress messages?
#' @return An object of class `localvs`.
#' @examples
#' \dontrun{
#' set.seed(2026)                       # adaptive-LASSO weights use random CV folds
#' fit <- fit_localvs(stl_crime)
#' fit
#' selected_covariates(fit)$LLI
#' ev <- evaluate_localvs(fit)
#' plot_performance(fit, ev)
#' }
#' @export
fit_localvs <- function(data, methods = names(method_labels), ntile = c(32, 32), J = 3,
                        sigma = 0.01, quad_weights = c("spatstat", "cell"), covariates = NULL,
                        tol = 1e-8, scad_args = list(), verbose = TRUE) {
  methods <- match.arg(methods, names(method_labels), several.ok = TRUE)
  data <- standardize_coords(data)
  X    <- if (is.null(covariates)) crime_covariates(data) else as.data.frame(covariates)
  if (nrow(X) != nrow(data)) stop("covariates must have one row per event.", call. = FALSE)

  say  <- function(...) if (verbose) message(...)
  quad <- bt_quadrature(as.matrix(data[c("longitude", "latitude")]), ntile = ntile,
                        weights = match.arg(quad_weights))
  say(sprintf("Quadrature: %d points (%d events, %d dummy)", nrow(quad$U), sum(quad$y), sum(quad$y == 0)))
  haar <- finest_level(haar_basis_2d(quad$U, J = J))
  say("Smoothing ", ncol(X), " covariates to the quadrature points")
  Xq   <- smooth_covariates(X, quad, sigma = sigma, dimyx = ntile)
  Z    <- localized_design(Xq, haar$B, haar$meta)
  keep_X <- ok_columns(Xq); keep_Z <- ok_columns(Z)
  say(sprintf("Design: %d global and %d localized candidate coefficients", sum(keep_X), sum(keep_Z)))

  y <- quad$y; off <- quad$offset
  Xg <- Xq[, keep_X, drop = FALSE]; Zc <- Z[, keep_Z, drop = FALSE]
  scad <- function(D) do.call(fit_scad_cbic, c(list(D, y, off, coef_tol = tol), scad_args))
  fitters <- list(GLI  = function() fit_lasso_cbic(Xg, y, off, tol = tol),
                  GSI  = function() scad(Xg),
                  GALI = function() fit_adaptive_lasso_cbic(Xg, y, off, tol = tol),
                  LLI  = function() fit_lasso_cbic(Zc, y, off, tol = tol),
                  LSI  = function() scad(Zc))
  fits <- list()
  for (m in methods) {
    say("Fitting ", m, " (", method_labels[[m]], ")")
    fits[[m]] <- fitters[[m]]()
  }

  coefs <- lapply(methods, function(m) {
    local <- m %in% local_methods
    keep  <- if (local) keep_Z else keep_X
    b <- stats::setNames(numeric(length(keep)), if (local) colnames(Z) else colnames(Xq))
    b[keep] <- fits[[m]]$beta
    b
  })
  names(coefs) <- methods

  selection <- lapply(methods, function(m) {
    b <- coefs[[m]]
    if (m %in% local_methods) return(map_selection(b, haar$meta, colnames(Xq), tol))
    nz <- which(abs(b) > tol)
    data.frame(covariate = names(b)[nz], coefficient = unname(b[nz]), stringsAsFactors = FALSE)
  })
  names(selection) <- methods

  structure(list(call = match.call(), data = data, covariates = X, quad = quad, B = haar$B,
                 meta = haar$meta, Xq = Xq, Z = Z, keep_X = keep_X, keep_Z = keep_Z, fits = fits,
                 coefficients = coefs, selection = selection, tol = tol),
            class = "localvs")
}

## Non-zero (covariate, atom) pairs of a localized coefficient vector
map_selection <- function(beta, meta, cov_names, tol = 1e-8) {
  G  <- matrix(beta, nrow = nrow(meta), ncol = length(cov_names))
  nz <- which(abs(G) > tol, arr.ind = TRUE)
  if (!nrow(nz)) {
    return(cbind(data.frame(covariate = character(0), coefficient = numeric(0)), meta[0, , drop = FALSE]))
  }
  out <- cbind(data.frame(covariate = cov_names[nz[, 2]], coefficient = G[nz], stringsAsFactors = FALSE),
               meta[nz[, 1], , drop = FALSE])
  out <- out[order(out$covariate, out$j, out$type), ]
  rownames(out) <- NULL
  out
}

#' Selected covariates per method
#'
#' @param object A `localvs` fit.
#' @return A named list of character vectors. For local methods a covariate
#'   counts as selected if it is active on at least one Haar atom.
#' @export
selected_covariates <- function(object) {
  lapply(object$selection, function(s) sort(unique(s$covariate)))
}

#' Covariates selected at a location (local methods)
#'
#' @param object A `localvs` fit.
#' @param lon,lat Location (same coordinates as the data).
#' @param method `"LLI"` or `"LSI"`.
#' @return Character vector of covariates whose active Haar atoms cover the location.
#' @export
selected_at <- function(object, lon, lat, method = "LLI") {
  method <- match.arg(method, local_methods)
  s <- object$selection[[method]]
  if (is.null(s)) stop(method, " was not fitted.", call. = FALSE)
  hit <- s$x_lo <= lon & s$x_hi >= lon & s$y_lo <= lat & s$y_hi >= lat
  sort(unique(s$covariate[hit]))
}

#' @export
print.localvs <- function(x, ...) {
  cat("Local variable selection for a spatial point process\n")
  print(x$quad)
  cat("Candidate covariates:", ncol(x$Xq), "| localized coefficients:", sum(x$keep_Z), "\n\n")
  sel <- selected_covariates(x)
  tab <- data.frame(
    Method     = names(x$fits),
    Label      = unname(method_labels[names(x$fits)]),
    lambda     = vapply(x$fits, `[[`, numeric(1), "lambda"),
    CBIC       = vapply(x$fits, `[[`, numeric(1), "cbic"),
    Nonzero    = vapply(x$coefficients, function(b) sum(abs(b) > x$tol), numeric(1)),
    Covariates = lengths(sel)[names(x$fits)],
    row.names  = NULL)
  print(tab, row.names = FALSE, digits = 4)
  invisible(x)
}
