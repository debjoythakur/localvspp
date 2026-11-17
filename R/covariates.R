#' Smooth event-level covariates to the quadrature points
#'
#' Each covariate is attached to the events as marks, smoothed with
#' [spatstat.explore::Smooth()] (Nadaraya-Watson) on a pixel grid and read
#' off at the nearest pixel of every quadrature point.
#'
#' @param X Data frame of covariates, one row per event.
#' @param quad A [bt_quadrature()] object built from the same events.
#' @param sigma Smoothing bandwidth (in coordinate units).
#' @param dimyx Pixel grid dimensions; defaults to the quadrature tiles.
#' @param standardize Center and scale each column?
#' @return A numeric matrix with one row per quadrature point.
#' @export
smooth_covariates <- function(X, quad, sigma = 0.01, dimyx = quad$ntile, standardize = TRUE) {
  pp <- quad$ppp
  at_points <- function(v) {
    spatstat.geom::marks(pp) <- v
    im <- spatstat.explore::Smooth(pp, sigma = sigma, at = "pixels", normalise = TRUE, dimyx = dimyx)
    np <- spatstat.geom::nearest.pixel(quad$U$x, quad$U$y, im)
    im$v[cbind(pmin(pmax(np$row, 1L), nrow(im$v)), pmin(pmax(np$col, 1L), ncol(im$v)))]
  }
  Xq <- vapply(as.data.frame(X), at_points, numeric(nrow(quad$U)))
  Xq <- matrix(Xq, nrow = nrow(quad$U), dimnames = list(NULL, names(X)))
  Xq[!is.finite(Xq)] <- 0
  if (standardize) {
    Xq <- scale(Xq)
    Xq[!is.finite(Xq)] <- 0
  }
  Xq
}

#' Haar-localized design matrix
#'
#' Multiplies every covariate column by every Haar atom, so each covariate
#' gets one coefficient per atom (a spatially varying coefficient).
#'
#' @param X Covariate matrix at the quadrature points.
#' @param B Haar basis matrix at the same points.
#' @param meta Atom description (rows match the columns of `B`).
#' @return A matrix with `ncol(X) * ncol(B)` columns named
#'   `<covariate>__<type>_j<j>_kx<kx>_ky<ky>`; covariate-major order.
#' @export
localized_design <- function(X, B, meta) {
  X <- as.matrix(X); B <- as.matrix(B)
  if (nrow(X) != nrow(B) || nrow(meta) != ncol(B)) stop("Incompatible X, B and meta dimensions.", call. = FALSE)
  basis_names <- paste0(meta$type, "_j", meta$j, "_kx", meta$kx, "_ky", meta$ky)
  Z <- do.call(cbind, lapply(seq_len(ncol(X)), function(p) B * X[, p]))
  colnames(Z) <- paste(rep(colnames(X), each = ncol(B)), basis_names, sep = "__")
  Z
}

## Columns that are finite and non-constant
ok_columns <- function(X) apply(X, 2, function(z) all(is.finite(z)) && isTRUE(stats::sd(z) > 0))
