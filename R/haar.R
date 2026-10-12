haar_phi <- function(j, k, x) 2^(j / 2) * (x >= k / 2^j & x < (k + 1) / 2^j)

haar_psi <- function(j, k, x) {
  m <- (k + 0.5) / 2^j
  2^(j / 2) * ((x >= k / 2^j & x < m) - (x >= m & x < (k + 1) / 2^j))
}

to_unit_square <- function(s) {
  s <- as.matrix(s)
  apply(s, 2, function(v) (v - min(v)) / max(diff(range(v)), 1e-12))
}

#' Two-dimensional Haar wavelet basis
#'
#' Evaluates the tensor-product Haar basis at the given locations after
#' rescaling them to the unit square: scaling functions (`phi`) at level `j0`
#' and wavelets (`psi_x`, `psi_y`, `psi_xy`) at levels `j0, ..., J - 1`.
#'
#' @param s Two-column matrix or data frame of locations.
#' @param J Number of resolution levels (finest wavelet level is `J - 1`).
#' @param j0 Coarsest level.
#' @param xy Locations used to map unit-square supports back to the original
#'   coordinates (`x_lo`, `x_hi`, `y_lo`, `y_hi` in `meta`). Defaults to `s`.
#' @return A list with the basis matrix `B` (one column per atom) and a data
#'   frame `meta` describing each atom (`j`, `type`, `kx`, `ky`, its support).
#' @examples
#' h <- haar_basis_2d(cbind(runif(50), runif(50)), J = 2)
#' dim(h$B)
#' head(h$meta)
#' @export
haar_basis_2d <- function(s, J = 3, j0 = 0, xy = s) {
  s01 <- to_unit_square(s)
  grid_j <- function(j, types) {
    g <- expand.grid(type = types, ky = 0:(2^j - 1), kx = 0:(2^j - 1), stringsAsFactors = FALSE)
    data.frame(j = j, type = g$type, kx = g$kx, ky = g$ky, stringsAsFactors = FALSE)
  }
  levels <- list(grid_j(j0, "phi"))
  if (J > j0) levels <- c(levels, lapply(j0:(J - 1), grid_j, types = c("psi_x", "psi_y", "psi_xy")))
  meta <- do.call(rbind, levels)
  meta$atom_index <- seq_len(nrow(meta))
  meta$x_lo01 <- meta$kx / 2^meta$j; meta$x_hi01 <- (meta$kx + 1) / 2^meta$j
  meta$y_lo01 <- meta$ky / 2^meta$j; meta$y_hi01 <- (meta$ky + 1) / 2^meta$j

  B <- matrix(0, nrow(s01), nrow(meta))
  for (r in seq_len(nrow(meta))) {
    fx <- if (meta$type[r] %in% c("psi_x", "psi_xy")) haar_psi else haar_phi
    fy <- if (meta$type[r] %in% c("psi_y", "psi_xy")) haar_psi else haar_phi
    B[, r] <- fx(meta$j[r], meta$kx[r], s01[, 1]) * fy(meta$j[r], meta$ky[r], s01[, 2])
  }
  if (!is.null(xy)) {
    xy <- as.matrix(xy); xr <- range(xy[, 1]); yr <- range(xy[, 2])
    meta$x_lo <- xr[1] + meta$x_lo01 * diff(xr); meta$x_hi <- xr[1] + meta$x_hi01 * diff(xr)
    meta$y_lo <- yr[1] + meta$y_lo01 * diff(yr); meta$y_hi <- yr[1] + meta$y_hi01 * diff(yr)
  }
  list(B = B, meta = meta)
}

#' Keep only the finest wavelet level
#'
#' @param haar Output of [haar_basis_2d()].
#' @return The same structure restricted to the finest wavelet atoms, with
#'   `atom_index` renumbered.
#' @export
finest_level <- function(haar) {
  wav  <- haar$meta$type != "phi"
  keep <- wav & haar$meta$j == max(haar$meta$j[wav])
  meta <- haar$meta[keep, , drop = FALSE]
  meta$atom_index <- seq_len(nrow(meta)); rownames(meta) <- NULL
  list(B = haar$B[, keep, drop = FALSE], meta = meta)
}
