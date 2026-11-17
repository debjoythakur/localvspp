haar_phi <- function(j, k, x) 2^(j / 2) * (x >= k / 2^j & x < (k + 1) / 2^j)

haar_psi <- function(j, k, x) {
  m <- (k + 0.5) / 2^j
  2^(j / 2) * ((x >= k / 2^j & x < m) - (x >= m & x < (k + 1) / 2^j))
}

## Rescale each column to [0, 1). Points on the upper edge are pulled just
## inside so that they fall in the last block instead of in no block at all.
to_unit_square <- function(s) {
  s <- as.matrix(s)
  u <- apply(s, 2, function(v) (v - min(v)) / max(diff(range(v)), 1e-12))
  pmin(u, 1 - 1e-12)
}

## Map unit-square supports back to the original coordinates
add_supports <- function(meta, xy) {
  xy <- as.matrix(xy); xr <- range(xy[, 1]); yr <- range(xy[, 2])
  meta$x_lo <- xr[1] + meta$x_lo01 * diff(xr); meta$x_hi <- xr[1] + meta$x_hi01 * diff(xr)
  meta$y_lo <- yr[1] + meta$y_lo01 * diff(yr); meta$y_hi <- yr[1] + meta$y_hi01 * diff(yr)
  meta
}

eval_atoms <- function(s01, meta) {
  B <- matrix(0, nrow(s01), nrow(meta))
  for (r in seq_len(nrow(meta))) {
    fx <- if (meta$type[r] %in% c("psi_x", "psi_xy")) haar_psi else haar_phi
    fy <- if (meta$type[r] %in% c("psi_y", "psi_xy")) haar_psi else haar_phi
    B[, r] <- fx(meta$j[r], meta$kx[r], s01[, 1]) * fy(meta$j[r], meta$ky[r], s01[, 2])
  }
  B
}

#' Local Haar basis at the finest level (the basis used in the paper)
#'
#' Evaluates the two-dimensional Haar atoms at level \eqn{j = J - 1} after
#' rescaling the locations to the unit square. Each of the \eqn{4^{J-1}}
#' blocks of side \eqn{2^{-(J-1)}} carries one scaling atom (`phi`, the level
#' of the coefficient on the block) and three detail atoms (`psi_x`, `psi_y`,
#' `psi_xy`: horizontal, vertical and diagonal contrasts), giving
#' \eqn{R = 4 \cdot 4^{J-1} = 4^J} atoms in total, each supported on a single
#' block. Together they span every surface that is constant on the
#' \eqn{2^J \times 2^J} sub-blocks (Eq. (5) of the paper).
#'
#' @param s Two-column matrix or data frame of locations.
#' @param J Number of resolution levels; atoms are taken at level `J - 1`.
#' @param xy Locations used to map unit-square supports back to the original
#'   coordinates (`x_lo`, `x_hi`, `y_lo`, `y_hi` in `meta`). Defaults to `s`.
#' @return A list with the basis matrix `B` (one column per atom) and a data
#'   frame `meta` describing each atom (`j`, `type`, `kx`, `ky`, `atom_index`
#'   and its support).
#' @examples
#' h <- local_haar_basis(cbind(runif(50), runif(50)), J = 3)
#' dim(h$B)            # 50 x 64
#' table(h$meta$type)  # 16 of each type
#' @export
local_haar_basis <- function(s, J = 3, xy = s) {
  if (J < 1) stop("J must be at least 1.", call. = FALSE)
  j <- J - 1
  g <- expand.grid(type = c("phi", "psi_x", "psi_y", "psi_xy"), ky = 0:(2^j - 1), kx = 0:(2^j - 1),
                   stringsAsFactors = FALSE)
  meta <- data.frame(j = j, type = g$type, kx = g$kx, ky = g$ky, stringsAsFactors = FALSE)
  meta$atom_index <- seq_len(nrow(meta))
  meta$x_lo01 <- meta$kx / 2^j; meta$x_hi01 <- (meta$kx + 1) / 2^j
  meta$y_lo01 <- meta$ky / 2^j; meta$y_hi01 <- (meta$ky + 1) / 2^j
  B <- eval_atoms(to_unit_square(s), meta)
  if (!is.null(xy)) meta <- add_supports(meta, xy)
  list(B = B, meta = meta)
}

#' Two-dimensional multi-resolution Haar wavelet basis
#'
#' Evaluates the full tensor-product Haar basis after rescaling the locations
#' to the unit square: scaling functions (`phi`) at level `j0` and wavelets
#' (`psi_x`, `psi_y`, `psi_xy`) at levels `j0, ..., J - 1`. The localized
#' fits in [fit_localvs()] use [local_haar_basis()] instead.
#'
#' @inheritParams local_haar_basis
#' @param j0 Coarsest level.
#' @return A list with the basis matrix `B` and the atom description `meta`.
#' @examples
#' h <- haar_basis_2d(cbind(runif(50), runif(50)), J = 2)
#' dim(h$B)
#' head(h$meta)
#' @export
haar_basis_2d <- function(s, J = 3, j0 = 0, xy = s) {
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
  B <- eval_atoms(to_unit_square(s), meta)
  if (!is.null(xy)) meta <- add_supports(meta, xy)
  list(B = B, meta = meta)
}

#' Finest-level detail atoms (deprecated)
#'
#' Keeps only the detail atoms of the finest level of [haar_basis_2d()],
#' without the scaling atoms. This was the localized basis in localvspp 0.1.0;
#' the paper's basis also includes the finest-level scaling atoms and is
#' built by [local_haar_basis()].
#'
#' @param haar Output of [haar_basis_2d()].
#' @return The same structure restricted to the finest wavelet atoms, with
#'   `atom_index` renumbered.
#' @export
finest_level <- function(haar) {
  .Deprecated("local_haar_basis", package = "localvspp")
  wav  <- haar$meta$type != "phi"
  keep <- wav & haar$meta$j == max(haar$meta$j[wav])
  meta <- haar$meta[keep, , drop = FALSE]
  meta$atom_index <- seq_len(nrow(meta)); rownames(meta) <- NULL
  list(B = haar$B[, keep, drop = FALSE], meta = meta)
}
