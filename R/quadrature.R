#' Berman-Turner quadrature scheme
#'
#' Builds a grid quadrature scheme on the (slightly padded) bounding
#' rectangle of the events. Data points get `y = 1`, dummy points `y = 0`,
#' and `log(w)` is the Poisson-regression offset.
#'
#' @param coords Two-column matrix of event coordinates (x, y).
#' @param ntile Number of grid tiles in x and y.
#' @param pad Fraction of each range added on both sides of the window.
#' @param weights `"spatstat"` uses the quadrature weights computed by
#'   [spatstat.geom::quadscheme()] (the standard Berman-Turner weights).
#'   `"cell"` gives every quadrature point the same weight, the grid-cell area.
#'   The original analysis script called `weights(Q)`, which returns `NULL`
#'   for a quadrature object and therefore fell back to `"cell"`; use
#'   `weights = "cell"` to reproduce it.
#' @return An object of class `bt_quad`: a list with the quadrature points `U`,
#'   indicator `y`, weights `w`, `offset`, `window`, the event pattern `ppp`,
#'   `ntile` and `cell_area`.
#' @export
bt_quadrature <- function(coords, ntile = c(32, 32), pad = 0.01,
                          weights = c("spatstat", "cell")) {
  weights <- match.arg(weights)
  coords  <- as.matrix(coords)
  padr <- function(r) r + c(-1, 1) * pad * diff(r)
  W  <- spatstat.geom::owin(padr(range(coords[, 1])), padr(range(coords[, 2])))
  X  <- spatstat.geom::ppp(coords[, 1], coords[, 2], window = W, checkdup = FALSE)
  Q  <- spatstat.geom::quadscheme(X, method = "grid", ntile = ntile)
  P  <- spatstat.geom::union.quad(Q)
  U  <- data.frame(x = P$x, y = P$y)
  y  <- as.integer(spatstat.geom::is.data(Q))
  cell_area <- diff(W$xrange) * diff(W$yrange) / prod(ntile)
  w  <- if (weights == "spatstat") as.numeric(spatstat.geom::w.quad(Q)) else rep(cell_area, nrow(U))
  w[!is.finite(w) | w <= 0] <- cell_area
  structure(list(U = U, y = y, w = w, offset = log(pmax(w, 1e-12)), window = W, ppp = X,
                 ntile = ntile, cell_area = cell_area, weights = weights),
            class = "bt_quad")
}

#' @export
print.bt_quad <- function(x, ...) {
  cat("Berman-Turner quadrature:", nrow(x$U), "points |", sum(x$y), "events |",
      sum(x$y == 0), "dummy | weights:", x$weights, "\n")
  invisible(x)
}

## Grid-cell index of each quadrature point (used for count comparisons)
grid_cell_id <- function(quad) {
  cell <- function(v, r, n) findInterval(v, seq(r[1], r[2], length.out = n + 1), all.inside = TRUE)
  W <- quad$window
  (cell(quad$U$y, W$yrange, quad$ntile[2]) - 1L) * quad$ntile[1] + cell(quad$U$x, W$xrange, quad$ntile[1])
}
