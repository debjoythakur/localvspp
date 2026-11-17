shape_theme <- function() {
  ggplot2::theme_bw(base_size = 13) + ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
    legend.position = "bottom", legend.title = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank())
}

#' Compare methods: log(MSE), Hellinger distance and selection proportion
#'
#' @param object A `localvs` fit.
#' @param evaluation Output of [evaluate_localvs()].
#' @param pad Each panel's y-range is the data range widened by `pad` times
#'   its width on both sides (`pad = 2` puts the points in the middle 20%).
#' @return A `ggplot` object.
#' @export
plot_performance <- function(object, evaluation, pad = 2) {
  perf <- evaluation$performance
  prop <- lengths(selected_covariates(object))[perf$Method] / ncol(object$Xq)
  measures <- c("log(MSE)", "Hellinger Distance", "Selection Proportion")
  d <- data.frame(Method = factor(rep(perf$Method, 3), levels = names(method_labels)),
                  Measure = factor(rep(measures, each = nrow(perf)), levels = measures),
                  Plot_Value = c(log(perf$Category_MSE), perf$Category_Hellinger, prop))
  lims <- function(x) { r <- range(x, na.rm = TRUE); r + c(-1, 1) * pad * diff(r) }
  scales_y <- lapply(measures, function(ms) ggplot2::scale_y_continuous(
    limits = lims(d$Plot_Value[d$Measure == ms]), breaks = scales::pretty_breaks(5), expand = c(0, 0)))

  ggplot2::ggplot(d, ggplot2::aes(Method, Plot_Value, shape = Method)) +
    ggplot2::geom_point(size = 4, stroke = 1.1) +
    ggplot2::facet_wrap(~ Measure, ncol = 3, scales = "free_y") +
    ggh4x::facetted_pos_scales(y = scales_y) +
    ggplot2::scale_shape_manual(values = method_shapes) +
    ggplot2::labs(x = NULL, y = NULL, shape = "Method") +
    shape_theme()
}

#' Map the fitted intensity surfaces
#'
#' Plots `log(1 + intensity)` of each penalized fit at the dummy quadrature
#' points, optionally clipped to a boundary, with the events overlaid.
#'
#' @param object A `localvs` fit.
#' @param evaluation Output of [evaluate_localvs()].
#' @param boundary Optional `sf` polygon (e.g. [stl_boundary]).
#' @param events Overlay the events?
#' @return A `ggplot` object.
#' @export
plot_intensity <- function(object, evaluation, boundary = NULL, events = TRUE) {
  U <- object$quad$U
  idx <- which(object$quad$y == 0)
  if (!is.null(boundary)) {
    pts <- sf::st_as_sf(U[idx, ], coords = c("x", "y"), crs = 4326)
    idx <- idx[lengths(sf::st_intersects(pts, sf::st_union(sf::st_transform(boundary, 4326)))) > 0]
  }
  int <- evaluation$intensity[idx, , drop = FALSE]
  d <- data.frame(longitude = rep(U$x[idx], ncol(int)), latitude = rep(U$y[idx], ncol(int)),
                  Method = factor(rep(names(int), each = length(idx)), levels = names(method_labels)),
                  Intensity = unlist(int, use.names = FALSE))
  rx <- range(d$longitude); ry <- range(d$latitude)

  p <- ggplot2::ggplot(d, ggplot2::aes(longitude, latitude, fill = log1p(Intensity))) +
    ggplot2::geom_tile()
  if (!is.null(boundary))
    p <- p + ggplot2::geom_sf(data = boundary, fill = NA, linewidth = 0.45, inherit.aes = FALSE)
  if (events)
    p <- p + ggplot2::geom_point(data = object$data, ggplot2::aes(longitude, latitude), shape = 3,
                                 size = 1, stroke = 0.5, inherit.aes = FALSE)
  p + ggplot2::facet_wrap(~ Method, ncol = 3) +
    ggplot2::scale_fill_viridis_c(option = "C", name = expression(log(1 + hat(pi)(s)))) +
    ggplot2::coord_sf(expand = FALSE, xlim = rx, ylim = ry, crs = 4326) +
    ggplot2::scale_x_continuous(breaks = pretty(rx, n = 3)) +
    ggplot2::scale_y_continuous(breaks = pretty(ry, n = 4)) +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
                   axis.text = ggplot2::element_text(size = 9),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                   panel.spacing = ggplot2::unit(0.8, "lines"))
}

#' Regularization path diagnostics
#'
#' Three panels: CBIC along the path with the selected value marked, the
#' coefficient paths, and the percentage of active coefficients.
#'
#' @param object A `localvs` fit.
#' @param method Fitted method to show (default `"LLI"`).
#' @return A `patchwork` object; the attribute `summary` holds the selected
#'   lambda, total/active coefficients and active percentage.
#' @export
plot_cbic_path <- function(object, method = "LLI") {
  f <- object$fits[[method]]
  if (is.null(f)) stop(method, " was not fitted.", call. = FALSE)
  if (!is.null(f$fit)) { Bp <- as.matrix(f$fit$beta); lam <- f$fit$lambda }
  else                 { Bp <- f$beta_path;            lam <- f$lambda_path }
  tol <- object$tol; log_l <- log(f$lambda)
  idx <- which.min(abs(lam - f$lambda))
  act <- 100 * colSums(abs(Bp) > tol) / nrow(Bp)

  thm  <- ggplot2::theme_bw(base_size = 13) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  vl   <- ggplot2::geom_vline(xintercept = log_l, linetype = "dashed", linewidth = 0.8, colour = "red")
  mark <- function(x, y) ggplot2::annotate("point", x = x, y = y, shape = 4, colour = "red", size = 5, stroke = 1.5)
  lx   <- ggplot2::labs(x = expression(log(lambda)))

  p1 <- ggplot2::ggplot(f$cbic_path, ggplot2::aes(log(lambda), CBIC)) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#0072B2") + mark(log_l, f$cbic) + lx + thm

  rows <- which(apply(abs(Bp), 1, max) > tol)
  dp <- data.frame(lambda = rep(lam, each = length(rows)), Coef = rep(rows, times = length(lam)),
                   delta = as.vector(Bp[rows, , drop = FALSE]))
  p2 <- ggplot2::ggplot(dp, ggplot2::aes(log(lambda), delta, group = Coef)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.5, colour = "#D55E00") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) + vl + lx +
    ggplot2::labs(y = "Local coefficients (scaling and detail)") + thm

  da <- data.frame(log_lambda = log(lam), active = act)
  p3 <- ggplot2::ggplot(da, ggplot2::aes(log_lambda, active)) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#009E73") + vl + mark(log(lam[idx]), act[idx]) +
    lx + ggplot2::labs(y = "Active coefficients (%)") + thm

  out <- patchwork::wrap_plots(p1, p2, p3, nrow = 1)
  attr(out, "summary") <- data.frame(Method = method, lambda = f$lambda, log_lambda = log_l, CBIC = f$cbic,
                                     Total = nrow(Bp), Active = sum(abs(Bp[, idx]) > tol),
                                     Active_pct = act[idx])
  out
}

#' Ripley's K-function against complete spatial randomness
#'
#' Projects the events and boundary, keeps events inside the boundary and
#' plots the estimated K-function (km^2) with the CSR reference \eqn{\pi d^2}.
#'
#' @param data Data frame of events with longitude/latitude.
#' @param boundary `sf` polygon.
#' @param crs Projected CRS in metres (default UTM 15N).
#' @param step Distance step in metres.
#' @param correction Edge correction: `"translation"`, `"isotropic"` or `"border"`.
#' @return Invisibly, a data frame of distances and K values.
#' @export
plot_ripleys_k <- function(data, boundary, crs = 32615, step = 50,
                           correction = c("translation", "isotropic", "border")) {
  correction <- match.arg(correction)
  pts <- sf::st_transform(sf::st_as_sf(standardize_coords(data), coords = c("longitude", "latitude"), crs = 4326), crs)
  win <- sf::st_make_valid(sf::st_union(sf::st_transform(boundary, crs)))
  ins <- sf::st_within(pts, win, sparse = FALSE)[, 1]
  if (!any(ins)) stop("No points were found inside the boundary.", call. = FALSE)
  xy <- sf::st_coordinates(pts[ins, ]); Wk <- spatstat.geom::as.owin(win)
  r  <- seq(0, 0.25 * min(diff(Wk$xrange), diff(Wk$yrange)), by = step)
  K  <- spatstat.explore::Kest(spatstat.geom::ppp(xy[, 1], xy[, 2], window = Wk, checkdup = FALSE),
                               r = r, correction = correction)
  col <- c(translation = "trans", isotropic = "iso", border = "border")[[correction]]
  res <- data.frame(r_km = K$r / 1000, empirical_K = K[[col]] / 1e6, theoretical_K = pi * K$r^2 / 1e6)
  res$difference <- res$empirical_K - res$theoretical_K
  graphics::plot(res$r_km, res$empirical_K, type = "l", lwd = 2, xlab = "Distance d (km)",
                 ylab = expression(K(d) ~ "(km"^2 * ")"))
  graphics::lines(res$r_km, res$theoretical_K, col = "red", lwd = 2, lty = 2)
  graphics::legend("topleft", c("Empirical K(d)", "Theoretical K(d)"), col = c("black", "red"),
                   lty = 1:2, lwd = 2, bty = "n", cex = 0.85, y.intersp = 1.5)
  invisible(res)
}
