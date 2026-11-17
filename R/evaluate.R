#' Evaluate fitted models
#'
#' * In-sample negative log-likelihoods of the penalized fits.
#' * Post-selection comparison: the selected design of each method is
#'   standardized and refitted by unpenalized Poisson regression; fitted
#'   intensities are aggregated to grid-cell counts and the frequency
#'   distribution of those counts is compared with the empirical one by
#'   mean squared error and Hellinger distance.
#'
#' @param object A `localvs` fit.
#' @return A list with `nll` (penalized fits), `performance` (post-selection
#'   MSE and Hellinger distance with ranks), `grid` (observed and predicted
#'   counts per grid cell) and `intensity` (penalized intensity at every
#'   quadrature point).
#' @export
evaluate_localvs <- function(object) {
  q <- object$quad; y <- q$y; off <- q$offset; m_names <- names(object$fits)
  design <- function(m) if (m %in% local_methods) object$Z else object$Xq

  eta <- lapply(m_names, function(m) object$fits[[m]]$intercept + drop(design(m) %*% object$coefficients[[m]]))
  names(eta) <- m_names
  intensity <- as.data.frame(lapply(eta, function(e) pmax(exp(e), 1e-12)))
  nll <- data.frame(
    Method            = m_names,
    Number_Selected   = vapply(object$coefficients, function(b) sum(abs(b) > object$tol), numeric(1)),
    Poisson_NLL       = vapply(eta, function(e) -sum(stats::dpois(y, pmax(exp(e + off), 1e-12), log = TRUE)), numeric(1)),
    Point_Process_NLL = vapply(intensity, function(l) -sum(y * log(l) - exp(off) * l), numeric(1)),
    row.names = NULL)

  ## Post-selection refit
  refit <- lapply(m_names, function(m) {
    X <- if (m %in% local_methods) selected_local_design(object, m)
         else object$Xq[, object$selection[[m]]$covariate, drop = FALSE]
    refit_intensity(rescale_design(X), y, off)
  })
  names(refit) <- m_names

  id <- grid_cell_id(q)
  counts <- c(list(Empirical_Count = y), lapply(refit, function(l) q$w * l))
  grid <- data.frame(grid_id = sort(unique(id)),
                     lapply(counts, function(v) as.numeric(tapply(v, id, sum))), check.names = FALSE)

  cats  <- sort(unique(c(grid$Empirical_Count, round(unlist(grid[m_names])))))
  freq  <- function(v) as.vector(table(factor(v, levels = cats)))
  n_obs <- freq(grid$Empirical_Count); p_obs <- n_obs / sum(n_obs)
  perf <- do.call(rbind, lapply(m_names, function(m) {
    n_pred <- freq(round(grid[[m]])); p_pred <- n_pred / sum(n_pred)
    data.frame(Method = m, Category_MSE = mean((n_obs - n_pred)^2),
               Category_Hellinger = sqrt(sum((sqrt(p_obs) - sqrt(p_pred))^2)) / sqrt(2))
  }))
  perf$Rank_MSE       <- rank(perf$Category_MSE, ties.method = "min")
  perf$Rank_Hellinger <- rank(perf$Category_Hellinger, ties.method = "min")

  list(nll = nll, performance = perf, grid = grid, intensity = intensity)
}

selected_local_design <- function(object, method) {
  s <- object$selection[[method]]
  if (!nrow(s)) return(matrix(numeric(0), nrow(object$Z), 0))
  n_atoms <- nrow(object$meta)
  p  <- match(s$covariate, colnames(object$Xq)); a <- as.integer(s$atom_index)
  ok <- is.finite(p) & is.finite(a) & a >= 1 & a <= n_atoms
  object$Z[, unique((p[ok] - 1L) * n_atoms + a[ok]), drop = FALSE]
}

rescale_design <- function(X) {
  if (!ncol(X)) return(X)
  X <- X[, apply(X, 2, function(z) all(is.finite(z)) && isTRUE(stats::sd(z) > 1e-12)), drop = FALSE]
  if (!ncol(X)) return(X)
  Xs <- scale(X); Xs[!is.finite(Xs)] <- 0
  Xs
}

refit_intensity <- function(X, y, offset) {
  Xr  <- cbind(1, X)
  fit <- suppressWarnings(stats::glm.fit(Xr, y, family = stats::poisson(), offset = offset,
                                         control = stats::glm.control(maxit = 200, epsilon = 1e-8)))
  b <- as.numeric(fit$coefficients); b[!is.finite(b)] <- 0
  exp(pmin(pmax(drop(Xr %*% b), -30), 30))
}

#' Microscale feature selection summary
#'
#' Counts selected covariates whose name ends in one of the given radii
#' (e.g. `buildings_50m`).
#'
#' @param object A `localvs` fit.
#' @param radii Radii (metres) treated as microscale.
#' @return A data frame per method with the attribute `variables`, a list of
#'   the selected microscale covariates.
#' @export
microscale_summary <- function(object, radii = c(50, 100, 200)) {
  sel <- selected_covariates(object)
  is_micro <- function(v) grepl("_[0-9]+m$", v) &
    suppressWarnings(as.numeric(sub(".*_([0-9]+)m$", "\\1", v))) %in% radii
  micro <- lapply(sel, function(v) v[is_micro(v)])
  out <- data.frame(Method = names(sel), Total_Selected = lengths(sel),
                    Microscale_Selected = lengths(micro), row.names = NULL)
  out$Microscale_Percentage <- 100 * out$Microscale_Selected / out$Total_Selected
  out <- out[order(-out$Microscale_Selected, -out$Microscale_Percentage), ]
  rownames(out) <- NULL
  attr(out, "variables") <- micro
  out
}
