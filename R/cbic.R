#' Composite BIC for a Poisson point-process fit
#'
#' \deqn{CBIC(\lambda) = -2\,\ell_{PP}(\hat\beta_\lambda) + K_0(\lambda)\log N(R),}
#' with \eqn{\ell_{PP} = \sum_m [y_m \eta_m - w_m \exp(\eta_m)]},
#' \eqn{w_m = \exp(\mathrm{offset}_m)}, \eqn{K_0} the number of non-zero
#' penalized coefficients (intercept excluded) and \eqn{N(R)} replaced by the
#' observed event count. Terms constant in \eqn{\lambda} are omitted.
#'
#' @param y 0/1 quadrature indicator.
#' @param X Design matrix.
#' @param intercept,beta Fitted intercept and coefficients.
#' @param offset Offset, `log(w)`.
#' @param tol Threshold for counting a coefficient as non-zero.
#' @param mu Event count used in the penalty.
#' @return A list with `cbic`, `loglik_pp`, `K0` and `mu_hat`.
#' @export
cbic_pp <- function(y, X, intercept, beta, offset, tol = 1e-8, mu = sum(y)) {
  if (!is.finite(mu) || mu <= 0) stop("CBIC requires at least one observed event.", call. = FALSE)
  eta <- pmin(pmax(drop(intercept + X %*% beta), -30), 30)
  ll  <- sum(y * eta - exp(offset) * exp(eta))
  K0  <- sum(abs(beta) > tol)
  list(cbic = -2 * ll + K0 * log(mu), loglik_pp = ll, K0 = K0, mu_hat = mu)
}

path_table <- function(res, lambda) {
  data.frame(lambda = lambda, CBIC = vapply(res, `[[`, numeric(1), "cbic"),
             K0 = vapply(res, `[[`, numeric(1), "K0"),
             LogLik_PP = vapply(res, `[[`, numeric(1), "loglik_pp"))
}

pois_glmnet <- function(X, y, offset, ...) {
  glmnet::glmnet(X, y, family = "poisson", alpha = 1, offset = offset, standardize = FALSE, ...)
}

check_xy <- function(X, y, offset) {
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y).", call. = FALSE)
  if (length(offset) != length(y)) stop("offset must have length equal to y.", call. = FALSE)
}

#' Penalized Poisson fits tuned by CBIC
#'
#' * `fit_lasso_cbic()`: LASSO (optionally weighted through `penalty.factor`).
#' * `fit_adaptive_lasso_cbic()`: adaptive LASSO; a cross-validated ridge fit
#'   only builds the weights \eqn{1/(|b_0| + \epsilon)^\gamma}. Its folds are
#'   random, so call [set.seed()] first for reproducible results.
#' * `fit_scad_cbic()`: SCAD by local linear approximation, warm-started at
#'   the LASSO solution for each \eqn{\lambda}, keeping the Berman-Turner
#'   offset. `glmnet` rescales `penalty.factor` to sum to the number of
#'   variables, so \eqn{\lambda} is rescaled by `sum(pf) / p`.
#'
#' @param X Design matrix.
#' @param y 0/1 quadrature indicator.
#' @param offset Offset, `log(w)`.
#' @param penalty.factor Per-coefficient penalty weights.
#' @param tol,coef_tol Threshold for counting a coefficient as non-zero.
#' @param gamma,eps Adaptive-weight exponent and stabilizer.
#' @param a SCAD shape parameter.
#' @param nlambda Length of the SCAD lambda path.
#' @param max_iter,conv_tol Maximum LLA iterations and convergence tolerance.
#' @param mc.cores Cores for [parallel::mclapply()] (forced to 1 on Windows).
#' @return A list with the selected `intercept`, `beta`, `lambda`, `cbic`,
#'   `loglik_pp`, `K0`, `mu_hat` and the whole `cbic_path`. LASSO fits also
#'   return the `glmnet` object as `fit`; SCAD fits return `beta_path`,
#'   `intercept_path` and `lambda_path`.
#' @name cbic_fits
NULL

#' @rdname cbic_fits
#' @export
fit_lasso_cbic <- function(X, y, offset, penalty.factor = rep(1, ncol(X)), tol = 1e-8) {
  X <- as.matrix(X); check_xy(X, y, offset)
  if (length(penalty.factor) != ncol(X)) stop("penalty.factor must have length ncol(X).", call. = FALSE)
  fit <- pois_glmnet(X, y, offset, penalty.factor = penalty.factor)
  cf  <- as.matrix(stats::coef(fit))
  res <- lapply(seq_along(fit$lambda), function(k) cbic_pp(y, X, cf[1, k], cf[-1, k], offset, tol))
  k   <- which.min(vapply(res, `[[`, numeric(1), "cbic"))
  c(list(intercept = unname(cf[1, k]), beta = as.numeric(cf[-1, k]), lambda = fit$lambda[k]),
    res[[k]], list(cbic_path = path_table(res, fit$lambda), fit = fit))
}

#' @rdname cbic_fits
#' @export
fit_adaptive_lasso_cbic <- function(X, y, offset, gamma = 1, eps = 1e-6, tol = 1e-8) {
  X <- as.matrix(X); check_xy(X, y, offset)
  ridge <- glmnet::cv.glmnet(X, y, family = "poisson", alpha = 0, offset = offset, standardize = FALSE)
  b0 <- as.numeric(stats::coef(ridge, s = "lambda.min"))[-1]
  fit_lasso_cbic(X, y, offset, penalty.factor = 1 / (abs(b0) + eps)^gamma, tol = tol)
}

#' @rdname cbic_fits
#' @export
fit_scad_cbic <- function(X, y, offset, a = 3.7, nlambda = 30, max_iter = 20, conv_tol = 1e-6,
                          coef_tol = 1e-8, mc.cores = max(1, parallel::detectCores() - 1)) {
  X <- as.matrix(X); check_xy(X, y, offset)
  p <- ncol(X)
  init <- pois_glmnet(X, y, offset, nlambda = nlambda)
  cf0  <- as.matrix(stats::coef(init)); lam_seq <- init$lambda

  one_lambda <- function(i) {
    lam <- lam_seq[i]; cf <- unname(cf0[, i]); conv <- FALSE; it <- 0L
    for (it in seq_len(max_iter)) {
      ab <- abs(cf[-1])
      pf <- ifelse(ab <= lam, 1, pmax(a * lam - ab, 0) / ((a - 1) * lam))  # p'_lambda(|b|) / lambda
      pf[!is.finite(pf)] <- 1
      new <- if (sum(pf) <= .Machine$double.eps) {
        b <- as.numeric(stats::glm.fit(cbind(1, X), y, family = stats::poisson(), offset = offset)$coefficients)
        replace(b, !is.finite(b), 0)
      } else {
        as.numeric(stats::coef(pois_glmnet(X, y, offset, lambda = lam * sum(pf) / p,
                                           penalty.factor = pf, thresh = 1e-7)))
      }
      d <- max(abs(new - cf), na.rm = TRUE); cf <- new
      if (is.finite(d) && d < conv_tol) { conv <- TRUE; break }
    }
    c(list(intercept = cf[1], beta = cf[-1], lambda = lam),
      cbic_pp(y, X, cf[1], cf[-1], offset, coef_tol), list(converged = conv, iterations = it))
  }

  if (.Platform$OS.type == "windows") mc.cores <- 1
  res  <- parallel::mclapply(seq_along(lam_seq), one_lambda, mc.cores = mc.cores)
  best <- res[[which.min(vapply(res, `[[`, numeric(1), "cbic"))]]
  best$cbic_path <- cbind(path_table(res, lam_seq),
                          Converged  = vapply(res, `[[`, logical(1), "converged"),
                          Iterations = vapply(res, `[[`, integer(1), "iterations"))
  best$beta_path      <- vapply(res, `[[`, numeric(p), "beta")
  best$intercept_path <- vapply(res, `[[`, numeric(1), "intercept")
  best$lambda_path    <- lam_seq
  best
}
