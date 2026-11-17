## Final real-data analysis script behind the revised manuscript.
## The package functions (fit_localvs() etc.) reproduce the same pipeline.

rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  for (p in c("spatstat.geom", "spatstat.explore", "spatstat.model", "glmnet", "tidyverse", "sf",
              "viridis", "leaflet", "htmltools", "htmlwidgets", "patchwork", "ggh4x", "tigris"))
    library(p, character.only = TRUE)
})
options(tigris_use_cache = TRUE)
set.seed(2026)  # cv.glmnet (adaptive-LASSO ridge weights) uses random folds

## ================================================================
## 0. Settings
## ================================================================
data_dir      <- "data-raw"   # run from the package root; raw CSV and boundary live in data-raw/
crime_file    <- file.path(data_dir, "stl_high_dimensional_crime_features.csv")
boundary_file <- file.path(data_dir, "stl_boundary.shp")
qtiles  <- c(32, 32)
haar_J  <- 3
tol0    <- 1e-8
use_inside_only <- FALSE   # original script analysed ALL points; set TRUE to drop points outside the boundary
methods <- c(GLI = "Global LASSO", GSI = "Global SCAD", GALI = "Adaptive LASSO",
             LLI = "Local LASSO",  LSI = "Local SCAD")
shapes  <- c(GLI = 1, GSI = 2, GALI = 3, LLI = 4, LSI = 5)

## ================================================================
## 1. Data and boundary
## ================================================================
std_coords <- function(d) {
  lat <- names(d)[tolower(names(d)) %in% c("latitude", "lat")][1]
  lon <- names(d)[tolower(names(d)) %in% c("longitude", "lon", "lng", "long")][1]
  if (is.na(lat) || is.na(lon)) stop("Latitude and Longitude columns not found.")
  d |> rename(latitude = all_of(lat), longitude = all_of(lon)) |>
    filter(is.finite(latitude), is.finite(longitude))
}
read_boundary <- function(path) {
  Sys.setenv(SHAPE_RESTORE_SHX = "YES")
  b <- st_read(path, quiet = TRUE)
  if (is.na(st_crs(b))) st_crs(b) <- 4326
  st_transform(st_make_valid(b), 4326)
}

df <- read.csv(crime_file, check.names = FALSE)
names(df) <- str_replace_all(str_trim(names(df)), "\\s+", "_")
df <- df[, -(164:168)] |> distinct(Longitude, Latitude, .keep_all = TRUE) |> std_coords()

stl_boundary <- read_boundary(boundary_file)
inside <- lengths(st_intersects(st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326),
                                st_union(stl_boundary))) > 0
cat("Points outside boundary:", sum(!inside), "\n"); print(df[!inside, c("longitude", "latitude")])
df_in <- df[inside, ]

## Leaflet helpers (reused by all maps)
add_plus <- function(m, lng, lat, size = 16)
  addLabelOnlyMarkers(m, lng = lng, lat = lat, label = "+", labelOptions = labelOptions(
    noHide = TRUE, direction = "center", textOnly = TRUE,
    style = list(color = "black", `font-size` = paste0(size, "px"), `font-weight` = "bold")))
add_city <- function(m, b, weight = 3)
  addPolygons(m, data = b, color = "black", weight = weight, opacity = 1,
              fill = TRUE, fillColor = "transparent", fillOpacity = 0)
fit_city <- function(m, b) {
  bb <- st_bbox(b); fitBounds(m, bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
}

overview_map <- leaflet() |> addProviderTiles(providers$OpenStreetMap) |>
  addPolygons(data = stl_boundary, color = "red", weight = 3, fill = FALSE) |>
  add_plus(df_in$longitude, df_in$latitude) |>
  setView(lng = -90.1994, lat = 38.6270, zoom = 11)
overview_map

if (use_inside_only) df <- df_in
cat("Observations used in the analysis:", nrow(df), "\n")

## Covariates: numeric, non-constant, non-finite values -> median
drop_cols <- c("Type", "Reference", "Road_Name", "Source", "Occurance_Time", "latitude",
               "longitude", "Latitude", "Longitude", "crime_id", "geometry")
Xdf <- df[setdiff(names(df), drop_cols)] |> select(where(is.numeric))
Xdf <- Xdf[vapply(Xdf, \(z) sum(is.finite(z)) > 5 && isTRUE(sd(z, na.rm = TRUE) > 0), logical(1))]
Xdf[] <- lapply(Xdf, \(z) replace(z, !is.finite(z), median(z[is.finite(z)])))

s_n <- as.matrix(df[c("longitude", "latitude")])
pad <- function(r, f = 0.01) r + c(-1, 1) * f * diff(r)
W   <- owin(pad(range(s_n[, 1])), pad(range(s_n[, 2])))
Xpp <- ppp(s_n[, 1], s_n[, 2], window = W, checkdup = FALSE)

## ================================================================
## 2. Ripley's K-function
## ================================================================
plot_ripleys_k <- function(df, boundary, crs = 32615, step = 50, correction = "translation") {
  pts <- st_as_sf(std_coords(df), coords = c("longitude", "latitude"), crs = 4326) |> st_transform(crs)
  win <- st_make_valid(st_union(st_transform(boundary, crs)))
  ins <- st_within(pts, win, sparse = FALSE)[, 1]
  cat("Points inside boundary:", sum(ins), "| outside:", sum(!ins), "\n")
  if (!any(ins)) stop("No points were found inside the specified boundary.")
  xy <- st_coordinates(pts[ins, ]); Wk <- as.owin(win)
  r  <- seq(0, 0.25 * min(diff(Wk$xrange), diff(Wk$yrange)), by = step)
  K  <- Kest(ppp(xy[, 1], xy[, 2], window = Wk, checkdup = FALSE), r = r, correction = correction)
  col <- switch(correction, translation = "trans", isotropic = "iso", border = "border", "trans")
  res <- data.frame(r_km = K$r / 1000, empirical_K = K[[col]] / 1e6, theoretical_K = pi * K$r^2 / 1e6)
  res$difference <- res$empirical_K - res$theoretical_K
  plot(res$r_km, res$empirical_K, type = "l", lwd = 2, xlab = "Distance d (km)",
       ylab = expression(K(d) ~ "(km"^2 * ")"))
  lines(res$r_km, res$theoretical_K, col = "red", lwd = 2, lty = 2)
  legend("topleft", c("Empirical K(d)", "Theoretical K(d)"), col = c("black", "red"),
         lty = 1:2, lwd = 2, bty = "n", cex = 0.85, y.intersp = 1.5)
  invisible(res)
}
plot_ripleys_k(df, stl_boundary)

## ================================================================
## 3. Berman-Turner quadrature
## ================================================================
Q   <- quadscheme(Xpp, method = "grid", ntile = qtiles)
Upp <- as.ppp(Q); U <- data.frame(x = Upp$x, y = Upp$y); M <- nrow(U)
yi  <- as.integer(as.logical(is.data(Q)))
wi  <- as.numeric(weights(Q))
cell_area <- diff(W$xrange) * diff(W$yrange) / prod(qtiles)
if (length(wi) != M) wi <- rep(cell_area, M)
if (length(yi) != M) yi <- c(rep(1L, npoints(Xpp)), rep(0L, M - npoints(Xpp)))
wi[!is.finite(wi) | wi <= 0] <- cell_area
off <- log(pmax(wi, 1e-12))
cat("Quadrature points:", M, "| Events:", sum(yi), "| Dummy:", sum(yi == 0), "\n")

## ================================================================
## 4. Local Haar basis at level j = J - 1
##    Each of the 4^(J-1) blocks carries one scaling atom (the level of beta_i on
##    the block) and three detail atoms (horizontal, vertical, diagonal contrasts):
##    4 * 4^(J-1) = 4^J atoms in total, every one supported on a single block.
##    Together they span all surfaces that are constant on the 2^J x 2^J sub-blocks.
## ================================================================
to_unit <- function(s) {
  u <- apply(as.matrix(s), 2, \(v) (v - min(v)) / max(diff(range(v)), 1e-12))
  pmin(u, 1 - 1e-12)                     # points on the upper edge stay in the last block
}
haar_phi <- function(j, k, x) 2^(j / 2) * (x >= k / 2^j & x < (k + 1) / 2^j)
haar_psi <- function(j, k, x) {
  m <- (k + 0.5) / 2^j
  2^(j / 2) * ((x >= k / 2^j & x < m) - (x >= m & x < (k + 1) / 2^j))
}

build_local_haar <- function(s01, J = 3, xy = NULL) {
  j <- J - 1
  meta <- expand.grid(type = c("phi", "psi_x", "psi_y", "psi_xy"), ky = 0:(2^j - 1), kx = 0:(2^j - 1),
                      stringsAsFactors = FALSE) |>
    mutate(j = j) |> select(j, type, kx, ky) |>
    mutate(atom_index = row_number(), x_lo01 = kx / 2^j, x_hi01 = (kx + 1) / 2^j,
           y_lo01 = ky / 2^j, y_hi01 = (ky + 1) / 2^j)
  B <- mapply(function(type, kx, ky) {
    fx <- if (type %in% c("psi_x", "psi_xy")) haar_psi else haar_phi
    fy <- if (type %in% c("psi_y", "psi_xy")) haar_psi else haar_phi
    fx(j, kx, s01[, 1]) * fy(j, ky, s01[, 2])
  }, meta$type, meta$kx, meta$ky)
  if (!is.null(xy)) {
    xr <- range(xy[, 1]); yr <- range(xy[, 2])
    meta <- mutate(meta, x_lo = xr[1] + x_lo01 * diff(xr), x_hi = xr[1] + x_hi01 * diff(xr),
                   y_lo = yr[1] + y_lo01 * diff(yr), y_hi = yr[1] + y_hi01 * diff(yr))
  }
  list(B = B, meta = meta)
}

haar      <- build_local_haar(to_unit(U), J = haar_J, xy = as.matrix(U))
B_fine    <- haar$B                      # names kept so that the rest of the script is unchanged
meta_fine <- haar$meta

## checks: 4^J atoms; every quadrature point lies in exactly one block (one scaling atom)
stopifnot(ncol(B_fine) == 4^haar_J,
          all(rowSums(B_fine[, meta_fine$type == "phi", drop = FALSE] != 0) == 1))
cat(sprintf("Haar basis: J = %d, %d blocks, %d atoms (%d scaling + %d detail)\n",
            haar_J, 4^(haar_J - 1), ncol(B_fine), sum(meta_fine$type == "phi"),
            sum(meta_fine$type != "phi")))

## ================================================================
## 5. Smooth covariates to quadrature points; localized design Z
## ================================================================
eval_im <- function(im, x, y) {
  np <- nearest.pixel(x, y, im)
  im$v[cbind(pmin(pmax(np$row, 1L), nrow(im$v)), pmin(pmax(np$col, 1L), ncol(im$v)))]
}
Xq_raw <- sapply(Xdf, \(v) eval_im(Smooth(`marks<-`(Xpp, value = v), sigma = 0.01, at = "pixels",
                                          normalise = TRUE, dimyx = qtiles), U$x, U$y))
Xq_raw[!is.finite(Xq_raw)] <- 0
Xq <- scale(Xq_raw); Xq[!is.finite(Xq)] <- 0

basis_names <- with(meta_fine, paste0(type, "_j", j, "_kx", kx, "_ky", ky))
Z <- do.call(cbind, lapply(seq_len(ncol(Xq)), \(p) B_fine * Xq[, p]))
colnames(Z) <- paste(rep(colnames(Xq), each = ncol(B_fine)), basis_names, sep = "__")
ok_col <- function(X) apply(X, 2, \(z) all(is.finite(z)) && sd(z) > 0)
keep_Z <- ok_col(Z); Zc <- Z[, keep_Z, drop = FALSE]
keep_X <- ok_col(Xq); Xg <- Xq[, keep_X, drop = FALSE]

## ================================================================
## 6. Penalized Poisson fits selected by CBIC
##    CBIC = -2 * ell_PP + K0 * log N(R),
##    ell_PP = sum[y*eta - w*exp(eta)], w = exp(offset); K0 excludes intercept.
## ================================================================
cbic_pp <- function(y, X, a0, beta, offset, tol = 1e-8, mu = sum(y)) {
  if (!is.finite(mu) || mu <= 0) stop("CBIC requires at least one observed event.")
  eta <- pmin(pmax(drop(a0 + X %*% beta), -30), 30)
  ll  <- sum(y * eta - exp(offset) * exp(eta)); K0 <- sum(abs(beta) > tol)
  list(cbic = -2 * ll + K0 * log(mu), loglik_pp = ll, K0 = K0, mu_hat = mu)
}
path_df <- function(res, lambda) data.frame(lambda = lambda, CBIC = sapply(res, `[[`, "cbic"),
                                            K0 = sapply(res, `[[`, "K0"), LogLik_PP = sapply(res, `[[`, "loglik_pp"))
pois_glmnet <- function(X, y, offset, ...)
  glmnet(X, y, family = "poisson", alpha = 1, offset = offset, standardize = FALSE, ...)

fit_lasso_cbic <- function(X, y, offset, penalty.factor = rep(1, ncol(X)), tol = 1e-8) {
  stopifnot(nrow(X) == length(y), length(offset) == length(y), length(penalty.factor) == ncol(X))
  fit <- pois_glmnet(X, y, offset, penalty.factor = penalty.factor)
  cf  <- as.matrix(coef(fit))
  res <- lapply(seq_along(fit$lambda), \(k) cbic_pp(y, X, cf[1, k], cf[-1, k], offset, tol))
  k   <- which.min(sapply(res, `[[`, "cbic"))
  c(list(intercept = cf[1, k], beta = as.numeric(cf[-1, k]), lambda = fit$lambda[k]), res[[k]],
    list(cbic_path = path_df(res, fit$lambda), fit = fit))
}

## Adaptive LASSO: ridge (CV) only builds the weights; final lambda chosen by CBIC
fit_adaptive_lasso_cbic <- function(X, y, offset, gamma = 1, eps = 1e-6, tol = 1e-8) {
  ridge <- cv.glmnet(X, y, family = "poisson", alpha = 0, offset = offset, standardize = FALSE)
  b0 <- as.numeric(coef(ridge, s = "lambda.min"))[-1]
  fit_lasso_cbic(X, y, offset, penalty.factor = 1 / (abs(b0) + eps)^gamma, tol = tol)
}

## SCAD via LLA with the BT offset (warm start = LASSO at the same lambda).
## glmnet rescales penalty.factor to sum to p, so lambda is rescaled by sum(pf)/p.
fit_scad_cbic <- function(X, y, offset, a = 3.7, nlambda = 30, max_iter = 20, tol = 1e-6,
                          coef_tol = 1e-8, mc.cores = max(1, parallel::detectCores() - 1)) {
  stopifnot(nrow(X) == length(y), length(offset) == length(y))
  p <- ncol(X)
  init <- pois_glmnet(X, y, offset, nlambda = nlambda)
  cf0 <- as.matrix(coef(init)); lam_seq <- init$lambda
  one_lambda <- function(i) {
    lam <- lam_seq[i]; cf <- unname(cf0[, i]); conv <- FALSE
    for (it in seq_len(max_iter)) {
      ab <- abs(cf[-1])
      pf <- ifelse(ab <= lam, 1, pmax(a * lam - ab, 0) / ((a - 1) * lam))  # p'_lambda(|b|) / lambda
      pf[!is.finite(pf)] <- 1
      new <- if (sum(pf) <= .Machine$double.eps) {
        b <- as.numeric(glm.fit(cbind(1, X), y, family = poisson(), offset = offset)$coefficients)
        replace(b, !is.finite(b), 0)
      } else as.numeric(coef(pois_glmnet(X, y, offset, lambda = lam * sum(pf) / p,
                                         penalty.factor = pf, thresh = 1e-7)))
      d <- max(abs(new - cf), na.rm = TRUE); cf <- new
      if (is.finite(d) && d < tol) { conv <- TRUE; break }
    }
    c(list(intercept = cf[1], beta = cf[-1], lambda = lam),
      cbic_pp(y, X, cf[1], cf[-1], offset, coef_tol), list(converged = conv, iterations = it))
  }
  if (.Platform$OS.type == "windows") mc.cores <- 1
  res  <- parallel::mclapply(seq_along(lam_seq), one_lambda, mc.cores = mc.cores)
  best <- res[[which.min(sapply(res, `[[`, "cbic"))]]
  best$cbic_path <- cbind(path_df(res, lam_seq), Converged = sapply(res, `[[`, "converged"),
                          Iterations = sapply(res, `[[`, "iterations"))
  best$beta_path      <- sapply(res, `[[`, "beta")
  best$intercept_path <- sapply(res, `[[`, "intercept")
  best$lambda_path    <- lam_seq
  best
}

## ================================================================
## 7. Fit all five methods
## ================================================================
fits <- list(LLI  = fit_lasso_cbic(Zc, yi, off),
             LSI  = fit_scad_cbic(Zc, yi, off),
             GLI  = fit_lasso_cbic(Xg, yi, off),
             GSI  = fit_scad_cbic(Xg, yi, off),
             GALI = fit_adaptive_lasso_cbic(Xg, yi, off))[names(methods)]
is_local <- setNames(names(fits) %in% c("LLI", "LSI"), names(fits))

## Coefficients on the full (unfiltered) design
coefs <- imap(fits, \(f, m) {
  keep <- if (is_local[m]) keep_Z else keep_X
  b <- setNames(numeric(length(keep)), if (!is_local[m]) colnames(Xq)); b[keep] <- f$beta; b
})

map_selection <- function(beta, meta, cov_names, tol = 1e-8) {
  G  <- matrix(beta, nrow = nrow(meta), ncol = length(cov_names))
  nz <- which(abs(G) > tol, arr.ind = TRUE)
  if (!nrow(nz)) return(data.frame())
  out <- cbind(covariate = cov_names[nz[, 2]], coefficient = G[nz], meta[nz[, 1], , drop = FALSE])
  out[order(out$covariate, out$j, out$type), ]
}
sel_local <- lapply(coefs[c("LLI", "LSI")], map_selection, meta = meta_fine, cov_names = colnames(Xq))
sel_vars  <- c(lapply(coefs[c("GLI", "GSI", "GALI")], \(b) names(b)[abs(b) > tol0]),
               lapply(sel_local, \(s) unique(s$covariate)))[names(methods)]

## ================================================================
## 8. Penalized in-sample evaluation
## ================================================================
eta <- imap(fits, \(f, m) f$intercept + drop((if (is_local[m]) Z else Xq) %*% coefs[[m]]))
intensity_penalized <- as.data.frame(lapply(eta, \(e) pmax(exp(e), 1e-12)))
nll_table <- tibble(
  Method            = names(fits),
  Number_Selected   = sapply(coefs, \(b) sum(abs(b) > tol0)),
  Poisson_NLL       = sapply(eta, \(e) -sum(dpois(yi, pmax(exp(e + off), 1e-12), log = TRUE))),
  Point_Process_NLL = sapply(intensity_penalized, \(l) -sum(yi * log(l) - exp(off) * l)))
print(nll_table)

## ================================================================
## 9. Post-selection refit and grid-count comparison (MSE, Hellinger)
## ================================================================
local_design <- function(sel, n_atoms = ncol(B_fine), pred = colnames(Xq)) {
  if (is.null(sel) || !nrow(sel)) return(matrix(numeric(0), nrow(Z), 0))
  p <- match(sel$covariate, pred); a <- as.integer(sel$atom_index)
  ok <- is.finite(p) & is.finite(a) & a >= 1 & a <= n_atoms
  Z[, unique((p[ok] - 1L) * n_atoms + a[ok]), drop = FALSE]
}
rescale_design <- function(X) {
  if (!ncol(X)) return(X)
  X <- X[, apply(X, 2, \(z) all(is.finite(z)) && isTRUE(sd(z) > 1e-12)), drop = FALSE]
  if (!ncol(X)) return(X)
  Xs <- scale(X); Xs[!is.finite(Xs)] <- 0; Xs
}
refit_intensity <- function(X, y, offset) {
  Xr <- cbind(1, X)
  fit <- suppressWarnings(glm.fit(Xr, y, family = poisson(), offset = offset,
                                  control = glm.control(maxit = 200, epsilon = 1e-8)))
  b <- as.numeric(fit$coefficients); b[!is.finite(b)] <- 0
  exp(pmin(pmax(drop(Xr %*% b), -30), 30))
}

X_designs <- imap(sel_vars, \(v, m) if (is_local[m]) local_design(sel_local[[m]])
                  else Xq[, v, drop = FALSE]) |> lapply(rescale_design)
lambdas <- lapply(X_designs, refit_intensity, y = yi, offset = off)

cell <- function(v, r, n) findInterval(v, seq(r[1], r[2], length.out = n + 1), all.inside = TRUE)
grid_pred <- data.frame(grid_id = (cell(U$y, W$yrange, qtiles[2]) - 1L) * qtiles[1] +
                          cell(U$x, W$xrange, qtiles[1]),
                        Empirical_Count = yi, lapply(lambdas, \(l) wi * l)) |>
  group_by(grid_id) |> summarise(across(everything(), sum))

all_cats <- sort(unique(c(grid_pred$Empirical_Count, round(unlist(grid_pred[names(fits)])))))
freq  <- function(v) as.vector(table(factor(v, levels = all_cats)))
n_obs <- freq(grid_pred$Empirical_Count); p_obs <- n_obs / sum(n_obs)

plot_table <- map_dfr(names(fits), \(m) {
  n_pred <- freq(round(grid_pred[[m]])); p_pred <- n_pred / sum(n_pred)
  tibble(Method = m, Category_MSE = mean((n_obs - n_pred)^2),
         Category_Hellinger = sqrt(sum((sqrt(p_obs) - sqrt(p_pred))^2)) / sqrt(2))
}) |> mutate(Rank_MSE = rank(Category_MSE, ties.method = "min"),
             Rank_Hellinger = rank(Category_Hellinger, ties.method = "min"))
print(as.data.frame(plot_table), row.names = FALSE)

## ================================================================
## 10. Microscale feature selection (50/100/200 m)
## ================================================================
micro_radii <- c(50, 100, 200)
is_micro <- function(v) suppressWarnings(as.numeric(sub(".*_([0-9]+)m$", "\\1", v))) %in% micro_radii &
  grepl("_[0-9]+m$", v)
micro_vars <- lapply(sel_vars, \(v) sort(v[is_micro(v)]))
micro_summary <- tibble(Method = names(sel_vars), Total_Selected = lengths(sel_vars),
                        Microscale_Selected = lengths(micro_vars)) |>
  mutate(Microscale_Percentage = 100 * Microscale_Selected / Total_Selected) |>
  arrange(desc(Microscale_Selected), desc(Microscale_Percentage))
print(as.data.frame(micro_summary), row.names = FALSE)
iwalk(micro_vars, \(v, m) cat(m, ": ", if (length(v)) paste(v, collapse = ", ") else "None", "\n", sep = ""))
## ================================================================
## 10b. Variable-selection table (check / cross), for verification
##      Paste this block right after Section 10 (it only needs sel_vars,
##      methods, data_dir and the tidyverse, all defined above).
##      A covariate counts as selected by a local method (LLI/LSI) if ANY of
##      its Haar atoms has a non-zero coefficient -- the same rule as sel_vars.
## ================================================================

## Category keywords in display order; within a category, rows follow the
## keyword order and then the radius (50, 100, 200, 500, 1000 m).
cat_keys <- list(
  "Road Types"           = c("intersection", "motorway", "trunk", "primary", "residential",
                             "secondary", "tertiary", "service", "road|street|highway|junction"),
  "Buildings / Land Use" = c("building", "commercial", "industrial", "retail", "\\bbars?\\b",
                             "\\bpubs?\\b", "nightclub", "restaurant", "fast food", "cafe",
                             "cinema", "convenience", "supermarket", "fuel"),
  "Transportation"       = c("\\bbus", "parking", "rail|tram|subway|transit"),
  "Education"            = c("school", "universit|college"),
  "Healthcare"           = c("hospital", "pharmac|clinic|doctor|dentist"),
  "Recreation"           = c("\\bparks?\\b", "playground"))

## If a raw column stem prints differently from the paper, map it here,
## e.g. c(fuel = "Fuel station", bus_stop = "Bus stops")
label_override <- character(0)

radius_of <- function(v) suppressWarnings(as.integer(
  ifelse(grepl("_[0-9]+m$", v), sub(".*_([0-9]+)m$", "\\1", v), NA)))
stem_of <- function(v) sub("_[0-9]+m$", "", v)
classify <- function(stem) {
  s <- tolower(gsub("_", " ", stem))
  for (ci in seq_along(cat_keys)) for (ki in seq_along(cat_keys[[ci]]))
    if (grepl(cat_keys[[ci]][ki], s, perl = TRUE)) return(c(ci, ki))
  c(length(cat_keys) + 1, 0)
}
make_label <- function(stem, radius) {
  s <- ifelse(stem %in% names(label_override), label_override[stem], gsub("_", " ", stem))
  s <- paste0(toupper(substr(s, 1, 1)), substring(s, 2))
  ifelse(is.na(radius), s, paste0(s, " within ", radius, "m"))
}

## ---- Build the table: one row per covariate selected by at least one method
all_sel <- sort(unique(unlist(sel_vars)))
cls <- t(vapply(stem_of(all_sel), classify, numeric(2)))
sel_table <- tibble(variable = all_sel, stem = stem_of(all_sel), radius = radius_of(all_sel),
                    cat_id = cls[, 1], key_id = cls[, 2]) |>
  mutate(Category = c(names(cat_keys), "Other")[cat_id],
         Predictor = make_label(stem, radius))
for (m in names(methods)) sel_table[[m]] <- sel_table$variable %in% sel_vars[[m]]
sel_table <- sel_table |>
  arrange(cat_id, key_id, stem, radius) |>
  select(Category, Predictor, variable, all_of(names(methods)))

cat("\nVariable-selection table (Y = selected, . = not selected)\n")
sel_show <- mutate(sel_table, across(all_of(names(methods)), \(x) ifelse(x, "Y", ".")))
print(as.data.frame(sel_show), row.names = FALSE)
cat("\nNumber selected per method (should equal micro_summary$Total_Selected):\n")
print(colSums(sel_table[names(methods)]))
if (any(sel_table$Category == "Other"))
  cat("\nUnclassified covariates (add a keyword to cat_keys):",
      paste(sel_table$variable[sel_table$Category == "Other"], collapse = ", "), "\n")

## ---- LaTeX in the same layout as the manuscript table
selection_latex <- function(tab, file = NULL) {
  ck <- function(x) ifelse(x, "$\\checkmark$", "$\\times$")
  blocks <- lapply(split(tab, factor(tab$Category, levels = unique(tab$Category))), function(g) {
    marks <- apply(as.matrix(as.data.frame(lapply(g[names(methods)], ck))), 1,
                   paste, collapse = " & ")
    c(sprintf("\\multirow{%d}{*}{%s}", nrow(g), g$Category[1]), "",
      paste0("& ", g$Predictor, "\n& ", marks, " \\\\\n"), "\\midrule", "")
  })
  body <- unlist(blocks, use.names = FALSE)
  body <- body[seq_len(length(body) - 2)]            # drop the trailing \midrule
  out <- c("\\begin{table*}[!htbp]", "\\centering",
           paste0("\\caption{Comparison of selected variables. A check mark ($\\checkmark$) ",
                  "indicates that the corresponding predictor was selected by the method, ",
                  "whereas ($\\times$) indicates that it was not selected.}"),
           "\\label{tab:variable_selection}", "",
           "\\renewcommand{\\arraystretch}{1.00}", "\\setlength{\\tabcolsep}{3pt}", "",
           "\\resizebox{0.6\\textwidth}{!}{%", "\\begin{tabular}{llccccc}", "\\toprule",
           paste0("\\textbf{Category} & \\textbf{Predictor} & ",
                  paste0("\\textbf{", names(methods), "}", collapse = " & "), " \\\\"),
           "\\midrule", "", body, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table*}")
  if (!is.null(file)) writeLines(out, file)
  invisible(out)
}

sel_tex_file <- file.path(data_dir, "variable_selection_table_from_code.tex")
sel_tex <- selection_latex(sel_table, sel_tex_file)
cat("\nLaTeX table written to:", sel_tex_file, "\n\n")
cat(sel_tex, sep = "\n")

## ---- Check against the manuscript table
## Save the table from the paper (as-is) to this file to run the comparison.
read_selection_latex <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt <- gsub("(?<!\\\\)%[^\n]*", "", txt, perl = TRUE)            # strip LaTeX comments
  rows <- strsplit(txt, "\\\\", fixed = TRUE)[[1]]                   # split on row ends (\\)
  rows <- gsub("\\\\multirow\\{[^}]*\\}\\{[^}]*\\}\\{[^}]*\\}", "", rows)
  bind_rows(lapply(rows, function(r) {
    f <- trimws(strsplit(r, "&", fixed = TRUE)[[1]])
    if (length(f) < 6) return(NULL)
    mk <- tail(f, 5)
    if (!all(grepl("checkmark|times", mk))) return(NULL)
    tibble(Predictor = f[length(f) - 5],
           !!!setNames(as.list(grepl("checkmark", mk)), names(methods)))
  }))
}

compare_selection <- function(code_tab, paper_tab) {
  key <- \(x) tolower(gsub("\\s+", " ", trimws(x)))
  a <- code_tab  |> transmute(key = key(Predictor), across(all_of(names(methods))))
  b <- paper_tab |> transmute(key = key(Predictor), across(all_of(names(methods))))
  j <- full_join(a, b, by = "key", suffix = c(".code", ".paper")) |>
    mutate(across(-key, \(x) coalesce(x, FALSE)))
  diffs <- map_dfr(names(methods), \(m) {
    cc <- paste0(m, ".code"); pp <- paste0(m, ".paper")
    j |> filter(.data[[cc]] != .data[[pp]]) |>
      transmute(Predictor = key, Method = m, Code = .data[[cc]], Paper = .data[[pp]])
  })
  summ <- tibble(Method = names(methods),
                 Code_selected  = colSums(a[names(methods)]),
                 Paper_selected = colSums(b[names(methods)]),
                 Mismatches     = sapply(names(methods), \(m) sum(diffs$Method == m)))
  only_code <- setdiff(a$key, b$key); only_paper <- setdiff(b$key, a$key)
  
  cat("\n==== Code vs. manuscript table ====\n")
  print(as.data.frame(summ), row.names = FALSE)
  cat("\nRows only in code table :", if (length(only_code))  paste(only_code,  collapse = "; ") else "none", "\n")
  cat("Rows only in paper table:", if (length(only_paper)) paste(only_paper, collapse = "; ") else "none", "\n")
  cat("(A row that appears in both lists under slightly different names is a label issue -- fix via label_override.)\n")
  if (nrow(diffs)) { cat("\nCell-level mismatches:\n"); print(as.data.frame(diffs), row.names = FALSE) }
  else cat("\nAll check marks agree.\n")
  invisible(list(summary = summ, diffs = diffs, only_code = only_code, only_paper = only_paper))
}

paper_tex_file <- file.path(data_dir, "variable_selection_table_paper.tex")
if (file.exists(paper_tex_file)) {
  paper_tab <- read_selection_latex(paper_tex_file)
  cat("\nParsed", nrow(paper_tab), "rows from the manuscript table.\n")
  sel_check <- compare_selection(sel_table, paper_tab)
} else {
  cat("\nTo compare with the manuscript, save its table to:\n  ", paper_tex_file, "\n")
}
## ================================================================
## 11. Figure 8: log(MSE), Hellinger, selection proportion
## ================================================================
measures <- c("log(MSE)", "Hellinger Distance", "Selection Proportion")
fig8_data <- bind_rows(
  plot_table |> transmute(Method, `log(MSE)` = log(Category_MSE), `Hellinger Distance` = Category_Hellinger) |>
    pivot_longer(-Method, names_to = "Measure", values_to = "Plot_Value"),
  micro_summary |> transmute(Method, Measure = "Selection Proportion",
                             Plot_Value = Total_Selected / ncol(Xq))) |>
  mutate(Method = factor(Method, levels = names(methods)), Measure = factor(Measure, levels = measures))

## y-limits: data range padded by 2x its width on each side (data occupy the middle 20%)
pad_limits <- function(x) { r <- range(x, na.rm = TRUE); r + c(-2, 2) * diff(r) }
shape_theme <- theme_bw(base_size = 13) + theme(
  strip.text = element_text(face = "bold", size = 12), axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), legend.position = "bottom",
  legend.title = element_text(face = "bold"), panel.grid = element_blank())

fig8 <- ggplot(fig8_data, aes(Method, Plot_Value, shape = Method)) +
  geom_point(size = 4, stroke = 1.1) +
  facet_wrap(~ Measure, ncol = 3, scales = "free_y") +
  facetted_pos_scales(y = lapply(measures, \(ms) scale_y_continuous(
    limits = pad_limits(fig8_data$Plot_Value[fig8_data$Measure == ms]),
    breaks = scales::pretty_breaks(5), expand = c(0, 0)))) +
  scale_shape_manual(values = shapes) +
  labs(x = NULL, y = NULL, shape = "Method") + shape_theme
print(fig8)

## ================================================================
## 12. Fitted intensity surfaces (penalized fits, dummy points inside city)
## ================================================================
int_df <- data.frame(longitude = U$x, latitude = U$y, intensity_penalized)[yi == 0, ]
int_df <- int_df[lengths(st_intersects(st_as_sf(int_df, coords = c("longitude", "latitude"), crs = 4326),
                                       stl_boundary)) > 0, ]
intensity_long <- int_df |>
  pivot_longer(-c(longitude, latitude), names_to = "Method", values_to = "Intensity") |>
  mutate(Method = factor(Method, levels = names(methods)))

rng <- \(v) range(v, na.rm = TRUE)
intensity_plot <- ggplot(intensity_long, aes(longitude, latitude, fill = log1p(Intensity))) +
  geom_tile() +
  geom_sf(data = stl_boundary, fill = NA, linewidth = 0.45, inherit.aes = FALSE) +
  geom_point(data = df, aes(longitude, latitude), shape = 3, size = 1, stroke = 0.5, inherit.aes = FALSE) +
  facet_wrap(~ Method, ncol = 3) +
  scale_fill_viridis_c(option = "C", name = expression(log(1 + hat(pi)(s)))) +
  coord_sf(expand = FALSE, xlim = rng(intensity_long$longitude), ylim = rng(intensity_long$latitude)) +
  scale_x_continuous(breaks = pretty(rng(intensity_long$longitude), n = 3)) +
  scale_y_continuous(breaks = pretty(rng(intensity_long$latitude), n = 4)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 9), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), panel.spacing = unit(0.8, "lines"))
print(intensity_plot)

## ================================================================
## 13. Two occurrence maps with observed features
## ================================================================
display_features <- c("primary_road_500m", "residential_road_500m", "service_road_500m",
                      "buildings_50m", "parking_500m")
fld <- function(pt, v, alt = "NA") if (v %in% names(pt)) as.character(pt[[v]]) else alt

occurrence_map <- function(i0, marker_color = "blue", tag = "(a)", features = display_features) {
  pt <- df[i0, , drop = FALSE]
  if (!is.finite(pt$longitude) || !is.finite(pt$latitude)) stop("Selected occurrence has invalid coordinates.")
  features <- intersect(features, names(pt))
  desc  <- paste0(fld(pt, "Type", "Unknown event"), " at ", fld(pt, "Road_Name", "Unknown road"))
  popup <- paste0("<b>Type:</b> ", fld(pt, "Type"), "<br><b>Road:</b> ", fld(pt, "Road_Name"),
                  "<br><b>Time:</b> ", fld(pt, "Occurance_Time"))
  map <- leaflet(width = "100%", height = 460) |> addProviderTiles(providers$CartoDB.Positron) |>
    add_city(stl_boundary) |>
    addAwesomeMarkers(pt$longitude, pt$latitude, popup = popup, label = desc,
                      icon = makeAwesomeIcon(icon = "star", markerColor = marker_color, iconColor = "white", library = "fa")) |>
    fit_city(stl_boundary)
  feat <- paste0("<p><b>Observed features:</b></p><ul style='font-size:14px;'>",
                 paste0("<li><b>", features, "</b>: ", unlist(pt[1, features]), "</li>", collapse = ""), "</ul>")
  tags$div(style = "border:1px solid #bdbdbd;border-radius:6px;padding:8px;",
           tags$h3(paste0(tag, " Occurrence ", i0, ": ", desc), style = "text-align:center;font-family:Arial;font-size:16px;"),
           map,
           tags$div(HTML(feat), style = "background:#f7f7f7;padding:10px;margin-top:8px;border-radius:6px;border:1px solid #ddd;"))
}

two_map_figure <- browsable(tags$div(
  style = "width:100%;display:grid;grid-template-columns:1fr 1fr;gap:20px;align-items:start;",
  occurrence_map(1, "blue", "(a)"), occurrence_map(10, "purple", "(b)")))
two_map_figure

## ================================================================
## 14. Four representative locations: locally selected predictors
##     (uses the TIGER/Line St. Louis City boundary, as in the original)
## ================================================================
## ================================================================
## 14. Four representative locations: one image per location
##     (uses the TIGER/Line St. Louis City boundary, as in the original)
##     Needs: webshot2 (install.packages("webshot2")); it uses a local Chrome/Chromium.
## ================================================================
stl_city <- counties(state = "MO", cb = TRUE, class = "sf") |> filter(GEOID == "29510") |> st_transform(4326)

choose_four_city_locations <- function(d) {
  lon_mid <- median(d$longitude, na.rm = TRUE); lat_mid <- median(d$latitude, na.rm = TRUE)
  d |> mutate(city_part = case_when(
    longitude <= lon_mid & latitude >= lat_mid ~ "Nearby University City",
    longitude >  lon_mid & latitude >= lat_mid ~ "Northeast St. Louis",
    longitude <= lon_mid & latitude <  lat_mid ~ "South St. Louis",
    TRUE                                       ~ "East St. Louis")) |>
    group_by(city_part) |>
    slice_min((longitude - median(longitude))^2 + (latitude - median(latitude))^2, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(city_part, Reference, Road_Name, Type, Occurance_Time, longitude, latitude)
}
selected_at <- function(lon, lat, sel) {
  v <- sel |> filter(x_lo <= lon, x_hi >= lon, y_lo <= lat, y_hi >= lat) |>
    distinct(covariate) |> arrange(covariate) |> pull(covariate)
  if (length(v)) v else "No selected predictors"
}
pretty_name <- function(x) tools::toTitleCase(gsub("([0-9]+)m$", "(\\1m)", gsub("_", " ", x)))

demo_locations <- choose_four_city_locations(df) |>
  mutate(lasso_selected = map2(longitude, latitude, \(x, y) selected_at(x, y, sel_local$LLI)),
         scad_selected  = map2(longitude, latitude, \(x, y) selected_at(x, y, sel_local$LSI)))
print(demo_locations)

## One stand-alone panel: map on the left, information box on the right
location_panel <- function(row, tag = NULL) {
  preds <- function(v) {
    if (identical(v, "No selected predictors")) return(tags$div(style = "font-size:18px;line-height:1.35;", v))
    tags$div(style = "display:grid;grid-template-columns:1fr 1fr;column-gap:24px;row-gap:6px;font-size:17px;line-height:1.25;",
             lapply(pretty_name(unique(v)), tags$div))
  }
  section <- function(label, col, v) tagList(
    tags$div(style = paste0("color:", col, ";font-weight:bold;font-size:21px;margin-bottom:8px;"), label), preds(v))
  title <- if (is.null(tag)) row$city_part else paste(tag, row$city_part)
  box <- tags$div(
    style = "width:560px;min-height:520px;padding:18px;border:2px solid #333;border-radius:8px;
             background-color:#f8f9fa;font-family:Arial;box-sizing:border-box;",
    tags$h3(title, style = "color:#2c3e50;font-weight:bold;margin:0 0 12px 0;font-size:27px;"),
    tags$div(style = "font-size:19px;line-height:1.35;margin-bottom:10px;",
             tags$b("Road: "), row$Road_Name, tags$br(), tags$b("Crime: "), row$Type, tags$br(),
             tags$b("Time: "), as.character(row$Occurance_Time), tags$br(),
             tags$b("Longitude: "), round(row$longitude, 5), tags$br(), tags$b("Latitude: "), round(row$latitude, 5)),
    tags$hr(style = "margin:12px 0;border-top:1px solid #888;"),
    section("Local LASSO selected predictors", "#008000", row$lasso_selected[[1]]),
    tags$hr(style = "margin:12px 0;border-top:1px solid #888;"),
    section("Local SCAD selected predictors",  "#0055cc", row$scad_selected[[1]]))
  m <- leaflet(width = "600px", height = "520px", options = leafletOptions(zoomControl = FALSE)) |>
    addProviderTiles(providers$CartoDB.Positron) |>
    add_city(stl_city, 4) |> add_plus(row$longitude, row$latitude, 44) |> fit_city(stl_city)
  tags$div(id = "panel",
           style = "display:flex;align-items:flex-start;gap:22px;padding:15px;width:1200px;
                    min-height:550px;background-color:white;box-sizing:border-box;",
           m, box)
}

## Build the four images (shown in the viewer) ...
tags4 <- c("(a)", "(b)", "(c)", "(d)")
location_plots <- lapply(seq_len(nrow(demo_locations)), \(i)
                         browsable(location_panel(demo_locations[i, ], tag = tags4[i]))) |>
  setNames(demo_locations$city_part)

## Four separate figures (no files are saved): print figure1, ..., figure4
tags4 <- c("(a)", "(b)", "(c)", "(d)")
figure1 <- browsable(location_panel(demo_locations[1, ], tag = tags4[1]))
figure2 <- browsable(location_panel(demo_locations[2, ], tag = tags4[2]))
figure3 <- browsable(location_panel(demo_locations[3, ], tag = tags4[3]))
figure4 <- browsable(location_panel(demo_locations[4, ], tag = tags4[4]))
cat(sprintf("figure%d: %s\n", 1:4, demo_locations$city_part), sep = "")

figure1
figure2
figure3
figure4

## ================================================================
## 15. Figure S3: LLI regularization path
## ================================================================
lf    <- fits$LLI
Bpath <- as.matrix(lf$fit$beta); lam <- lf$fit$lambda; log_l <- log(lf$lambda)
idx   <- which.min(abs(lam - lf$lambda))
act_pct <- 100 * colSums(abs(Bpath) > tol0) / nrow(Bpath)
thm   <- theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
vline <- geom_vline(xintercept = log_l, linetype = "dashed", linewidth = 0.8, color = "red")
xmark <- function(x, y) annotate("point", x = x, y = y, shape = 4, color = "red", size = 5, stroke = 1.5)
lab_x <- labs(x = expression(log(lambda)))

p1 <- ggplot(lf$cbic_path, aes(log(lambda), CBIC)) + geom_line(linewidth = 0.9, color = "#0072B2") +
  xmark(log_l, lf$cbic) + lab_x + thm
p2 <- as.data.frame(t(Bpath[apply(abs(Bpath), 1, max) > tol0, , drop = FALSE])) |>
  mutate(lambda = lam) |> pivot_longer(-lambda, names_to = "Coef", values_to = "delta") |>
  ggplot(aes(log(lambda), delta, group = Coef)) +
  geom_line(linewidth = 0.35, alpha = 0.5, color = "#D55E00") + geom_hline(yintercept = 0, linewidth = 0.3) +
  vline + lab_x + labs(y = "Local coefficients (scaling and detail)") + thm
p3 <- ggplot(data.frame(log_lam = log(lam), act = act_pct), aes(log_lam, act)) +
  geom_line(linewidth = 0.9, color = "#009E73") + vline + xmark(log(lam[idx]), act_pct[idx]) +
  lab_x + labs(y = "Active coefficients (%)") + thm
print(p1 | p2 | p3)

K0 <- sum(abs(Bpath[, idx]) > tol0); P <- nrow(Bpath)
cat(sprintf(paste0("\nLLI at CBIC-selected lambda\n---------------------------------\n",
                   "Selected lambda     = %g\nlog(lambda)         = %g\nTotal coefficients  = %d\nActive coefficients = %d\n",
                   "Zero coefficients   = %d\nActive percentage   = %.2f%%\nSparsity percentage = %.2f%%\n"),
            lf$lambda, log_l, P, K0, P - K0, act_pct[idx], 100 - act_pct[idx]))
print(data.frame(Method = "LLI", Lambda = lf$lambda, Log_Lambda = log_l, CBIC = lf$cbic))