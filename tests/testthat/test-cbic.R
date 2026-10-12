test_that("CBIC-tuned LASSO and SCAD pick the CBIC minimum", {
  set.seed(42)
  n <- 400; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  off <- rep(log(0.5), n)
  y <- rbinom(n, 1, plogis(-1 + X[, 1] - X[, 2]))

  fl <- fit_lasso_cbic(X, y, off)
  expect_equal(fl$cbic, min(fl$cbic_path$CBIC))
  expect_length(fl$beta, p)

  fs <- fit_scad_cbic(X, y, off, nlambda = 10, mc.cores = 1)
  expect_equal(fs$cbic, min(fs$cbic_path$CBIC))
  expect_equal(dim(fs$beta_path), c(p, 10))

  fa <- fit_adaptive_lasso_cbic(X, y, off)
  expect_length(fa$beta, p)
})

test_that("quadrature marks data and dummy points", {
  set.seed(3)
  xy <- cbind(runif(60), runif(60))
  q <- bt_quadrature(xy, ntile = c(8, 8))
  expect_equal(sum(q$y), 60)
  expect_equal(length(q$w), nrow(q$U))
  expect_true(all(q$w > 0))
  qc <- bt_quadrature(xy, ntile = c(8, 8), weights = "cell")
  expect_true(all(qc$w == qc$cell_area))
})
