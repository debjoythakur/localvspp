test_that("multi-resolution Haar basis has the expected atoms and values", {
  set.seed(1)
  s <- cbind(runif(200), runif(200))
  h <- haar_basis_2d(s, J = 3)
  expect_equal(ncol(h$B), 1 + 3 * (1 + 4 + 16))
  expect_equal(nrow(h$meta), ncol(h$B))
  allowed <- c(0, 1, 2, 4)  # |phi * psi| = 2^j at level j
  expect_true(all(vapply(abs(as.vector(h$B)), function(v) min(abs(v - allowed)) < 1e-12, logical(1))))

  f <- suppressWarnings(finest_level(h))
  expect_warning(finest_level(h), "deprecated")
  expect_equal(ncol(f$B), 3 * 16)
  expect_true(all(f$meta$j == 2))
})

test_that("local Haar basis has 4^J atoms and every point lies in one block", {
  set.seed(2)
  s <- cbind(runif(300), runif(300))
  for (J in 1:3) {
    h <- local_haar_basis(s, J = J)
    expect_equal(ncol(h$B), 4^J)
    expect_equal(h$meta$atom_index, seq_len(4^J))
    expect_true(all(h$meta$j == J - 1))
    expect_equal(as.vector(table(h$meta$type)), rep(4^(J - 1), 4))
    phi <- h$B[, h$meta$type == "phi", drop = FALSE]
    expect_true(all(rowSums(phi != 0) == 1))   # one scaling atom per point
  }
})

test_that("points on the upper edge fall in the last block", {
  s <- rbind(c(0, 0), c(1, 1), c(0.5, 0.2))
  h <- local_haar_basis(s, J = 2)
  phi <- h$B[, h$meta$type == "phi"]
  expect_true(all(rowSums(phi != 0) == 1))
  top <- h$meta[h$meta$type == "phi", ][phi[2, ] != 0, ]
  expect_equal(c(top$kx, top$ky), c(1, 1))
})

test_that("local Haar atoms span block-constant surfaces", {
  set.seed(3)
  s <- rbind(c(0, 0), c(1, 1), cbind(runif(500), runif(500)))  # range exactly [0, 1]
  h <- local_haar_basis(s, J = 2)
  ## a surface constant on the 4 x 4 sub-blocks is reproduced exactly
  cell <- (floor(4 * pmin(s[, 1], 1 - 1e-12)) * 4 + floor(4 * pmin(s[, 2], 1 - 1e-12)))
  f <- c(rnorm(16))[cell + 1]
  fit <- stats::lm.fit(h$B, f)
  expect_lt(max(abs(fit$residuals)), 1e-8)
})

test_that("localized design is covariate-major", {
  B <- matrix(1:6, 3, 2)
  X <- matrix(c(1, 1, 1, 2, 2, 2), 3, 2, dimnames = list(NULL, c("a", "b")))
  meta <- data.frame(type = "psi_x", j = 1, kx = 0:1, ky = 0)
  Z <- localized_design(X, B, meta)
  expect_equal(dim(Z), c(3, 4))
  expect_equal(Z[, 3], 2 * B[, 1])
  expect_match(colnames(Z)[1], "^a__psi_x_j1_kx0_ky0$")
})
