test_that("Haar basis has the expected atoms and values", {
  set.seed(1)
  s <- cbind(runif(200), runif(200))
  h <- haar_basis_2d(s, J = 3)
  expect_equal(ncol(h$B), 1 + 3 * (1 + 4 + 16))
  expect_equal(nrow(h$meta), ncol(h$B))
  allowed <- c(0, 1, 2, 4)  # |phi * psi| = 2^j at level j
  expect_true(all(vapply(abs(as.vector(h$B)), function(v) min(abs(v - allowed)) < 1e-12, logical(1))))

  f <- finest_level(h)
  expect_equal(ncol(f$B), 3 * 16)
  expect_true(all(f$meta$j == 2))
  expect_equal(f$meta$atom_index, seq_len(48))
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
