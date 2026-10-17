fake_fit <- function() {
  sel <- list(
    GLI = data.frame(covariate = c("bars_50m", "primary_road_200m"), coefficient = c(0.1, -0.2)),
    LLI = data.frame(covariate = c("bars_50m", "bars_50m", "parks_1000m"), coefficient = c(0.3, 0.1, 0.2)))
  structure(list(selection = sel), class = "localvs")
}

test_that("selection table groups, orders and marks covariates", {
  tab <- selection_table(fake_fit())
  expect_equal(tab$variable, c("primary_road_200m", "bars_50m", "parks_1000m"))
  expect_equal(tab$Category, c("Road Types", "Buildings / Land Use", "Recreation"))
  expect_equal(tab$Predictor[1], "Primary road within 200m")
  expect_equal(tab$GLI, c(TRUE, TRUE, FALSE))
  expect_equal(tab$LLI, c(FALSE, TRUE, TRUE))
})

test_that("LaTeX table round-trips through read_selection_latex", {
  tab <- selection_table(fake_fit())
  f <- tempfile(fileext = ".tex")
  tex <- selection_latex(tab, f)
  expect_true(any(grepl("\\\\begin\\{tabular\\}\\{llcc\\}", tex)))
  back <- read_selection_latex(f, methods = c("GLI", "LLI"))
  expect_equal(back$Predictor, tab$Predictor)
  cmp <- compare_selection(tab, back, verbose = FALSE)
  expect_equal(nrow(cmp$diffs), 0)
  back$LLI[1] <- TRUE
  cmp <- compare_selection(tab, back, verbose = FALSE)
  expect_equal(nrow(cmp$diffs), 1)
})
