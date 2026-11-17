test_that("built-in data are available", {
  skip_if_not(exists("stl_crime"), "stl_crime has not been built yet (run data-raw/stl_crime.R)")
  expect_true(all(c("longitude", "latitude") %in% names(stl_crime)))
  expect_s3_class(stl_boundary, "sf")
  expect_equal(nrow(stl_crime), 405)
  expect_equal(ncol(crime_covariates(stl_crime)), 140)
  expect_true(mean(inside_boundary(stl_crime, stl_boundary)) > 0.5)
})
