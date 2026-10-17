#' localvspp: Multi-Resolution Local Variable Selection for Spatial Point Processes
#'
#' The main entry point is [fit_localvs()], which runs the full pipeline:
#' Berman-Turner quadrature ([bt_quadrature()]), covariate smoothing
#' ([smooth_covariates()]), the Haar-localized design ([local_haar_basis()],
#' [localized_design()]) and five CBIC-tuned penalized fits:
#'
#' * `GLI`  global LASSO
#' * `GSI`  global SCAD
#' * `GALI` global adaptive LASSO
#' * `LLI`  local (Haar-localized) LASSO
#' * `LSI`  local (Haar-localized) SCAD
#'
#' Use [evaluate_localvs()] for in-sample and post-selection comparisons,
#' [plot_performance()], [plot_intensity()], [plot_cbic_path()],
#' [occurrence_panel()] and [location_panel()] for figures, and
#' [selection_table()] / [selection_latex()] for the variable-selection table.
#' Example data: [stl_crime] and [stl_boundary].
#'
#' @keywords internal
#' @importFrom utils globalVariables
"_PACKAGE"

## Method codes used throughout the package
method_labels <- c(GLI = "Global LASSO", GSI = "Global SCAD", GALI = "Adaptive LASSO",
                   LLI = "Local LASSO",  LSI = "Local SCAD")
local_methods <- c("LLI", "LSI")
method_shapes <- c(GLI = 1, GSI = 2, GALI = 3, LLI = 4, LSI = 5)

utils::globalVariables(c("Method", "Plot_Value", "longitude", "latitude", "Intensity",
                         "lambda", "CBIC", "delta", "Coef", "log_lambda", "active"))