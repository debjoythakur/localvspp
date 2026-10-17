# localvspp: Multi-Resolution Variable Selection for Spatial Point Processes

`localvspp` is an R package for **localized variable selection** in the
intensity of spatial point processes. Using Haar wavelets and penalized
likelihood, it identifies **which predictors matter, where, and at what
spatial scale**. The motivating application, crime occurrences in
St. Louis, Missouri, ships with the package as example data.

## Installation

```r
# install.packages("remotes")
remotes::install_github("debjoythakur/localvspp")
```

## Quick start

```r
library(localvspp)

dim(stl_crime)               # built-in crime events + built-environment features
plot(sf::st_geometry(stl_boundary))

set.seed(2026)               # the adaptive-LASSO weights use random CV folds
fit <- fit_localvs(stl_crime)  # defaults reproduce the paper (J = 3, cell weights)
fit                          # lambda, CBIC and number of selected terms per method

selected_covariates(fit)$LLI                       # covariates active anywhere (local LASSO)
selected_at(fit, lon = -90.25, lat = 38.63, "LLI") # covariates active at a location

ev <- evaluate_localvs(fit)
ev$performance                                     # post-selection MSE and Hellinger distance
plot_performance(fit, ev)                          # log(MSE), Hellinger, selection proportion
plot_intensity(fit, ev, boundary = stl_boundary)   # fitted intensity surfaces
plot_cbic_path(fit, "LLI")                         # CBIC path, coefficient paths, sparsity
plot_ripleys_k(stl_crime, stl_boundary)            # clustering diagnostic
microscale_summary(fit)                            # 50/100/200 m features selected

tab <- selection_table(fit)                        # Table S1: selected variables per method
selection_latex(tab, "variable_selection.tex")     # ... as LaTeX
# compare_selection(tab, read_selection_latex("table_from_paper.tex"))

loc <- representative_locations(stl_crime)         # Figures 6-9: local selections at 4 events
htmltools::browsable(location_panel(fit, loc[1, ], stl_boundary, large = TRUE))
htmltools::browsable(occurrence_panel(stl_crime, 1, stl_boundary))  # Figure 1
```

The fit takes a while on the full data: the localized design has one column
per covariate × Haar atom, and SCAD runs one LLA fit per lambda (in parallel
on macOS/Linux; set `scad_args = list(mc.cores = 1)` to disable).

## Data

* `stl_crime`: 405 crime events at unique locations (St. Louis Police,
  1 March to 30 April 2026) with 155 built-environment covariates: 31
  OpenStreetMap feature types counted within 50, 100, 200, 500 and 1000 m.
  `crime_covariates(stl_crime)` returns the 140 non-constant candidate
  predictors.
* `stl_boundary`: St. Louis city-limit polygon (`sf`, EPSG:4326).

`data-raw/stl_crime.R` rebuilds both from the raw files in `data-raw/`.

## Methods

| Code | Method | Design |
|---|---|---|
| `GLI`  | Global LASSO          | smoothed covariates |
| `GSI`  | Global SCAD (LLA)     | smoothed covariates |
| `GALI` | Global adaptive LASSO | smoothed covariates |
| `LLI`  | Local LASSO           | covariates × Haar atoms at level J − 1 |
| `LSI`  | Local SCAD (LLA)      | covariates × Haar atoms at level J − 1 |

Pipeline (`fit_localvs()`):

1. **Berman–Turner quadrature** (`bt_quadrature()`): grid dummy points turn the
   Poisson point-process likelihood into a weighted Poisson regression with
   offset `log(w)`.
2. **Covariate smoothing** (`smooth_covariates()`): event-level covariates are
   kernel-smoothed to the quadrature points and standardized.
3. **Haar localization** (`local_haar_basis()`, `localized_design()`): each
   covariate is multiplied by every 2-D Haar atom at level `J - 1`, one
   scaling and three detail atoms per block (4^J atoms; 64 for J = 3), giving
   coefficients that are constant on the 2^J × 2^J sub-blocks.
4. **Penalized estimation** (`fit_lasso_cbic()`, `fit_adaptive_lasso_cbic()`,
   `fit_scad_cbic()`): the tuning parameter minimizes the composite BIC
   `CBIC = -2 ℓ_PP + K0 log N(R)`.
5. **Post-selection evaluation** (`evaluate_localvs()`): unpenalized refit on
   the selected design, compared through grid-cell count distributions.

`bt_quadrature()` uses equal cell-area weights by default, which is what the
paper's analysis used. Pass `quad_weights = "spatstat"` to `fit_localvs()`
for the standard Berman–Turner quadrature weights.

## Abstract

Lasso, adaptive Lasso and SCAD are standard approaches to variable selection
with a large number of predictors, and many researchers have used these
penalized methods to estimate the intensity function of spatial point
processes with a diverging number of predictors. We take a multi-resolution
perspective on variable selection for spatial point process data. Its
advantage is twofold: it selects the predictors efficiently, and it shows
which points are responsible for selecting a predictor at a specific
resolution. The research is motivated by crime and accident occurrences in
St. Louis and its neighborhoods, where it is more relevant to select
predictors at the local level: we learn which set of predictors is relevant
to crimes or accidents in which parts of St. Louis. Theoretical and
simulation results justify the accuracy of local-level variable selection
during intensity function estimation.

## Repository structure

```
├── DESCRIPTION, NAMESPACE, LICENSE   package metadata
├── R/                                package functions
├── data/                             built-in datasets (stl_crime, stl_boundary)
├── man/                              help pages (generated by roxygen2)
├── tests/testthat/                   unit tests
├── inst/scripts/                     original analysis and simulation scripts
│   ├── real-data/                    final_analysis.R = script behind the revised paper
│   └── simulation/
├── data-raw/                         raw inputs and the script that builds data/
│   └── gis/                          QGIS project, GeoPackage, bounding rectangle
└── Results/                          manuscript results (not part of the package)
```

The original scripts are installed with the package:
`system.file("scripts", package = "localvspp")`.

## Development

```r
source("data-raw/stl_crime.R")   # rebuild data/ after changing the raw files
devtools::document()             # regenerate man/ and NAMESPACE
devtools::test()
devtools::check()
```

## Citation

Thakur, D. *Multi-resolution variable selection for spatial point processes.*
(Manuscript in preparation.)
