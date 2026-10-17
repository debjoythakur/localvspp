# localvspp 0.2.0

Brings the package in line with the final analysis script for the revised
manuscript.

* **Localized basis changed.** `fit_localvs()` now uses the new
  `local_haar_basis()`: one scaling atom plus three detail atoms per block at
  level `J - 1`, i.e. `4^J` atoms (64 for `J = 3`), as in Eq. (5) of the paper.
  Version 0.1.0 used only the 48 detail atoms (`finest_level()`, now
  deprecated). Results from 0.1.0 and 0.2.0 therefore differ.
* Points on the upper edge of the study window now fall in the last Haar
  block; before, they fell in no block and got an all-zero row.
* `bt_quadrature()` and `fit_localvs()` default to equal cell-area weights
  (`"cell"`), as in the paper's analysis. Use `"spatstat"` for the standard
  Berman–Turner weights.
* New `selection_table()`, `selection_latex()`, `read_selection_latex()` and
  `compare_selection()` build the variable-selection table (Table S1), write
  it as LaTeX and check it against the manuscript.
* New `occurrence_panel()` draws the event maps of Figure 1.
* `plot_cbic_path()` labels the coefficient-path panel as in Figure S3.
* `inst/scripts/real-data/final_analysis.R` is the script behind the revised
  manuscript.

# localvspp 0.1.0

* First release.
