# pmxTools (development version)

* Added `sample_nhanes()`, which simulates virtual pediatric populations (ages 2-17) with realistic body weight and/or height by MEC-weighted resampling of children from four NHANES releases (2013-2023), bundled as `nhanes_peds`. Resampling keeps each child's weight and height together; `method = "smooth"` adds correlated kernel noise for continuous values. `compare_nhanes()` and `plot_nhanes()` check the simulated population against the reference. A new vignette, `vignette("nhanes", package = "pmxTools")`, walks through them.

* Added `het_cor()`, which computes a correlation matrix for mixed continuous and ordinal variables (Pearson, polyserial or polychoric for each pair) and can draw it with `plot_cor_triangle()`. The pairwise estimators are available as `polychoric()` and `polyserial()`.

## Breaking changes

* `patchwork` has moved from `Depends` to `Suggests`, so `library(pmxTools)` no longer attaches it. Call `library(patchwork)` to combine plots with `+`, `/` and `plot_layout()`.

* The PK calculation functions now stop with an error for inputs outside their domain instead of returning `NaN`, `Inf` or wrong values: non-positive clearances, volumes, `ka`, `dur`, `tinf`, `tau` or `step`; negative `tlag`, `dose` or times; a `ka` equal to a disposition rate constant (where the closed-form solution is singular); and, for the steady-state functions, an infusion or zero-order absorption (`tinf`, `dur`) longer than `tau`, which the steady-state solutions do not cover.

* `calc_derived()` and the `calc_derived_*()` functions now stop when vector inputs have different lengths (other than 1) rather than recycling them.

* `read_nmcov()` now labels the rows and columns of the variance-covariance matrix with the parameter names from the `.cov` file (`THETA1`, `SIGMA(1,1)`, `OMEGA(2,1)`, ...). Rows were previously numbered and column names were altered (`SIGMA.1.1.`). As a result, the columns of `sample_uncert()` output are named after the parameters instead of `V1`, `V2`, ...

## Bug fixes

* `pk_curve()` now superposes additional doses by the time elapsed since each dose, so irregular time grids, grids not starting at zero and `addl` doses beyond the observed window give correct concentrations. Previously doses were added by shifting rows, which was only correct for regular grids aligned with `ii`.

* `pk_curve(ss = TRUE)` now gives a steady-state curve (steady state reached before the dose at time 0, as for `SS=1` in NONMEM, followed by any `addl` doses and then washout). The argument was previously ignored. `pk_curve()` also validates `model`, `t`, `dose`, `ii`, `addl` and `ss`.

* `calc_derived()` and the `calc_derived_*()` functions calculate `Cmax`, `Tmax`, `Cmaxss`, `AUCtau` and their dose-normalized forms separately for each parameter set when given vector inputs, returning one value per set. Previously one value was returned from the pooled profiles.

* `AUCtau` is now integrated over the whole dosing interval when `step` does not divide `tau` exactly.

* `calc_sd_1cmt()`, `calc_sd_2cmt()`, `calc_sd_3cmt()`, `calc_ss_1cmt()`, `calc_ss_2cmt()` and `calc_ss_3cmt()` are now implemented. They choose the regimen-specific function from the parameters given (`tinf`, `dur`, `ka`, `tlag`); previously they stopped with "Not yet implemented".

* `read_nmtables()` merges per-subject (`FIRSTONLY`) tables onto per-record tables correctly when there are three or more table files, by `ID`, keeping all per-record rows in their original order. It warns about subjects missing from the `FIRSTONLY` tables and stops if there is no shared `ID` column. Columns repeated across tables (such as `ID` and `TIME`) are now kept once.

* `read_nmtables()` reads tables with exactly three columns and NSIM > 1 output correctly; both previously lost the header row. `read_nmtables()` and `read_nm_multi_table()` no longer turn a one-row table into a one-column table.

* `read_nm()` now parses the XML document instead of requiring its last line to be `</nm:output>`, so files with trailing blank lines are read. Unparseable files still give a warning and `NULL`.

* `read_scm()` and `plot_scm()` no longer pick up objects with the same names from the global environment, and `rnm()` no longer deletes `THETA`, `ETAD` or `EPSD` objects from the global environment.

* `read_scm()` records the model directory of each step for any path style; previously only Windows paths with an upper-case drive letter were recognized and others were `NA`.

* `plot_scm()` now draws SCMs that stopped after a single forward step.

* `calc_derived()` and the `calc_derived_*()` functions now accept several output names in `type`, as documented, returning a list of those outputs.

* `read_nmtables(quiet = TRUE)` no longer reports table files that were not found.

* `ggplot2` is now imported selectively instead of in full.

# pmxTools 1.6

* Added `cut_quantile()`, which creates quantile-based bins for one or more continuous variables, with per-variable control over the number of groups.

* Added `datamap()`, which plots individual dosing and observation records over time.

* `dgr_table()` now takes `fields` as a named character vector mapping column names to display labels. The separate `names` argument is deprecated, but unnamed `fields` with `names` continues to work.

* Added dependencies `purrr` and `tidyr`.

* Fixed `read_nm()`, which passed the file extension to `xml2::read_xml()` as an encoding rather than appending it to the file name, and read the unresolved file name instead of the one located on disk. The `directory` argument and file names given without the `.xml` extension were therefore ignored, and macOS reported an `Unsupported encoding: .xml` warning (#36).

* Removed the blanket `stats` and `utils` imports, which caused a `replacing previous import 'stats::filter' by 'dplyr::filter'` warning when the package was loaded (#36).

* `calc_derived()` and the `calc_derived_*()` functions now return `AUCtau` and `AUCtau_dose_normalized` as plain numbers.  PKNCA 0.12.1.9000 attaches a `method` attribute to the AUC it returns, which propagated into the derived-parameter list and broke comparisons against plain numeric values (#37).

* Test snapshots are now committed, the CI workflows have been refreshed, and the documentation has been regenerated with roxygen2 8.1.0 (#36).

* `calc_derived()` and the `calc_derived_*()` functions now work with vector inputs (#29).

# pmxTools 1.5

* `plot_dist` has been rewritten to adapt to `ggplot2` 4.0.0 changes, and to make it more intuitive to interpret

* Dependency `gghalves` has been removed

* Many fixes and updates to documentation

# pmxTools 1.4

* `dgr_table()` now provides geometric means by default

* `gm()` now provides options to strip `NA` values (`na.rm`) and non-positive values (`neg.rm`) before computation.

# pmxTools 1.3

* Added NCA parameter estimation to `calc_derived_1cpt()`, `calc_derived_2cpt()` and `calc_derived_3cpt()`, if dose and other required information (e.g. `tinf`, `dur`, `tau`) is provided.

* Added a warning about flip-flop kinetics and their potential effects on derived half-lives to `calc_derived_1cpt()`, `calc_derived_2cpt()`, and `calc_derived_3cpt()`.

* Added `blq_trans()` and `blq_log_trans()` for assistance with visualization of measurements below a limit of quantification.

* Added a test for valid XML to `read_nm()`.

* Small fixes and enhancements to `read_nmtable_single()` and `read_nm_multi_table()` to ensure that stray text is properly handled while reading in NONMEM output files.

* `count_na()` now throws a warning if `NaN` values are included in the `NA` count, and indicates how many of them there are.

* Changed the significant digits functionality for 95% confidence intervals for `get_theta()`, `get_omega()` and `get_sigma()`, and all functions that apply these.

* Replaced `gridExtra` with `patchwork`.

* Numerous small fixes.

* Removed some unnecessary dependencies.

* Added many new unit tests, increasing test coverage considerably.

# pmxTools 1.2.4

* Fixed some minor documentation issues. Thanks to Julien Grassot for spotting these.

# pmxTools 1.2.3

* Fixed another error in `sample_uncert()` which was still crashing the function.

* Added helper functions `count_na()`, `dgr_table()`, `gcv()`, `fmt_signif()`.

* Added distribution plotting function `plot_dist()`. 

# pmxTools 1.2.2

* Fixed an error in `sample_uncert()` which was crashing the function.

# pmxTools 1.2.1

* Added vignette describing PK curves.

* Fixed a systematic error in `calc_sd_1cmt_linear_oral_0_lag()`, `calc_ss_1cmt_linear_oral_0_lag()`, `calc_sd_2cmt_linear_oral_0_lag()` which resulted in incorrect curves being plotted.

# pmxTools 1.2

* Rewrote `plot_scm` function to produce tree diagrams via `DiagrammeR`.

* Fixed a rare bug in `get_auc` in which measurements from different individuals could be erroneously mixed.

* Amended (or added) `read_nm()`, `read_nm_all()`, `read_nmext()`, `read_nmcov()`, `read_nmtables()` functions to allow reading of multiple NONMEM estimation steps.

* Clarified documentation for `calc_derived()`.

* Now using `xml2` to read NONMEM-generated XML files.

* Fixed partial argument match warnings encountered during tests.

* Add function for calculting geometric CV.
