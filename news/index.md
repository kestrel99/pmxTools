# Changelog

## pmxTools (development version)

- Added
  [`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md),
  which simulates virtual pediatric populations (ages 2-17) with
  realistic body weight and/or height by MEC-weighted resampling of
  children from four NHANES releases (2013-2023), bundled as
  `nhanes_peds`. Resampling keeps each child’s weight and height
  together; `method = "smooth"` adds correlated kernel noise for
  continuous values.
  [`compare_nhanes()`](https://kestrel99.github.io/pmxTools/reference/compare_nhanes.md)
  and
  [`plot_nhanes()`](https://kestrel99.github.io/pmxTools/reference/plot_nhanes.md)
  check the simulated population against the reference. A new vignette,
  [`vignette("nhanes", package = "pmxTools")`](https://kestrel99.github.io/pmxTools/articles/nhanes.md),
  walks through them.

- Added
  [`het_cor()`](https://kestrel99.github.io/pmxTools/reference/het_cor.md),
  which computes a correlation matrix for mixed continuous and ordinal
  variables (Pearson, polyserial or polychoric for each pair) and can
  draw it with
  [`plot_cor_triangle()`](https://kestrel99.github.io/pmxTools/reference/plot_cor_triangle.md).
  The pairwise estimators are available as
  [`polychoric()`](https://kestrel99.github.io/pmxTools/reference/polychoric.md)
  and
  [`polyserial()`](https://kestrel99.github.io/pmxTools/reference/polyserial.md).

### Breaking changes

- `patchwork` has moved from `Depends` to `Suggests`, so
  [`library(pmxTools)`](https://github.com/kestrel99/pmxTools) no longer
  attaches it. Call
  [`library(patchwork)`](https://patchwork.data-imaginist.com) to
  combine plots with `+`, `/` and
  [`plot_layout()`](https://patchwork.data-imaginist.com/reference/plot_layout.html).

- The PK calculation functions now stop with an error for inputs outside
  their domain instead of returning `NaN`, `Inf` or wrong values:
  missing or infinite values; non-positive clearances, volumes, `ka`,
  `dur`, `tinf`, `tau` or `step`; negative `tlag`, `dose` or times; a
  `ka` equal to a disposition rate constant (where the closed-form
  solution is singular); and, for the steady-state functions, an
  infusion or zero-order absorption (`tinf`, `dur`) longer than `tau`,
  which the steady-state solutions do not cover.

- [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and the `calc_derived_*()` functions now stop when vector inputs have
  different lengths (other than 1) rather than recycling them.

- [`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md)
  now labels the rows and columns of the variance-covariance matrix with
  the parameter names from the `.cov` file (`THETA1`, `SIGMA(1,1)`,
  `OMEGA(2,1)`, …). Rows were previously numbered and column names were
  altered (`SIGMA.1.1.`). As a result, the columns of
  [`sample_uncert()`](https://kestrel99.github.io/pmxTools/reference/sample_uncert.md)
  output are named after the parameters instead of `V1`, `V2`, …

- [`pk_curve()`](https://kestrel99.github.io/pmxTools/reference/pk_curve.md)
  requires a single value for each parameter in `params`, since every
  dose uses the same parameter set. Vectors previously mixed parameter
  sets across doses.

### Bug fixes

- [`pk_curve()`](https://kestrel99.github.io/pmxTools/reference/pk_curve.md)
  now superposes additional doses by the time elapsed since each dose,
  so irregular time grids, grids not starting at zero and `addl` doses
  beyond the observed window give correct concentrations. A dose
  scheduled at an observation time is counted there even when rounding
  makes the computed dose time slightly later (e.g. `3 * 0.1 > 0.3`).
  Previously doses were added by shifting rows, which was only correct
  for regular grids aligned with `ii`.

- `pk_curve(ss = TRUE)` now gives a steady-state curve (steady state
  reached before the dose at time 0, as for `SS=1` in NONMEM, followed
  by any `addl` doses and then washout). The argument was previously
  ignored.
  [`pk_curve()`](https://kestrel99.github.io/pmxTools/reference/pk_curve.md)
  also validates `model`, `t`, `dose`, `ii`, `addl` and `ss`.

- [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and the `calc_derived_*()` functions calculate `Cmax`, `Tmax`,
  `Cmaxss`, `AUCtau` and their dose-normalized forms separately for each
  parameter set when given vector inputs, returning one value per set.
  Previously one value was returned from the pooled profiles.

- `AUCtau` is now integrated over the whole dosing interval when `step`
  does not divide `tau` exactly.

- `Cmax` and `Tmax` from
  [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and the `calc_derived_*()` functions are now found when the peak comes
  late, e.g. after a long `tlag`. The search window now allows for
  `tlag` and is extended until it contains the peak; previously it could
  end before absorption started and report `Cmax = 0`. The search also
  checks the end of an infusion or of zero-order absorption, so a peak
  there is found exactly rather than at the nearest `step`. A zero dose
  gives `Cmax = 0` and `Tmax = 0` directly.

- The steady-state lag-time functions
  (`calc_ss_*cmt_linear_oral_1_lag()` and
  `calc_ss_*cmt_linear_oral_0_lag()`) are now correct when `tlag`
  exceeds the time after dose by more than one dosing interval, and when
  zero-order absorption runs past the end of the interval. Both
  previously gave wrong, sometimes negative, concentrations.

- [`sample_uncert()`](https://kestrel99.github.io/pmxTools/reference/sample_uncert.md)
  matches the parameter estimates to the covariance matrix by name
  rather than by position, and returns a one-row data frame when
  `n = 1`.

- [`calc_sd_1cmt()`](https://kestrel99.github.io/pmxTools/reference/calc_sd_1cmt.md),
  [`calc_sd_2cmt()`](https://kestrel99.github.io/pmxTools/reference/calc_sd_2cmt.md),
  [`calc_sd_3cmt()`](https://kestrel99.github.io/pmxTools/reference/calc_sd_3cmt.md),
  [`calc_ss_1cmt()`](https://kestrel99.github.io/pmxTools/reference/calc_ss_1cmt.md),
  [`calc_ss_2cmt()`](https://kestrel99.github.io/pmxTools/reference/calc_ss_2cmt.md)
  and
  [`calc_ss_3cmt()`](https://kestrel99.github.io/pmxTools/reference/calc_ss_3cmt.md)
  are now implemented. They choose the regimen-specific function from
  the parameters given (`tinf`, `dur`, `ka`, `tlag`); previously they
  stopped with “Not yet implemented”.

- [`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md)
  joins each per-subject (`FIRSTONLY`) table onto the per-record tables
  by `ID`, whatever the order of its rows and however many table files
  there are, keeping all per-record rows in their original order. It
  warns about subjects missing from a `FIRSTONLY` table and stops if
  there is no shared `ID` column. Per-record tables are combined row by
  row only after checking that their shared columns agree (`ID` exactly,
  other columns such as `TIME` allowing for different `$TABLE`
  precision), and those columns are kept once. Previously tables could
  be combined in the wrong row order without warning.

- [`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md)
  reads tables with exactly three columns and NSIM \> 1 output
  correctly; both previously lost the header row.
  [`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md)
  and
  [`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md)
  no longer turn a one-row table into a one-column table.

- [`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md)
  now parses the XML document instead of requiring its last line to be
  `</nm:output>`, so files with trailing blank lines are read.
  Unparseable files still give a warning and `NULL`.

- [`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)
  and
  [`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md)
  no longer pick up objects with the same names from the global
  environment, and
  [`rnm()`](https://kestrel99.github.io/pmxTools/reference/rnm.md) no
  longer deletes `THETA`, `ETAD` or `EPSD` objects from the global
  environment.

- [`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)
  records the model directory of each step for any path style;
  previously only Windows paths with an upper-case drive letter were
  recognized and others were `NA`.

- [`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md)
  now draws SCMs that stopped after a single forward step.

- [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and the `calc_derived_*()` functions now accept several output names
  in `type`, as documented, returning a list of those outputs.

- `read_nmtables(quiet = TRUE)` no longer reports table files that were
  not found.

- `ggplot2` is now imported selectively instead of in full.

## pmxTools 1.6

CRAN release: 2026-10-07

- Added
  [`cut_quantile()`](https://kestrel99.github.io/pmxTools/reference/cut_quantile.md),
  which creates quantile-based bins for one or more continuous
  variables, with per-variable control over the number of groups.

- Added
  [`datamap()`](https://kestrel99.github.io/pmxTools/reference/datamap.md),
  which plots individual dosing and observation records over time.

- [`dgr_table()`](https://kestrel99.github.io/pmxTools/reference/dgr_table.md)
  now takes `fields` as a named character vector mapping column names to
  display labels. The separate `names` argument is deprecated, but
  unnamed `fields` with `names` continues to work.

- Added dependencies `purrr` and `tidyr`.

- Fixed
  [`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
  which passed the file extension to
  [`xml2::read_xml()`](http://xml2.r-lib.org/reference/read_xml.md) as
  an encoding rather than appending it to the file name, and read the
  unresolved file name instead of the one located on disk. The
  `directory` argument and file names given without the `.xml` extension
  were therefore ignored, and macOS reported an
  `Unsupported encoding: .xml` warning
  ([\#36](https://github.com/kestrel99/pmxTools/issues/36)).

- Removed the blanket `stats` and `utils` imports, which caused a
  `replacing previous import 'stats::filter' by 'dplyr::filter'` warning
  when the package was loaded
  ([\#36](https://github.com/kestrel99/pmxTools/issues/36)).

- [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and the `calc_derived_*()` functions now return `AUCtau` and
  `AUCtau_dose_normalized` as plain numbers. PKNCA 0.12.1.9000 attaches
  a `method` attribute to the AUC it returns, which propagated into the
  derived-parameter list and broke comparisons against plain numeric
  values ([\#37](https://github.com/kestrel99/pmxTools/issues/37)).

- Test snapshots are now committed, the CI workflows have been
  refreshed, and the documentation has been regenerated with roxygen2
  8.1.0 ([\#36](https://github.com/kestrel99/pmxTools/issues/36)).

- [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and the `calc_derived_*()` functions now work with vector inputs
  ([\#29](https://github.com/kestrel99/pmxTools/issues/29)).

## pmxTools 1.5

CRAN release: 2025-08-25

- `plot_dist` has been rewritten to adapt to `ggplot2` 4.0.0 changes,
  and to make it more intuitive to interpret

- Dependency `gghalves` has been removed

- Many fixes and updates to documentation

## pmxTools 1.4

- [`dgr_table()`](https://kestrel99.github.io/pmxTools/reference/dgr_table.md)
  now provides geometric means by default

- [`gm()`](https://kestrel99.github.io/pmxTools/reference/gm.md) now
  provides options to strip `NA` values (`na.rm`) and non-positive
  values (`neg.rm`) before computation.

## pmxTools 1.3

CRAN release: 2023-02-21

- Added NCA parameter estimation to
  [`calc_derived_1cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md),
  [`calc_derived_2cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)
  and
  [`calc_derived_3cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md),
  if dose and other required information (e.g. `tinf`, `dur`, `tau`) is
  provided.

- Added a warning about flip-flop kinetics and their potential effects
  on derived half-lives to
  [`calc_derived_1cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md),
  [`calc_derived_2cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md),
  and
  [`calc_derived_3cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md).

- Added
  [`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md)
  and
  [`blq_log_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md)
  for assistance with visualization of measurements below a limit of
  quantification.

- Added a test for valid XML to
  [`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md).

- Small fixes and enhancements to `read_nmtable_single()` and
  [`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md)
  to ensure that stray text is properly handled while reading in NONMEM
  output files.

- [`count_na()`](https://kestrel99.github.io/pmxTools/reference/count_na.md)
  now throws a warning if `NaN` values are included in the `NA` count,
  and indicates how many of them there are.

- Changed the significant digits functionality for 95% confidence
  intervals for
  [`get_theta()`](https://kestrel99.github.io/pmxTools/reference/get_theta.md),
  [`get_omega()`](https://kestrel99.github.io/pmxTools/reference/get_omega.md)
  and
  [`get_sigma()`](https://kestrel99.github.io/pmxTools/reference/get_sigma.md),
  and all functions that apply these.

- Replaced `gridExtra` with `patchwork`.

- Numerous small fixes.

- Removed some unnecessary dependencies.

- Added many new unit tests, increasing test coverage considerably.

## pmxTools 1.2.4

- Fixed some minor documentation issues. Thanks to Julien Grassot for
  spotting these.

## pmxTools 1.2.3

CRAN release: 2022-04-06

- Fixed another error in
  [`sample_uncert()`](https://kestrel99.github.io/pmxTools/reference/sample_uncert.md)
  which was still crashing the function.

- Added helper functions
  [`count_na()`](https://kestrel99.github.io/pmxTools/reference/count_na.md),
  [`dgr_table()`](https://kestrel99.github.io/pmxTools/reference/dgr_table.md),
  [`gcv()`](https://kestrel99.github.io/pmxTools/reference/gcv.md),
  [`fmt_signif()`](https://kestrel99.github.io/pmxTools/reference/fmt_signif.md).

- Added distribution plotting function
  [`plot_dist()`](https://kestrel99.github.io/pmxTools/reference/plot_dist.md).

## pmxTools 1.2.2

- Fixed an error in
  [`sample_uncert()`](https://kestrel99.github.io/pmxTools/reference/sample_uncert.md)
  which was crashing the function.

## pmxTools 1.2.1

CRAN release: 2020-08-26

- Added vignette describing PK curves.

- Fixed a systematic error in
  [`calc_sd_1cmt_linear_oral_0_lag()`](https://kestrel99.github.io/pmxTools/reference/calc_sd_1cmt.md),
  [`calc_ss_1cmt_linear_oral_0_lag()`](https://kestrel99.github.io/pmxTools/reference/calc_ss_1cmt.md),
  [`calc_sd_2cmt_linear_oral_0_lag()`](https://kestrel99.github.io/pmxTools/reference/calc_sd_2cmt.md)
  which resulted in incorrect curves being plotted.

## pmxTools 1.2

- Rewrote `plot_scm` function to produce tree diagrams via `DiagrammeR`.

- Fixed a rare bug in `get_auc` in which measurements from different
  individuals could be erroneously mixed.

- Amended (or added)
  [`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
  [`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
  [`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
  [`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
  [`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md)
  functions to allow reading of multiple NONMEM estimation steps.

- Clarified documentation for
  [`calc_derived()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md).

- Now using `xml2` to read NONMEM-generated XML files.

- Fixed partial argument match warnings encountered during tests.

- Add function for calculting geometric CV.
