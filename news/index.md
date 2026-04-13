# Changelog

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
