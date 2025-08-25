# pmxTools 1.5

* `plot_dist` has been rewritten to adapt to `ggplot2` 4.0.0 changes

* Dependency `gghalves` has been removed

* Many fixes to documentation

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
