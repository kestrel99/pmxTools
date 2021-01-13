# pmxTools >1.2.1

* Added `blq_trans()` and `blq_log_trans()` for assistance with visualization of measurements below a limit of quantification.

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
