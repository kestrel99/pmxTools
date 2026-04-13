# Read NONMEM output into a list.

`read_nmext` returns a summary of a given NONMEM run, including
termination messages, parameter estimates, and precision estimates.
Minimally, the NONMEM output and '.ext' files must be available.

## Usage

``` r
read_nmext(
  fileName,
  fileExt = ".lst",
  directory = NULL,
  quiet = FALSE,
  estNo = NULL,
  ...
)
```

## Arguments

- fileName:

  A NONMEM output file prefix, without extension (e.g. "run315").

- fileExt:

  The file extension for NONMEM output, set to ".lst" by default.

- directory:

  The directory to look for files within. If NULL, uses the current
  directory.

- quiet:

  Flag for displaying intermediate output.

- estNo:

  The estimation number to report (by default, if only one estimation
  step is present, that will be reported; if multiple are reported, the
  last will be reported by default).

- ...:

  Passed to each of the read functions (ignored in the functions).

## Value

A list of lists, containing 'Termination' (summary of NONMEM's
termination output, including shrinkages and ETABAR estimates), 'OFV'
(the objective function value), 'Thetas' (a vector of structural
parameter estimates, or THETAs), 'Omega', a list of lists containing the
OMEGA matrix, 'Sigma', a list of lists containing the SIGMA matrix,
'seThetas', a vector of standard errors for THETAs, 'seOmega', a list of
lists containing standard errors for the OMEGA matrix, and 'seSigma', a
list of lists containing standard errors for the SIGMA matrix.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
read_nmext("run315")
read_nmext("run315", ".nmlst")
} # }
```
