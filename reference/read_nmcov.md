# Read in the NONMEM variance-covariance matrix.

Read in the NONMEM variance-covariance matrix.

## Usage

``` r
read_nmcov(fileName, quiet = FALSE, directory = NULL, ...)
```

## Arguments

- fileName:

  Root filename for the NONMEM run (e.g. "run315").

  This function reads the ".cov" NONMEM output table, and will return an
  error if this is missing.

- quiet:

  Flag for displaying intermediate output.

- directory:

  The directory to look for files within. If NULL, uses the current
  directory.

- ...:

  Passed to each of the read functions (ignored in the functions).

## Value

A symmetrical variance-covariance matrix covering all model parameters.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
nmVcov <- read_nmcov("run315")
} # }
```
