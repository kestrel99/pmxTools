# Reads NONMEM output tables.

Reads NONMEM output tables.

## Usage

``` r
read_nmtables(
  tableFiles = NULL,
  runNo = NULL,
  tabSuffix = "",
  tableNames = c("sdtab", "mutab", "patab", "catab", "cotab", "mytab", "extra", "xptab"),
  quiet = FALSE,
  directory = NULL,
  output_type = c("data.frame", "list"),
  ...
)
```

## Arguments

- tableFiles:

  NONMEM table files to be read.

- runNo:

  Run number.

- tabSuffix:

  Table file suffix.

- tableNames:

  List of root table names, using the Xpose naming convention as the
  default.

- quiet:

  Flag for displaying intermediate output.

- directory:

  The directory to look for files within. If NULL, uses the current
  directory.

- output_type:

  Should output be a "data.frame" where all results are merged or a
  "list" of data.frames.

- ...:

  Passed to each of the read functions (ignored in the functions).

## Value

A data.frame or list of data.frames depending on the `output_type`
argument.

## Note

Adapted from Xpose 4 (<https://CRAN.R-project.org/package=xpose4>).

## References

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

Jonsson EN, Karlsson MO. Xpose–an S-PLUS based population
pharmacokinetic/pharmacodynamic model building aid for NONMEM. Comput
Methods Programs Biomed. 1999 Jan;58(1):51-64

## See also

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)

## Author

Bill Denney, Justin Wilkins, Niclas Jonsson, Andrew Hooker

## Examples

``` r
if (FALSE) { # \dontrun{
tables <- read_nmtables(runNo=315)
} # }
```
