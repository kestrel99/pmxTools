# Read (single or) multiple NONMEM tables from a single file

Read (single or) multiple NONMEM tables from a single file

## Usage

``` r
read_nm_multi_table(
  fileName,
  header = TRUE,
  ...,
  simplify = TRUE,
  table_start_pattern = "^TABLE NO"
)
```

## Arguments

- fileName:

  The filename to read from

- header, ...:

  Arguments passed to read.table

- simplify:

  If a single table is present, return a data.frame instead of a list of
  data.frames?

- table_start_pattern:

  What should be found to start a new table?

## Value

A list of data.frames, or if only one is present and simplify=TRUE, a
data.frame.

## See also

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)

## Author

Bill Denney

## Examples

``` r
if (FALSE) { # \dontrun{
read_nm_multi_table("run1.cov", row.names=1)
} # }
```
