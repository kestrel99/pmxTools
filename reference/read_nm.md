# Read NONMEM 7.2+ output into a list of lists.

Read NONMEM 7.2+ output into a list of lists.

## Usage

``` r
read_nm(fileName, directory = NULL, quiet = FALSE, ...)
```

## Arguments

- fileName:

  A NONMEM XML output file (e.g. "run315.xml").

- directory:

  The directory to look for files within. If NULL, uses the current
  directory.

- quiet:

  Flag for displaying intermediate output.

- ...:

  Passed to each of the read functions (ignored in the functions).

## Value

A list of lists corresponding to a NONMEM output object.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
nmOutput <- read_nm("run315.xml")
} # }
```
