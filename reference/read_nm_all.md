# Read all NONMEM files for a single NONMEM run.

Read all NONMEM files for a single NONMEM run.

## Usage

``` r
read_nm_all(runNo, run_prefix = "run", directory = NULL, quiet = FALSE, ...)
```

## Arguments

- runNo:

  Run number.

- run_prefix:

  The start to the name of the run.

- directory:

  The directory to look for files within. If NULL, uses the current
  directory.

- quiet:

  Flag for displaying intermediate output.

- ...:

  Passed to each of the read functions (ignored in the functions).

## Details

The filename for loading is constructed as `paste(run_prefix, runNo)`.
To load a nonstandard file, simply set one of those values to `NULL`.

## See also

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)
