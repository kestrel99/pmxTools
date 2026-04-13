# Read a standard NONMEM extension file

Read a standard NONMEM extension file

## Usage

``` r
read_nm_std_ext(fileName, extension, directory = NULL, ...)
```

## Arguments

- fileName:

  The filename (with directory name, if applicable) to read (with or
  without the extension)

- extension:

  The file extension to optionally append (preferably starting with a
  ".")

- directory:

  The directory to look for files within. If NULL, uses the current
  directory.

- ...:

  Passed to
  [`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md)

## Value

NULL if the file does not exist or the value of
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md)
if it does exist.

## Examples

``` r
if (FALSE) { # \dontrun{
read_nm_std_ext("run1", "phi")
} # }
```
