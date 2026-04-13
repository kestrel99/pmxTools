# Read NONMEM output into a list.

`table_rtf` generates an RTF table from a data frame.

## Usage

``` r
table_rtf(
  df,
  outFile = NULL,
  rtfFile = TRUE,
  boldHeader = TRUE,
  rowNames = FALSE,
  ...
)
```

## Arguments

- df:

  A data frame.

- outFile:

  A filename for writing the table to. If `NULL`, writes to console.

- rtfFile:

  If `TRUE` (the default), then add RTF tabs to generate a fully
  formatted RTF file.

- boldHeader:

  If `TRUE`, make the header bold.

- rowNames:

  If `TRUE`, include row names in the table. Default is `FALSE`.

- ...:

  Other formatting options for the table body.

## Value

An RTF table based on the data frame provided.

## References

<https://www.r-bloggers.com/2008/10/another-solution-to-the-r-to-word-table-problem/>

## Author

John Johnson, <johndjohnson@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
scm <- read_scm("E:/DrugX/ModelDevelopment/scm310")
myRTF <- table_rtf(scm$forwardSummary)
} # }
```
