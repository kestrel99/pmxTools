# Extract problem and estimation information from a NONMEM output object.

Extract problem and estimation information from a NONMEM output object.

## Usage

``` r
get_probinfo(x, sigdig = 6, est.step = NULL)
```

## Arguments

- x:

  A NONMEM output object generated using
  [`read_nm`](https://kestrel99.github.io/pmxTools/reference/read_nm.md).

- sigdig:

  Specifies the number of significant digits to be provided (default=6).

- est.step:

  Specifies which estimation step to return parameters from (default is
  the last).

## Value

A list with elements `prob_title`, `prog_info`, `lic_info`, `prob_times`
(data frame of start, stop and elapsed times), `prob_info` (data frame
of data set and model dimensions) and `est_info` (data frame of
estimation results for `est.step`).

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Examples

``` r
if (FALSE) { # \dontrun{
 nmOutput <- read_nm("run315.xml")
 probInfo <- get_probinfo(nmOutput)
} # }
```
