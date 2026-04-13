# Extract structural model parameter estimates and associated information from a NONMEM output object.

Extract structural model parameter estimates and associated information
from a NONMEM output object.

## Usage

``` r
get_theta(x, output = "est", sigdig = 6, sep = "-", est.step = NULL)
```

## Arguments

- x:

  A NONMEM output object generated using
  [`read_nm`](https://kestrel99.github.io/pmxTools/reference/read_nm.md).

- output:

  A flag specifying the matrix or matrices to be output. Valid flag
  values are `est` (the default), `se`, `rse`, `95ci`, or `all`.

- sigdig:

  Specifies the number of significant digits to be provided (default=6).

- sep:

  Specifies the separator character to use for 95% confidence intervals
  (default="-").

- est.step:

  Specifies which estimation step to return parameters from (default is
  the last).

## Value

A named vector of NONMEM model parameter estimates, or in the case of
`all`, a list of named vectors.

## Details

Output options are as follows: *est* returns a vector of `THETA` values.
*se* returns a vector of `THETA` standard errors. *rse* returns a vector
of `THETA` relative standard errors (`se/est*100`). *95ci* returns a
vector of the asymptotic 95% confidence intervals for the elements of
`THETA` (`est +/- 1.96*se`). *all* returns all available `THETA`
information as a list of named vectors.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 nmOutput <- read_nm("run315.xml")
 thetas <- get_theta(nmOutput)
} # }
```
