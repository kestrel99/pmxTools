# Create a table of model parameter estimates from a NONMEM output object.

Create a table of model parameter estimates from a NONMEM output object.

## Usage

``` r
get_est_table(
  x,
  thetaLabels = c(),
  omegaLabels = c(),
  sigmaLabels = c(),
  sigdig = 3
)
```

## Arguments

- x:

  A NONMEM output object generated using
  [`read_nm`](https://kestrel99.github.io/pmxTools/reference/read_nm.md).

- thetaLabels:

  A vector containing labels for THETA parameters.

- omegaLabels:

  A vector containing labels for OMEGA parameters.

- sigmaLabels:

  A vector containing labels for SIGMA parameters.

- sigdig:

  The desired number of significant digits to display.

## Value

A named vector of NONMEM model parameter estimates.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 nmOutput <- read_nm("run315.xml")
 estTab   <- get_est_table(nmOutput)
} # }
```
