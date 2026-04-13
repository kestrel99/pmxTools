# Extract variability parameter estimates from a NONMEM output object.

Extract variability parameter estimates from a NONMEM output object.

## Usage

``` r
get_omega(x, output = "est", sigdig = 6, sep = "-", est.step = NULL)
```

## Arguments

- x:

  A NONMEM output object generated using
  [`read_nm`](https://kestrel99.github.io/pmxTools/reference/read_nm.md).

- output:

  A flag specifying the matrix or matrices to be output. Valid flag
  values are `est` (the default), `se`, `rse`, `cor`, `cse`, `95ci`, or
  `all`.

- sigdig:

  Specifies the number of significant digits to be provided (default=6).

- sep:

  Specifies the separator character to use for 95% confidence intervals
  (default="-").

- est.step:

  Specifies which estimation step to return parameters from (default is
  the last).

## Value

A symmetrical matrix, or a list of symmetrical matrices if `all` is
specified.

## Details

Output options are as follows:

- *est* returns the estimated `OMEGA` variance-covariance matrix.

- *se* returns the standard errors for the estimated `OMEGA`
  variance-covariance matrix.

- *rse* returns the relative standard errors for the estimated `OMEGA`
  variance-covariance matrix (`se/est*100`).

- *cor* returns the correlation matrix matrix.

- *cse* returns the standard errors for the correlation matrix.

- *95ci* returns the asymptotic 95% confidence intervals for the
  elements of the `OMEGA` variance-covariance matrix
  (`est +/- 1.96*se`).

- *all* returns all available `OMEGA` matrices.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 nmOutput  <- read_nm("run315.xml")
 omegas    <- get_omega(nmOutput)
 omegaRSEs <- get_omega(nmOutput, "rse")
} # }
```
