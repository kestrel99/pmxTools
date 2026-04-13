# Sample from the multivariate normal distribution to generate new sets of parameters from NONMEM output.

Sample from the multivariate normal distribution to generate new sets of
parameters from NONMEM output.

## Usage

``` r
sample_uncert(nmRun, n, seed)
```

## Arguments

- nmRun:

  Root filename for the NONMEM run (e.g. "run315.xml").

- n:

  Number of samples required.

- seed:

  Random seed.

## Value

A data frame containing `n` samples from the multivariate normal
distribution, using NONMEM typical parameter estimates the NONMEM
variance-covariance matrix (from the \*.cov file). This provides `n`
sets of parameter estimates sampled from the uncertainty distribution,
suitable for simulation under model uncertainty.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 nmMatrix <- sample_uncert("run315.xml", 5000, seed=740727)
} # }
```
