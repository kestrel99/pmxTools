# Sample from the multivariate normal distribution using the SIGMA variance-covariance matrix to generate new sets of simulated EPSILONs from NONMEM output.

Sample from the multivariate normal distribution using the SIGMA
variance-covariance matrix to generate new sets of simulated EPSILONs
from NONMEM output.

## Usage

``` r
sample_sigma(nmRun, n, seed)
```

## Arguments

- nmRun:

  Root filename for the NONMEM run (e.g. "run315").

- n:

  Number of samples required.

- seed:

  Random seed.

## Value

A data frame containing `n` samples from the multivariate normal
distribution, using the estimated NONMEM SIGMA variance-covariance
matrix. This provides `n` sets of EPSILON estimates suitable for
simulation of new datasets.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 sigDist <- sample_sigma("run315", 5000, seed=740727)
} # }
```
