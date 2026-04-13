# Sample from the multivariate normal distribution using the OMEGA variance-covariance matrix to generate new sets of simulated ETAs from NONMEM output.

Sample from the multivariate normal distribution using the OMEGA
variance-covariance matrix to generate new sets of simulated ETAs from
NONMEM output.

## Usage

``` r
sample_omega(nmRun, n, seed)
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
distribution, using the estimated NONMEM OMEGA variance-covariance
matrix. This provides `n` sets of ETA estimates suitable for simulation
of new patients.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 omDist <- sample_omega("run315", 5000, seed=740727)
} # }
```
