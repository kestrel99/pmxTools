# Simulate pediatric body weight and height from NHANES

Simulates a virtual pediatric population with realistic body weight
and/or standing height by smoothed resampling of children in the
National Health and Nutrition Examination Survey (NHANES), see
[nhanes_peds](https://kestrel99.github.io/pmxTools/reference/nhanes_peds.md).

## Usage

``` r
sample_nhanes(
  n = 500,
  ages = 2:17,
  sex = c("Male", "Female"),
  vars = c("WT", "HT"),
  cycles = NULL,
  bandwidth_factor = 1,
  seed = NULL,
  data = pmxTools::nhanes_peds
)
```

## Arguments

- n:

  Number of children to simulate per age x sex stratum.

- ages:

  Ages in whole years (2 to 17 for
  [nhanes_peds](https://kestrel99.github.io/pmxTools/reference/nhanes_peds.md)).

- sex:

  `"Male"`, `"Female"` or both.

- vars:

  Measures to simulate: `"WT"` (body weight, kg), `"HT"` (standing
  height, cm) or both.

- cycles:

  NHANES releases to use (values of `data$CYCLE`); `NULL` uses all of
  them.

- bandwidth_factor:

  Multiplier for the kernel bandwidth. `0` gives plain weighted
  resampling of real children; values above 1 smooth more.

- seed:

  Optional random seed. The previous random number generator state is
  restored afterwards.

- data:

  Reference data with the columns of
  [nhanes_peds](https://kestrel99.github.io/pmxTools/reference/nhanes_peds.md).

## Value

A tibble with one row per simulated child and columns `ID`, `AGE`,
`SEXN` (1 = male, 2 = female), `SEX`, the requested `vars`,
`SOURCE_CYCLE` and `SOURCE_SEQN` (the NHANES release and sequence number
of the donor child). Attributes: `"kernels"`, a tibble of per-stratum
diagnostics (`N_NHANES` donors, `N_EFF` effective sample size, `H_WT`
and/or `H_HT` bandwidth SDs on the log scale, and `RHO`, the weighted
log weight-height correlation, when both measures are simulated);
`"cycles"` and `"vars"`, as used.

## Details

Each age (whole years) x sex stratum is simulated separately:

1.  Donors are the stratum's children with every requested measure in
    `vars`. A child missing height is therefore still used when only
    weight is requested.

2.  Each donor's sampling probability is its MEC exam weight, normalised
    within its NHANES release, times an equal share per release, so that
    no release dominates because of a larger total weight. The equal
    share is a modelling choice, not an official pooled NHANES weight.

3.  `n` donors are drawn with these probabilities and Gaussian noise is
    added on the log scale, so values stay positive and continuous. The
    noise covariance is \$\$H = (b \cdot 0.9)^2 \\ n\_{eff}^{-2/(d+4)}
    \\ S R S\$\$ where `b` is `bandwidth_factor`, `d` is the number of
    `vars`, \\n\_{eff}\\ is Kish's effective sample size, `S` holds each
    log measure's robust scale \\\min(SD, IQR/1.349)\\ and `R` is the
    weighted correlation matrix of the log measures. With one measure
    this is the robust Silverman rule \\h = 0.9 \\ \sigma \\
    n\_{eff}^{-1/5}\\; with both, the correlated noise preserves the
    weight-height relationship.

Use
[`compare_nhanes()`](https://kestrel99.github.io/pmxTools/reference/compare_nhanes.md)
and
[`plot_nhanes()`](https://kestrel99.github.io/pmxTools/reference/plot_nhanes.md)
to check the simulated population against the reference.

## See also

[nhanes_peds](https://kestrel99.github.io/pmxTools/reference/nhanes_peds.md),
[`compare_nhanes()`](https://kestrel99.github.io/pmxTools/reference/compare_nhanes.md),
[`plot_nhanes()`](https://kestrel99.github.io/pmxTools/reference/plot_nhanes.md)

## Examples

``` r
sim <- sample_nhanes(n = 100, ages = c(2, 8, 14), seed = 20261008)
head(sim)
#> # A tibble: 6 × 8
#>      ID   AGE  SEXN SEX      WT    HT SOURCE_CYCLE SOURCE_SEQN
#>   <int> <dbl> <int> <chr> <dbl> <dbl> <chr>              <dbl>
#> 1     1     2     1 Male   14.3  91.7 2013-14            81515
#> 2     2     2     1 Male   13.3  91.7 2015-16            90822
#> 3     3     2     1 Male   12.5  84.4 2021-23           136484
#> 4     4     2     1 Male   15.8  96.2 2015-16            84344
#> 5     5     2     1 Male   10.6  79.8 2013-14            73710
#> 6     6     2     1 Male   13.0  85.9 2015-16            85690
attr(sim, "kernels")
#> # A tibble: 6 × 7
#>     AGE SEX    N_NHANES N_EFF   H_WT   H_HT   RHO
#>   <dbl> <chr>     <int> <dbl>  <dbl>  <dbl> <dbl>
#> 1     2 Male        366  263. 0.0429 0.0171 0.735
#> 2     2 Female      381  238. 0.0474 0.0168 0.752
#> 3     8 Male        408  275. 0.0835 0.0153 0.703
#> 4     8 Female      364  245. 0.0909 0.0176 0.690
#> 5    14 Male        327  208. 0.104  0.0188 0.637
#> 6    14 Female      333  204. 0.0814 0.0137 0.435

# Weight only, using every child with a measured weight
wt <- sample_nhanes(n = 100, ages = 10, vars = "WT", seed = 1)
summary(wt$WT)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   21.25   32.58   39.43   42.18   48.75  116.05 
```
