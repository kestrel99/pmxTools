# Compare a simulated pediatric population with NHANES

Summarises body weight, height and (when both were simulated) body mass
index by age and sex for a population from
[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md)
and for the NHANES reference it was drawn from. The reference uses the
same releases and donor eligibility as the simulation, weighted with the
cycle-balanced MEC weights.

## Usage

``` r
compare_nhanes(sim, data = pmxTools::nhanes_peds)
```

## Arguments

- sim:

  Output of
  [`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md),
  with its attributes intact.

- data:

  Reference data with the columns of
  [nhanes_peds](https://kestrel99.github.io/pmxTools/reference/nhanes_peds.md);
  normally the data used for the simulation.

## Value

A tibble with one row per `AGE`, `SEX` and `VARIABLE` (`WT`, `HT` and
`BMI` in kg/m^2, as simulated). Columns `N`, `Mean`, `SD`, `P05`, `Q1`,
`Median`, `Q3` and `P95` appear with suffix `_REF` (weighted NHANES
reference) and `_SIM` (simulated), followed by the ratios `MEAN_RATIO`,
`MEDIAN_RATIO`, `P05_RATIO` and `P95_RATIO` (simulated / reference).
When both weight and height were simulated, `RHO_REF` and `RHO_SIM` give
the correlation of log weight and log height.

## See also

[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md),
[`plot_nhanes()`](https://kestrel99.github.io/pmxTools/reference/plot_nhanes.md)

## Examples

``` r
sim <- sample_nhanes(n = 200, ages = c(4, 12), seed = 1)
compare_nhanes(sim)
#> # A tibble: 12 × 25
#>      AGE SEX    VARIABLE N_REF Mean_REF SD_REF P05_REF Q1_REF Median_REF Q3_REF
#>    <dbl> <chr>  <chr>    <dbl>    <dbl>  <dbl>   <dbl>  <dbl>      <dbl>  <dbl>
#>  1     4 Male   WT         364     18.8   3.24    14.8   16.7       18.1   20.2
#>  2     4 Male   HT         364    106.    4.85    98.7  103.       106.   109. 
#>  3     4 Male   BMI        364     16.5   1.95    14.3   15.3       16.1   17.1
#>  4     4 Female WT         357     18.2   3.32    14.1   15.9       17.6   19.7
#>  5     4 Female HT         357    105.    5.18    97.1  102.       105.   109. 
#>  6     4 Female BMI        357     16.4   1.99    13.9   15.0       16.1   17.0
#>  7    12 Male   WT         344     51.2  15.0     32.9   40.5       47.8   58.6
#>  8    12 Male   HT         344    155.    8.78   141.   149.       155.   161. 
#>  9    12 Male   BMI        344     20.9   4.66    15.4   17.5       19.8   23.6
#> 10    12 Female WT         319     52.6  13.9     32.9   42.3       50.6   60.0
#> 11    12 Female HT         319    155.    7.30   142.   150.       154.   159. 
#> 12    12 Female BMI        319     21.9   5.06    15.5   18.0       21.0   25.0
#> # ℹ 15 more variables: P95_REF <dbl>, N_SIM <dbl>, Mean_SIM <dbl>,
#> #   SD_SIM <dbl>, P05_SIM <dbl>, Q1_SIM <dbl>, Median_SIM <dbl>, Q3_SIM <dbl>,
#> #   P95_SIM <dbl>, MEAN_RATIO <dbl>, MEDIAN_RATIO <dbl>, P05_RATIO <dbl>,
#> #   P95_RATIO <dbl>, RHO_REF <dbl>, RHO_SIM <dbl>
```
