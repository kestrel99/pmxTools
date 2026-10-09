# Plot simulated against NHANES percentiles

Plots the 5th, 50th and 95th percentiles of each simulated measure
against age, for the simulated population and the NHANES reference, by
sex.

## Usage

``` r
plot_nhanes(comparison)
```

## Arguments

- comparison:

  Output of
  [`compare_nhanes()`](https://kestrel99.github.io/pmxTools/reference/compare_nhanes.md).

## Value

A ggplot object.

## See also

[`compare_nhanes()`](https://kestrel99.github.io/pmxTools/reference/compare_nhanes.md),
[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md)

## Examples

``` r
sim <- sample_nhanes(n = 200, ages = 2:17, seed = 1)
plot_nhanes(compare_nhanes(sim))
```
