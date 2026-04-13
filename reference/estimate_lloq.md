# Estimate the lower limit of quantification (LLOQ) from a vector

Nonnegative values are considered to be above the LLOQ. `NA` values are
ignored.

## Usage

``` r
estimate_lloq(x)
```

## Arguments

- x:

  The numeric vector to use for estimation of the LLOQ

## Value

The lowest, nonzero value from `x`. If all are `NA` or zero, 1 is
returned, and a warning is issued.

## See also

Other BLQ Transformation:
[`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md),
[`breaks_blq_general()`](https://kestrel99.github.io/pmxTools/reference/breaks_blq_general.md),
[`ftrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/ftrans_blq_linear.md),
[`itrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/itrans_blq_linear.md),
[`label_blq()`](https://kestrel99.github.io/pmxTools/reference/label_blq.md)

## Examples

``` r
estimate_lloq(c(NA, 0, 2, 5))
#> [1] 2
```
