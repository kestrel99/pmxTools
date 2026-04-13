# Generate breaks for measurements below the limit of quantification

Breaks that are `< lloq` are removed. If the lowest break is removed if
it is too close to the lloq.

## Usage

``` r
breaks_blq_general(lloq, breakfun, trans = identity, ...)
```

## Arguments

- lloq:

  The value of the lower limit of quantification as a numeric scalar

- breakfun:

  The function used for normal scale breaks if the `lloq` were not
  present.

- trans:

  A parameter translation function (typically either `identity` for
  linear scale or `log` for log scale).

- ...:

  passed as `breakfun(n=n, ...)`

## Value

A function for calculating breaks with arguments `x` and `n`

## Details

For ggplot2 scales. This is not usually used directly. See
[`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md)
and `blq_log10_trans()` for the functions that are more commonly used.

## See also

Other BLQ Transformation:
[`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md),
[`estimate_lloq()`](https://kestrel99.github.io/pmxTools/reference/estimate_lloq.md),
[`ftrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/ftrans_blq_linear.md),
[`itrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/itrans_blq_linear.md),
[`label_blq()`](https://kestrel99.github.io/pmxTools/reference/label_blq.md)

## Examples

``` r
breaks_blq_general(lloq=3, breakfun=scales::breaks_extended)(1:100, n=5)
#> [1]   3  25  50  75 100
```
