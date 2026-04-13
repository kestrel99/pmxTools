# Forward transformation for linear BLQ data

For ggplot2 scales.

## Usage

``` r
ftrans_blq_linear(lloq, multiplier)

ftrans_blq_log(lloq, multiplier, base = 10)
```

## Arguments

- lloq:

  The value of the lower limit of quantification as a numeric scalar

- multiplier:

  When data are `< lloq`, they are replaced by `lloq*multiplier` for
  display.

- base:

  The base for the logarithm

## Value

A function of `x` that replaces `x < lloq` with `lloq*multiplier`

## Functions

- `ftrans_blq_log()`: Log-scale transformation

## See also

Other BLQ Transformation:
[`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md),
[`breaks_blq_general()`](https://kestrel99.github.io/pmxTools/reference/breaks_blq_general.md),
[`estimate_lloq()`](https://kestrel99.github.io/pmxTools/reference/estimate_lloq.md),
[`itrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/itrans_blq_linear.md),
[`label_blq()`](https://kestrel99.github.io/pmxTools/reference/label_blq.md)
