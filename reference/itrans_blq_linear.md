# Inverse transformation for linear BLQ data

For ggplot2 scales.

## Usage

``` r
itrans_blq_linear(lloq)

itrans_blq_log(lloq, base)
```

## Arguments

- lloq:

  The value of the lower limit of quantification as a numeric scalar

- base:

  The base for the logarithm

## Value

A function of `x` that replaces `x < lloq` with `lloq`

## Functions

- `itrans_blq_log()`: Log-scale inverse transform

## See also

Other BLQ Transformation:
[`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md),
[`breaks_blq_general()`](https://kestrel99.github.io/pmxTools/reference/breaks_blq_general.md),
[`estimate_lloq()`](https://kestrel99.github.io/pmxTools/reference/estimate_lloq.md),
[`ftrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/ftrans_blq_linear.md),
[`label_blq()`](https://kestrel99.github.io/pmxTools/reference/label_blq.md)
