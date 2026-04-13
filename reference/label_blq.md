# Label axes with censoring labels for BLQ

For ggplot2 scales.

## Usage

``` r
label_blq(lloq, lloq_text)
```

## Arguments

- lloq:

  The value of the lower limit of quantification as a numeric scalar

- lloq_text:

  The text to use on the axis to indicate values `< lloq`. It will be
  automatically set to `paste0("<", lloq)` if missing.

## Value

A function of `x` which returns the formatted values.

## See also

Other BLQ Transformation:
[`blq_trans()`](https://kestrel99.github.io/pmxTools/reference/blq_trans.md),
[`breaks_blq_general()`](https://kestrel99.github.io/pmxTools/reference/breaks_blq_general.md),
[`estimate_lloq()`](https://kestrel99.github.io/pmxTools/reference/estimate_lloq.md),
[`ftrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/ftrans_blq_linear.md),
[`itrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/itrans_blq_linear.md)
