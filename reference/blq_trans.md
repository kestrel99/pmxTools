# A transform for ggplot2 with data that may be below the lower limit of quantification

If the `lloq` is not provided, it will be estimated from the data as the
minimum value above zero.

## Usage

``` r
blq_trans(lloq, x, multiplier = 0.5, lloq_text)

blq_log_trans(lloq, x, multiplier = 0.5, base = 10, lloq_text)
```

## Arguments

- lloq:

  The value of the lower limit of quantification as a numeric scalar

- x:

  (only used if `lloq` is missing), the data for `lloq` estimation.

- multiplier:

  When data are `< lloq`, they are replaced by `lloq*multiplier` for
  display.

- lloq_text:

  The text to use on the axis to indicate values `< lloq`. It will be
  automatically set to `paste0("<", lloq)` if missing.

- base:

  The base for the logarithm

## Value

A "trans" object based on the `scales` package for BLQ data.

## Functions

- `blq_log_trans()`: Log-scale transformation with BLQ

## See also

Other BLQ Transformation:
[`breaks_blq_general()`](https://kestrel99.github.io/pmxTools/reference/breaks_blq_general.md),
[`estimate_lloq()`](https://kestrel99.github.io/pmxTools/reference/estimate_lloq.md),
[`ftrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/ftrans_blq_linear.md),
[`itrans_blq_linear()`](https://kestrel99.github.io/pmxTools/reference/itrans_blq_linear.md),
[`label_blq()`](https://kestrel99.github.io/pmxTools/reference/label_blq.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(data=data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
  geom_point()

ggplot(data=data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
  geom_point() +
  scale_x_continuous(trans=blq_trans(lloq=3))

ggplot(data=data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
  geom_point() +
  scale_x_continuous(trans=blq_log10_trans(lloq=3))
} # }
```
