# Polyserial correlation

Estimates the correlation between a continuous variable and an ordinal
variable, assuming the ordinal variable arises from cutting an
underlying normal variable into categories. Uses the two-step estimator:
thresholds from the marginal proportions, then maximum likelihood for
the correlation (Olsson, Drasgow and Dorans 1982).

## Usage

``` r
polyserial(x, y, max_cor = 0.9999)
```

## Arguments

- x:

  A numeric vector.

- y:

  A factor of the same length as `x`. Its level order is used as the
  ordinal scale; an unordered factor gives a warning.

- max_cor:

  Largest absolute correlation searched, between 0 and 1.

## Value

The estimated correlation, using cases complete for both variables, or
`NA` if `y` has fewer than two observed levels or `x` is constant.

## References

Olsson U, Drasgow F, Dorans NJ (1982). The polyserial correlation
coefficient. *Psychometrika*, 47, 337-347.

## See also

[`het_cor()`](https://kestrel99.github.io/pmxTools/reference/het_cor.md),
[`polychoric()`](https://kestrel99.github.io/pmxTools/reference/polychoric.md)

## Examples

``` r
set.seed(1)
x <- rnorm(500)
latent <- 0.6 * x + sqrt(1 - 0.6^2) * rnorm(500)
y <- cut(latent, c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
polyserial(x, y)
#> [1] 0.5373007
```
