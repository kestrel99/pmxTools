# Polychoric correlation

Estimates the correlation between two ordinal variables, assuming each
arises from cutting an underlying normal variable into categories. Uses
the two-step estimator: thresholds from the marginal proportions, then
maximum likelihood for the correlation (Olsson 1979).

## Usage

``` r
polychoric(x, y, max_cor = 0.9999)
```

## Arguments

- x, y:

  Factors of the same length. Their level order is used as the ordinal
  scale; unordered factors give a warning.

- max_cor:

  Largest absolute correlation searched, between 0 and 1.

## Value

The estimated correlation, using cases complete for both variables, or
`NA` if either has fewer than two observed levels.

## References

Olsson U (1979). Maximum likelihood estimation of the polychoric
correlation coefficient. *Psychometrika*, 44, 443-460.

## See also

[`het_cor()`](https://kestrel99.github.io/pmxTools/reference/het_cor.md),
[`polyserial()`](https://kestrel99.github.io/pmxTools/reference/polyserial.md)

## Examples

``` r
set.seed(1)
z1 <- rnorm(500)
z2 <- 0.6 * z1 + sqrt(1 - 0.6^2) * rnorm(500)
x <- cut(z1, c(-Inf, -0.5, 0.5, Inf), ordered_result = TRUE)
y <- cut(z2, c(-Inf, 0, 1, Inf), ordered_result = TRUE)
polychoric(x, y)
#> [1] 0.5784943
```
