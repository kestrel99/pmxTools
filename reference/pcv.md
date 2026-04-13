# Calculate percentage coefficient of variation

Calculate percentage coefficient of variation

## Usage

``` r
pcv(x, na.rm = FALSE)
```

## Arguments

- x:

  Numeric vector.

- na.rm:

  A logical value indicating whether `NA` values should be stripped
  before the computation proceeds.

## Value

The percentage coefficient of variation.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
pcv(rnorm(50, 5, 7.56))
#> [1] 189.4419
```
