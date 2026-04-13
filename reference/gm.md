# Calculate geometric mean

Calculate geometric mean

## Usage

``` r
gm(x, na.rm = FALSE, neg.rm = FALSE)
```

## Arguments

- x:

  Numeric vector.

- na.rm:

  Flag for removing `NA` values (defaults to `FALSE`).

- neg.rm:

  Flag for removing negative or zero values (defaults to `FALSE`).

## Value

The geometric mean. `NA` is returned if there are any non-positive
elements in `x`.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
gm(c(0.5, 7, 8, 5))
#> [1] 3.439791
```
