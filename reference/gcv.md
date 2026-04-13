# Calculate a geometric coefficient of variation.

Calculate a geometric coefficient of variation.

## Usage

``` r
gcv(x, na.rm = F, neg.rm = F)
```

## Arguments

- x:

  A vector.

- na.rm:

  Flag for removing `NA` values (defaults to `FALSE`).

- neg.rm:

  Flag for removing negative or zero values (defaults to `FALSE`).

## Value

The geometric coefficient of variation of the input vector. If `neg.rm`
is `FALSE` and values \<= 0 are present, `NA` will be returned.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 gcv(myvector)
} # }
```
