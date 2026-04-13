# Format a number with the correct number of significant digits and trailing zeroes.

Format a number with the correct number of significant digits and
trailing zeroes.

## Usage

``` r
fmt_signif(x, digits = 3)
```

## Arguments

- x:

  A vector of numeric values.

- digits:

  The number of significant digits values should have (defaults to 3).

## Value

A string containing the properly-formatted number.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 fmt_signif(c(36.44, 0.0002, 3336.7), digits=3)
} # }
```
