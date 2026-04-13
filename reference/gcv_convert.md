# Convert geometric variance or standard deviation to a geometric coefficient of variation

The equation used is: `100*sqrt(exp(gvar)-1)`

## Usage

``` r
gcv_convert(gvar = gsd^2, gsd)
```

## Arguments

- gvar:

  The geometric variance (note that this is the variance not a vector of
  values to compute the gcv from)

- gsd:

  The geometric standard deviation

## Value

Geometric coefficient of variation

## References

[http://onbiostatistics.blogspot.com/2008/07/geometric-statistics-geometric-cv-vs.html](http://onbiostatistics.blogspot.com/2008/07/geometric-statistics-geometric-cv-vs.md)

## Author

Bill Denney

## Examples

``` r
gcv_convert(0.2)
#> [1] 47.05345
gcv_convert(gsd=0.2)
#> [1] 20.20168
```
