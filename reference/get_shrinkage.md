# Extract shrinkage estimates from a NONMEM output object.

Extract shrinkage estimates from a NONMEM output object.

## Usage

``` r
get_shrinkage(x, output = "eta", type = "sd", sigdig = 3, est.step = NULL)
```

## Arguments

- x:

  A NONMEM output object generated using
  [`read_nm`](https://kestrel99.github.io/pmxTools/reference/read_nm.md).

- output:

  A flag specifying the shrinkage estimates to be output. Valid flag
  values are `eta` (the default), `epsilon`, or `all`.

- type:

  Specifies the type of shrinkage to report. Valid values are `sd`
  (standard deviation, the default) or `vr` (variance, if present in the
  XML output).

- sigdig:

  Specifies the number of significant digits to be provided (default=3).

- est.step:

  Specifies which estimation step to return parameters from (default is
  the last).

## Value

A named vector of NONMEM shrinkage estimates, or in the case of `all`, a
list of named vectors.

`eta` returns a vector of ETA shrinkages, as reported by NONMEM.
`epsilon` returns EPSILON shrinkage, as reported by NONMEM. `all`
returns both ETA and EPSILON shrinkage estimates as a list of vectors.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 nmOutput <- read_nm("run315.xml")
 shr <- get_shrinkage(nmOutput, output="all")
} # }
```
