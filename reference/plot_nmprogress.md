# Plot NONMEM parameter estimation by iteration.

`plot_nmprogress` returns a plot or set of plots showing the evolution
of parameter estimates by iteration.

## Usage

``` r
plot_nmprogress(
  fileName,
  fileExt = ".lst",
  metric = "perc",
  lineCol = "#902C10",
  idlineCol = "black"
)
```

## Arguments

- fileName:

  A NONMEM output file prefix, without extension (e.g. 'run315').

- fileExt:

  The file extension for NONMEM output, set to '.lst' by default.

- metric:

  What to show in the plot. Allowed options are 'est' (the actual
  estimate) or 'perc' (the percentage change in the estimated or OFV
  since estimation began). Default is 'perc'.

- lineCol:

  Line color. Default is '#902C10'.

- idlineCol:

  Identity line color (only used if 'perc' metric is selected). Default
  is black.

## Value

A set of plots.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
plot_nmprogress("run315")
plot_nmprogress("run315", ".nmlst")
} # }
```
