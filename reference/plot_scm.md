# Visualize PsN SCM output.

`plot_scm` returns a visualization of a Perl-speaks-NONMEM (PsN,
<https://uupharmacometrics.github.io/PsN/>) SCM (stepwise covariate
modeling) procedure. It depends on the presence of `scmlog.txt` and
`short_scmlog.txt` files in the specified directory.

## Usage

``` r
plot_scm(
  dir,
  startPhase = "forward",
  fwdSuccessCol = "#66C2A5",
  fwdFailCol = "black",
  bwdSuccessCol = "#FC8D62",
  bwdFailCol = "black",
  defCol = "black",
  fwdSuccessFillCol = "#B3E2CD",
  fwdFailFillCol = "white",
  bwdSuccessFillCol = "#FDCDAC",
  bwdFailFillCol = "white",
  defFillCol = "white",
  fwdSuccessFontCol = "black",
  fwdFailFontCol = "black",
  bwdSuccessFontCol = "black",
  bwdFailFontCol = "black",
  defFontCol = "black",
  fullFwdCol = "#8DA0CB",
  finalCol = "#E78AC3",
  fullFwdFillCol = "#CBD5E8",
  finalFillCol = "#F4CAE4",
  fullFwdFontCol = "black",
  finalFontCol = "black",
  fullFwdWidth = "2px",
  finalWidth = "2px",
  defWidth = "1px",
  nodeStyle = "filled,rounded",
  nodeShape = "box",
  fontname = "helvetica",
  rankdir = "TB",
  layout = "dot",
  lookupDF = NULL,
  ...
)
```

## Arguments

- dir:

  A PsN SCM folder (containing `scmlog.txt` and `short_scmlog.txt`).

- startPhase:

  Where to start collating the output; can be `"forward"` (the default)
  or `"backward"`.

- fwdSuccessCol:

  Node outline color for a model fit matching the forward inclusion
  criterion.

- fwdFailCol:

  Node outline color for a model fit not matching the forward inclusion
  criterion.

- bwdSuccessCol:

  Node outline color for a model fit matching the backward elimination
  criterion.

- bwdFailCol:

  Node outline color for a model fit not matching the backward
  elimination criterion.

- defCol:

  Default node outline color.

- fwdSuccessFillCol:

  Node fill color for a model fit matching the forward inclusion
  criterion.

- fwdFailFillCol:

  Node fill color for a model fit not matching the forward inclusion
  criterion.

- bwdSuccessFillCol:

  Node fill color for a model fit matching the backward elimination
  criterion.

- bwdFailFillCol:

  Node fill color for a model fit not matching the backward elimination
  criterion.

- defFillCol:

  Default node fill color.

- fwdSuccessFontCol:

  Node font color for a model fit matching the forward inclusion
  criterion.

- fwdFailFontCol:

  Node font color for a model fit not matching the forward inclusion
  criterion.

- bwdSuccessFontCol:

  Node font color for a model fit matching the backward elimination
  criterion.

- bwdFailFontCol:

  Node font color for a model fit not matching the backward elimination
  criterion.

- defFontCol:

  Default node font color.

- fullFwdCol:

  Node outline color for the full forward model (i.e. the final model
  before the backward elimination procedure in SCM).

- finalCol:

  Node outline color for the final reduced model (i.e. the final model
  reached after the backward elimination procedure in SCM).

- fullFwdFillCol:

  Node fill color for the full forward model (i.e. the final model
  before the backward elimination procedure in SCM).

- finalFillCol:

  Node fill color for the final reduced model (i.e. the final model
  reached after the backward elimination procedure in SCM).

- fullFwdFontCol:

  Node font color for the full forward model (i.e. the final model
  before the backward elimination procedure in SCM).

- finalFontCol:

  Node font color for the final reduced model (i.e. the final model
  reached after the backward elimination procedure in SCM).

- fullFwdWidth:

  Node outline width for the full forward model (i.e. the final model
  before the backward elimination procedure in SCM).

- finalWidth:

  Node outline width for the final reduced model (i.e. the final model
  reached after the backward elimination procedure in SCM).

- defWidth:

  Default node outline width.

- nodeStyle:

  Node style. A string containing a comma-separated list of options
  (which include "filled", "striped", "wedged", "diagonals" and
  "rounded"). See the GraphViz documentation for further details.

- nodeShape:

  Node shape. Options include "box" (the default), "oval", "diamond",
  "egg", "plaintext", "point", "square", "triangle" and many more. See
  the GraphViz documentation for further details.

- fontname:

  Font for nodes. Options depend heavily on the local system - see the
  GraphViz documentation for further details.

- rankdir:

  Direction of graph layout. Possible values are "TB" (the default),
  "LR", "BT", "RL", corresponding to directed graphs drawn from top to
  bottom, from left to right, from bottom to top, and from right to
  left, respectively.

- layout:

  Graph layout. Possible values are "dot" (the default), "neato",
  "twopi", and "circo". Note that of these, "dot" is the easiest to
  interpret and the others may produce odd results.

- lookupDF:

  A data frame containing a lookup table for node labels. By default,
  this function will use the PSN model names. If a lookup table
  containing the fields `Model` and `Alias` is provided, model names in
  `Model` will be replaced in the output plots by attaching labels in
  `Alias`.

- ...:

  Additional parameters passed to the underlying
  [`SetNodeStyle`](https://rdrr.io/pkg/data.tree/man/ToDiagrammeRGraph.html)
  and
  [`SetEdgeStyle`](https://rdrr.io/pkg/data.tree/man/ToDiagrammeRGraph.html)
  functions, which in turn rely on `DiagrammeR`.

## Value

A `grViz` object.

## Details

This function parses PsN SCM output and displays it as a GraphViz graph
(effectively, an HTML widget). It is built on
[`plot.Node`](https://rdrr.io/pkg/data.tree/man/ToDiagrammeRGraph.html) -
please refer to documentation for this function for a more detailed
overview of what is possible (a lot). For more specific details, see
<https://rich-iannone.github.io/DiagrammeR/>.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

GraphViz (<https://graphviz.org/Documentation.php>)

Lindbom L, Ribbing J & Jonsson EN (2004). Perl-speaks-NONMEM (PsN) - A
Perl module for NONMEM related programming. Computer Methods and
Programs in Biomedicine, 75(2), 85-94.
[doi:10.1016/j.cmpb.2003.11.003](https://doi.org/10.1016/j.cmpb.2003.11.003)

Lindbom L, Pihlgren P & Jonsson N (2005). PsN-Toolkit - A collection of
computer intensive statistical methods for non-linear mixed effect
modeling using NONMEM. Computer Methods and Programs in Biomedicine,
79(3), 241-257.
[doi:10.1016/j.cmpb.2005.04.005](https://doi.org/10.1016/j.cmpb.2005.04.005)

Other NONMEM reading:
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md),
[`read_scm()`](https://kestrel99.github.io/pmxTools/reference/read_scm.md)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
scm <- plot_scm("E:/DrugX/ModelDevelopment/scm310")
} # }
```
