# Plot a distribution as a hybrid containing a halfeye, a boxplot and jittered points.

Plot a distribution as a hybrid containing a halfeye, a boxplot and
jittered points.

## Usage

``` r
plot_dist(
  dat,
  yvar,
  xvar = NULL,
  ylim = NULL,
  xlb = "",
  ylb = "",
  identity_line = FALSE,
  identity_value = 0,
  he_adjust = 0.5,
  he_width = 0.6,
  he_justification = -0.2,
  he_col = "black",
  he_fill = "#F8766D",
  he_alpha = 0.9,
  he_slab_type = "pdf",
  he_breaks = "Sturges",
  he_outline_bars = FALSE,
  he_point_interval = "median_qi",
  he_point_alpha = 0.9,
  he_point_fill = "#F8766D",
  he_point_colour = "#F8766D",
  he_point_size = 2,
  bxp_width = 0.12,
  bxp_outlier_col = NA,
  bxp_outlier_fill = NA,
  bxp_outlier_shape = 19,
  bxp_outlier_size = 1.5,
  bxp_col = "black",
  bxp_fill = "#F8766D",
  bxp_alpha = 0.9,
  bxp_notch = FALSE,
  bxp_notchwidth = 0.5,
  hp_alpha = 0.25,
  hp_col = "#F8766D",
  hp_size = 1,
  hp_shape = 16,
  na.rm = FALSE
)
```

## Arguments

- dat:

  A data frame.

- yvar:

  The name of the field containing values to be plotted.

- xvar:

  The name of the field containing the grouping variable (defaults to
  `NULL`).

- ylim:

  Limits for the y-axis. Defaults to `NULL`. If provided, should be a
  2-element vector containing the upper and lower limits.

- xlb:

  Label for the x-axis.

- ylb:

  Label for the y-axis.

- identity_line:

  Show a line of identity? Default `FALSE`.

- identity_value:

  If an identity line is shown, it will be drawn horizontally at this
  y-value (default 0).

- he_adjust:

  If `he_slab_type` is `"pdf"`, bandwidth for the density estimator is
  adjusted by multiplying it by this value.

- he_width:

  Width of the halfeye component of the plot (default 0.6).

- he_justification:

  Justification of the halfeye component of the plot (default -0.2).

- he_col:

  Color for the halfeye component of the plot.

- he_fill:

  Fill color for the halfeye component of the plot.

- he_alpha:

  Alpha for the halfeye component of the plot (default 0.9).

- he_slab_type:

  The type of slab function to calculate for the halfeye component of
  the plot: probability density (or mass) function (`"pdf"`, the
  default), cumulative distribution function (`"cdf"`), complementary
  CDF (`"ccdf"`) or histogram (`"histogram"`).

- he_breaks:

  If slab_type is `"histogram"`, the breaks parameter that is passed to
  [`hist()`](https://rdrr.io/r/graphics/hist.html) to determine where to
  put breaks in the histogram.

- he_outline_bars:

  If slab_type is `"histogram"`, determines if outlines in between the
  bars are drawn when the slab_color aesthetic is used. If `FALSE` (the
  default), the outline is drawn only along the tops of the bars; if
  `TRUE`, outlines in between bars are also drawn.

- he_point_interval:

  A function from the
  [`ggdist::point_interval`](https://mjskay.github.io/ggdist/reference/point_interval.html)
  family (e.g., `median_qi`, `mean_qi`, `mode_hdi`, etc), or a string
  giving the name of a function from that family (e.g., `"median_qi"`,
  `"mean_qi"`, `"mode_hdi"`, etc. This function determines the point
  summary (typically mean, median, or mode) and interval type (quantile
  interval, `qi`; highest-density interval, `hdi`; or highest-density
  continuous interval, `hdci`). Output will be converted to the
  appropriate x- or y-based aesthetics depending on the value of
  orientation.

- he_point_alpha:

  Alpha for the point.

- he_point_fill:

  Fill colour for the point.

- he_point_colour:

  Colour for the point.

- he_point_size:

  Size for the point.

- bxp_width:

  Width of the boxplot component (default 0.12).

- bxp_outlier_col:

  Color for outliers in the boxplot component.

- bxp_outlier_fill:

  Fill color for outliers in the boxplot component.

- bxp_outlier_shape:

  Shape for outliers in the boxplot component.

- bxp_outlier_size:

  Size for outliers in the boxplot component.

- bxp_col:

  Color for the boxplot component.

- bxp_fill:

  Fill color for the boxplot component.

- bxp_alpha:

  Alpha for the boxplot component.

- bxp_notch:

  If `FALSE` (default) make a standard box plot. If `TRUE`, make a
  notched box plot. Notches are used to compare groups; if the notches
  of two boxes do not overlap, this suggests that the medians are
  significantly different.

- bxp_notchwidth:

  For a notched box plot, width of the notch relative to the body
  (default 0.5).

- hp_alpha:

  Alpha for the jitter.

- hp_col:

  Color for the jitter.

- hp_size:

  Size for the jitter.

- hp_shape:

  Shape for the jitter.

- na.rm:

  If `FALSE`, the default, missing values are removed with a warning. If
  `TRUE`, missing values are silently removed.

## Value

A plot containing jittered points, a boxplot and a density plot or
histogram illustrating the distribution of every group of the data under
evaluation.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 plot_dist(dat, "ETA1", identity_line = T, he_slab_type = "histogram", he_breaks = 30)
} # }
```
