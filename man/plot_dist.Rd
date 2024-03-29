% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dist.R
\name{plot_dist}
\alias{plot_dist}
\title{Plot a distribution as a hybrid containing a halfeye, a boxplot and jittered points.}
\usage{
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
  hp_range_scale = 0.4,
  hp_alpha = 0.25,
  hp_col = "#F8766D",
  hp_transformation = position_jitter(),
  na.rm = FALSE
)
}
\arguments{
\item{dat}{A data frame.}

\item{yvar}{The name of the field containing values to be plotted.}

\item{xvar}{The name of the field containing the grouping variable (defaults to `NULL`).}

\item{ylim}{Limits for the y-axis. Defaults to \code{NULL}. If provided, should be a 2-element vector containing the upper and lower limits.}

\item{xlb}{Label for the x-axis.}

\item{ylb}{Label for the y-axis.}

\item{identity_line}{Show a line of identity? Default \code{FALSE}.}

\item{identity_value}{If an identity line is shown, it will be drawn horizontally at this y-value (default 0).}

\item{he_adjust}{If \code{he_slab_type} is \code{"pdf"}, bandwidth for the density estimator is adjusted by multiplying it by this value.}

\item{he_width}{Width of the halfeye component of the plot (default 0.6).}

\item{he_justification}{Justification of the halfeye component of the plot (default -0.2).}

\item{he_col}{Color for the halfeye component of the plot.}

\item{he_fill}{Fill color for the halfeye component of the plot.}

\item{he_alpha}{Alpha for the halfeye component of the plot (default 0.9).}

\item{he_slab_type}{The type of slab function to calculate for the halfeye component of the plot: probability density (or mass) function (\code{"pdf"}, the default), cumulative distribution function (\code{"cdf"}), complementary CDF (\code{"ccdf"}) or histogram (\code{"histogram"}).}

\item{he_breaks}{If slab_type is \code{"histogram"}, the breaks parameter that is passed to \code{hist()} to determine where to put breaks in the histogram.}

\item{he_outline_bars}{If slab_type is \code{"histogram"}, determines if outlines in between the bars are drawn when the slab_color aesthetic is used. If \code{FALSE} (the default), the outline is drawn only along the tops of the bars; if \code{TRUE}, outlines in between bars are also drawn.}

\item{he_point_interval}{A function from the \code{\link[ggdist:point_interval]{ggdist::point_interval}} family (e.g., \code{median_qi}, \code{mean_qi}, \code{mode_hdi}, etc), or a string giving the name of a function from that family (e.g., \code{"median_qi"}, \code{"mean_qi"}, \code{"mode_hdi"}, etc. This function determines the point summary (typically mean, median, or mode) and interval type (quantile interval, \code{qi}; highest-density interval, \code{hdi}; or highest-density continuous interval, \code{hdci}). Output will be converted to the appropriate x- or y-based aesthetics depending on the value of orientation.}

\item{bxp_width}{Width of the boxplot component (default 0.12).}

\item{bxp_outlier_col}{Color for outliers in the boxplot component.}

\item{bxp_outlier_fill}{Fill color for outliers in the boxplot component.}

\item{bxp_outlier_shape}{Shape for outliers in the boxplot component.}

\item{bxp_outlier_size}{Size for outliers in the boxplot component.}

\item{bxp_col}{Color for the boxplot component.}

\item{bxp_fill}{Fill color for the boxplot component.}

\item{bxp_alpha}{Alpha for the boxplot component.}

\item{bxp_notch}{If \code{FALSE} (default) make a standard box plot. If \code{TRUE}, make a notched box plot. Notches are used to compare groups; if the notches of two boxes do not overlap, this suggests that the medians are significantly different.}

\item{bxp_notchwidth}{For a notched box plot, width of the notch relative to the body (default 0.5).}

\item{hp_range_scale}{If no `width` argument is specified in \code{hp_transformation}, used to determine the width of the jitter. Defaults to 0.75, which is half of the allotted space for the jitter-points, whereas 1 would use all of the alotted space.}

\item{hp_alpha}{Alpha for the jitter.}

\item{hp_col}{Color for the jitter.}

\item{hp_transformation}{An evaluated \code{position_*()} function yielding a `Position` object with specified parameters to calculate the transformation of the points. Defaults to \code{\link[ggplot2:position_jitter]{ggplot2::position_jitter}}.}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with a warning. If \code{TRUE}, missing values are silently removed.}
}
\value{
A plot containing jittered points, a boxplot and a density plot or histogram illustrating the distribution of every group of the data under evaluation.
}
\description{
Plot a distribution as a hybrid containing a halfeye, a boxplot and jittered points.
}
\examples{
\dontrun{
 plot_dist(dat, "ETA1", identity_line = T, he_slab_type = "histogram", he_breaks = 30)
}

}
\author{
Justin Wilkins, \email{justin.wilkins@occams.com}
}
