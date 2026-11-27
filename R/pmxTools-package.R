#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

# Only the ggplot2 functions used are imported, to avoid masking other
# packages' exports (`.data` and `sym` are ggplot2's re-exports from rlang).
#' @importFrom ggplot2 .data aes after_stat coord_cartesian coord_flip
#'   element_blank element_line element_text facet_grid facet_wrap geom_boxplot
#'   geom_hline geom_jitter geom_line geom_point ggplot guide_legend guides labs
#'   position_nudge scale_color_manual scale_shape_manual scale_x_continuous
#'   scale_x_discrete scale_y_continuous sym theme theme_bw theme_light
NULL

utils::globalVariables("pdf")
utils::globalVariables("cdf")

NULL
