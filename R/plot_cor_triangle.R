#' Plot a correlation matrix as a lower triangle
#'
#' Draws the lower triangle of a correlation matrix as a heatmap, with blues
#' for negative and oranges/reds for positive correlations, each cell labelled
#' with its value.
#'
#' @param R A square correlation matrix with column names, such as the output
#'   of [het_cor()].
#' @param types Optional matrix of correlation types (`"Pearson"`,
#'   `"Polyserial"`, `"Polychoric"`) with the same dimensions as `R`. When
#'   given, polyserial cells are marked `*` and polychoric cells `**`, and a
#'   footnote explains the marks. [het_cor()] supplies this when
#'   `plot = TRUE`.
#' @param labels Optional display labels, named by variable, e.g.
#'   `c(WT = "Body weight (kg)")`. Unlabelled variables keep their names.
#' @param digits Number of decimal places shown in each cell.
#' @param text_size Font size of the cell values, in points.
#' @return A ggplot object.
#' @seealso [het_cor()]
#' @examples
#' R <- cor(mtcars[c("mpg", "disp", "hp", "wt")])
#' plot_cor_triangle(R, labels = c(wt = "Weight (1000 lb)"))
#' @export
plot_cor_triangle <- function(R, types = NULL, labels = character(),
                              digits = 2, text_size = 11) {
  if (!is.matrix(R) || nrow(R) != ncol(R) || is.null(colnames(R))) {
    stop("`R` must be a square matrix with column names", call. = FALSE)
  }
  vars <- colnames(R)
  bad_labels <- is.null(names(labels)) || !all(names(labels) %in% vars)
  if (length(labels) && bad_labels) {
    stop("`labels` must be named by variables in R", call. = FALSE)
  }
  relabel <- function(v) unname(ifelse(is.na(labels[v]), v, labels[v]))

  df <- as.data.frame(as.table(R), stringsAsFactors = FALSE)
  names(df) <- c("row", "col", "r")
  df <- df[match(df$row, vars) > match(df$col, vars), ]

  marks <- c(Pearson = "", Polyserial = "*", Polychoric = "**")
  notes <- c(Pearson = "Unmarked: Pearson", Polyserial = "* Polyserial",
             Polychoric = "** Polychoric")
  df$type <- if (is.null(types)) "Pearson" else types[cbind(df$row, df$col)]
  caption <- if (is.null(types)) NULL else
    paste(notes[intersect(names(notes), df$type)], collapse = "    ")

  df$row <- factor(df$row, levels = rev(vars))
  df$col <- factor(df$col, levels = vars)
  value <- formatC(df$r, digits = digits, format = "f")
  df$label <- ifelse(is.na(df$r), "", paste0(value, marks[df$type]))
  df$dark <- !is.na(df$r) & abs(df$r) > 0.6

  ggplot2::ggplot(df, ggplot2::aes(.data$col, .data$row, fill = .data$r)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$dark),
                       size = text_size / ggplot2::.pt) +
    ggplot2::scale_fill_gradientn(
      colours = c("#08306B", "#6BAED6", "#F7F7F7", "#FD8D3C", "#B10026"),
      limits = c(-1, 1), na.value = "grey85", name = NULL
    ) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"),
                                 guide = "none") +
    ggplot2::scale_x_discrete(labels = relabel) +
    ggplot2::scale_y_discrete(labels = relabel) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, caption = caption) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.caption.position = "plot")
}
