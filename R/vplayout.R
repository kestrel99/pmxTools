#' Helper function for arranging plots on a page.
#'
#' @param x Row.
#' @param y Column.
#'
#' @return A layout for the viewport.
#'
#' @seealso \code{\link{grid::grid.newpage}, \code{\link{grid::viewport}, \code{\link{grid::pushViewport}, \code{\link{print}}.
#'
#' @examples
#' library(ggplot)
#' df <- data.frame(gp = factor(rep(letters[1:3], each = 10)), y = rnorm(30))
#' plot1 <- ggplot(df, aes(gp, y)) + geom_point()
#' plot2 <- ggplot(df, aes(y, gp)) + geom_point()
#'
#' grid.newpage()
#' pushViewport(viewport(layout = grid.layout(1, 2)))
#' print(plot1, vp = vplayout(1, 1))
#' print(plot2, vp = vplayout(1, 2))

vplayout <- function(x, y) {
  viewport(layout.pos.row=x, layout.pos.col=y)
}
