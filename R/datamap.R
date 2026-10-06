#' Create a data map showing individual dosing and observation records over time
#'
#' @param dat A data frame containing PK/PD data with columns for `id`, `time`, `event`, and optionally `dose` and `group`.
#' @param id Column name for subject ID (defaults to "ID").
#' @param time Column name for time variable (defaults to "TIME").
#' @param event Column name for event type identifier (defaults to "TYPE").
#' @param events A named list defining event types. Each element should have:
#'   \describe{
#'     \item{value}{Numeric value in the event column identifying this event type}
#'     \item{label}{Display label for the legend}
#'     \item{color}{Point color (optional, defaults to grey for obs, red for dose)}
#'     \item{shape}{Point shape (optional, defaults vary by type)}
#'     \item{size}{Point size (optional, defaults vary by type)}
#'     \item{alpha}{Point alpha/transparency (optional, defaults vary by type)}
#'     \item{line}{Logical, whether to draw connecting lines for this event type (optional)}
#'   }
#'   Examples of common event types:
#'   \describe{
#'     \item{dose}{Dose events (TYPE = 0)}
#'     \item{pk}{PK observations (TYPE = 1)}
#'     \item{pd}{PD observations (TYPE = 2)}
#'   }
#'   If NULL, uses default NONMEM EVID events with doses (1) and observations (0).
#' @param group Column name for group identifier (optional, used for faceting).
#' @param ysize Font size for y-axis text (defaults to 6).
#' @param max_time Maximum time to display on y-axis (defaults to 500).
#' @param time_unit Unit for time axis label ("d" for days, "h" for hours, defaults to "d").
#' @param line_color Color for connecting lines (defaults to "#B1B1B1").
#' @param show_lines Logical, whether to show connecting lines (defaults to TRUE).
#'
#' @return A ggplot2 object showing individual data records over time.
#'
#' @examples
#' # Create sample data
#' sample_data <- data.frame(
#'   ID = rep(1:3, each = 4),
#'   TIME = rep(c(0, 1, 12, 24), 3),
#'   TYPE = rep(c(1, 0, 0, 0), 3),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Default: dose and PK observations
#' \dontrun{
#' datamap(sample_data)
#' }
#'
#' # With custom events
#' events <- list(
#'   dose = list(value = 1, label = "Dose", color = "#902C10", shape = 3, size = 1.5),
#'   pk = list(value = 0, label = "PK", color = "#333333", shape = 1, size = 1)
#' )
#' \dontrun{
#' datamap(sample_data, events = events)
#' }
#'
#' @import ggplot2
#' @export
datamap <- function(dat,
                    id = "ID",
                    time = "TIME",
                    event = "TYPE",
                    events = NULL,
                    group = NULL,
                    ysize = 6,
                    max_time = 500,
                    time_unit = "d",
                    line_color = "#B1B1B1",
                    show_lines = TRUE) {

  # Default event definitions
  if (is.null(events)) {
    events <- list(
      dose = list(value = 1, label = "Dose", color = "#902C10",
                  shape = 3, size = 1.5, alpha = 1, line = FALSE),
      pk = list(value = 0, label = "PK", color = "#333333",
                shape = 1, size = 1, alpha = 0.75, line = TRUE)
    )
  }

  # Validate and fill in event defaults
  default_event <- list(
    color = "#666666",
    shape = 1,
    size = 1,
    alpha = 0.75,
    line = TRUE
  )

  for (name in names(events)) {
    events[[name]] <- modifyList(default_event, events[[name]])
    if (is.null(events[[name]]$value)) {
      stop("Event '", name, "' must have a 'value' specified")
    }
    if (is.null(events[[name]]$label)) {
      events[[name]]$label <- name
    }
  }

  # Extract relevant columns
  dat$TIME <- dat[[time]]
  dat$ID <- dat[[id]]
  dat$EVENT <- dat[[event]]

  # Get unique event values in data
  present_events <- unique(dat$EVENT)

  # Filter events to those present in data
  events_present <- events[sapply(events, function(e) e$value %in% present_events)]

  if (length(events_present) == 0) {
    stop("No events found in data. Check 'event' column and 'events' values.")
  }

  # Time formatting
  time_divisor <- if (time_unit == "h") 1 else 24
  time_label <- if (time_unit == "h") "Time (h)" else "Time (d)"
  time_breaks <- if (time_unit == "h") seq(0, max_time * 24, by = 168) else seq(0, max_time, by = 28)

  # Build ggplot
  p <- ggplot(dat, aes(factor(ID), TIME / time_divisor))

  # Add connecting lines if requested and any event has line = TRUE
  if (show_lines) {
    for (evt in events_present) {
      if (isTRUE(evt$line)) {
        p <- p + geom_line(
          data = subset(dat, EVENT == evt$value),
          color = line_color,
          linetype = 1
        )
        break  # Only need one line layer
      }
    }
  }

  # Add point layers for each event type
  colors <- c()
  shapes <- c()
  labels <- c()

  for (evt_name in names(events_present)) {
    evt <- events_present[[evt_name]]
    evt_data <- subset(dat, EVENT == evt$value)

    if (nrow(evt_data) > 0) {
      p <- p + geom_point(
        data = evt_data,
        aes(
          x = factor(ID),
          y = TIME / time_divisor,
          color = !!evt_name,
          shape = !!evt_name
        ),
        size = evt$size,
        alpha = evt$alpha,
        show.legend = TRUE
      )
      colors <- c(colors, evt$color)
      shapes <- c(shapes, evt$shape)
      labels <- c(labels, evt$label)
    }
  }

  names(colors) <- names(events_present)
  names(shapes) <- names(events_present)
  names(labels) <- names(events_present)

  p <- p +
    coord_flip() +
    scale_x_discrete("Subject") +
    scale_y_continuous(time_label, breaks = time_breaks) +
    scale_color_manual(
      name = "Event type",
      values = colors,
      labels = labels,
      breaks = names(events_present)
    ) +
    scale_shape_manual(
      name = "Event type",
      values = shapes,
      labels = labels,
      breaks = names(events_present)
    ) +
    guides(color = guide_legend(override.aes = list(shape = shapes))) +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      panel.grid.major.x = element_line(color = "lightgrey", linetype = 2),
      axis.text.y = element_text(size = ysize)
    )

  if (!is.null(group)) {
    dat$GROUP <- dat[[group]]
    p <- p + facet_wrap(~GROUP, ncol = 1, scales = "free_y")
  }

  p
}
