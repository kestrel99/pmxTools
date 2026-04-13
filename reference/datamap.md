# Create a data map showing individual dosing and observation records over time

Create a data map showing individual dosing and observation records over
time

## Usage

``` r
datamap(
  dat,
  id = "ID",
  time = "TIME",
  event = "TYPE",
  events = NULL,
  group = NULL,
  ysize = 6,
  max_time = 500,
  time_unit = "d",
  line_color = "#B1B1B1",
  show_lines = TRUE
)
```

## Arguments

- dat:

  A data frame containing PK/PD data with columns for `id`, `time`,
  `event`, and optionally `dose` and `group`.

- id:

  Column name for subject ID (defaults to "ID").

- time:

  Column name for time variable (defaults to "TIME").

- event:

  Column name for event type identifier (defaults to "TYPE").

- events:

  A named list defining event types. Each element should have:

  value

  :   Numeric value in the event column identifying this event type

  label

  :   Display label for the legend

  color

  :   Point color (optional, defaults to grey for obs, red for dose)

  shape

  :   Point shape (optional, defaults vary by type)

  size

  :   Point size (optional, defaults vary by type)

  alpha

  :   Point alpha/transparency (optional, defaults vary by type)

  line

  :   Logical, whether to draw connecting lines for this event type
      (optional)

  Examples of common event types:

  dose

  :   Dose events (TYPE = 0)

  pk

  :   PK observations (TYPE = 1)

  pd

  :   PD observations (TYPE = 2)

  If NULL, uses default NONMEM EVID events with doses (1) and
  observations (0).

- group:

  Column name for group identifier (optional, used for faceting).

- ysize:

  Font size for y-axis text (defaults to 6).

- max_time:

  Maximum time to display on y-axis (defaults to 500).

- time_unit:

  Unit for time axis label ("d" for days, "h" for hours, defaults to
  "d").

- line_color:

  Color for connecting lines (defaults to "#B1B1B1").

- show_lines:

  Logical, whether to show connecting lines (defaults to TRUE).

## Value

A ggplot2 object showing individual data records over time.

## Examples

``` r
# Create sample data
sample_data <- data.frame(
  ID = rep(1:3, each = 4),
  TIME = rep(c(0, 1, 12, 24), 3),
  TYPE = rep(c(1, 0, 0, 0), 3),
  stringsAsFactors = FALSE
)

# Default: dose and PK observations
if (FALSE) { # \dontrun{
datamap(sample_data)
} # }

# With custom events
events <- list(
  dose = list(value = 1, label = "Dose", color = "#902C10", shape = 3, size = 1.5),
  pk = list(value = 0, label = "PK", color = "#333333", shape = 1, size = 1)
)
if (FALSE) { # \dontrun{
datamap(sample_data, events = events)
} # }
```
