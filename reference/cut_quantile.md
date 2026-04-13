# Create quantile-based bins for continuous variables

Create quantile-based bins for continuous variables

## Usage

``` r
cut_quantile(
  dat,
  var,
  n_groups = 4,
  missing_codes = c(-99, -999),
  blq_label = "BLQ",
  unit = NULL,
  id = NULL,
  verbose = FALSE
)
```

## Arguments

- dat:

  A data frame containing the variables to bin.

- var:

  Variable(s) to bin: single name, character vector, or named list.
  Named list allows different quantile cuts per variable, e.g.
  `list(AGE = c(4, 3), WT = 4)` creates both quartiles and tertiles for
  AGE.

- n_groups:

  Number of quantile groups (2-5). Can be:

  - Single value applied to all variables

  - Vector applied in order

  - Named list for per-variable settings

- missing_codes:

  Values to treat as missing (defaults to c(-99, -999)).

- blq_label:

  Label for BLQ/zero values (defaults to "BLQ").

- unit:

  Unit for display, appended to interval labels. Can be single value or
  named list per variable.

- id:

  Optional subject ID column. If provided, quantile calculation uses 1
  row per subject. Stops if a subject has conflicting values.

- verbose:

  If TRUE, prints summary and returns list with data and summary.

## Value

If `verbose = FALSE` (default): returns modified data frame invisibly.
If `verbose = TRUE`: returns a list with:

- `data`: modified data frame

- `summary`: tibble with cut details per variable

- `skipped`: tibble of any skipped cuts (due to zero-range bins)

Output columns added (for var = "CONC", n_groups = 4):

- `CONCQ4Q`: numeric factor (1, 2, 3, 4)

- `CONCQ4C`: character ("Q1", "Q2", "Q3", "Q4", "BLQ")

- `CONCQ4CC`: continuous factor with intervals

## Examples

``` r
if (FALSE) { # \dontrun{
# Single variable, quartiles
dat <- cut_quantile(dat, "AGE", n_groups = 4)

# Multiple cuts on same variable
dat <- cut_quantile(dat, "AGE", n_groups = c(4, 3))

# With longitudinal data
dat <- cut_quantile(dat, list(CONC = 4, AGE = 4), id = "SUBJID")

# Verbose output
result <- cut_quantile(dat, list(CONC = c(4, 3), AGE = 4), verbose = TRUE)
} # }
```
