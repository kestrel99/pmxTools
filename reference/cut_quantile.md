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
set.seed(1)
dat <- data.frame(
  SUBJID = rep(1:20, each = 3),
  AGE = rep(round(runif(20, 20, 80)), each = 3),
  CONC = c(0, round(rlnorm(59, 2, 1), 2))
)

# Single variable, quartiles
head(cut_quantile(dat, "AGE", n_groups = 4))
#>   SUBJID AGE  CONC AGEQ4Q AGEQ4C   AGEQ4CC
#> 1      1  36  0.00      1     Q1 [24,40.5)
#> 2      1  36 33.51      1     Q1 [24,40.5)
#> 3      1  36 10.91      1     Q1 [24,40.5)
#> 4      2  42  3.97      2     Q2 [40.5,56)
#> 5      2  42  0.81      2     Q2 [40.5,56)
#> 6      2  42 22.76      2     Q2 [40.5,56)

# Multiple cuts on same variable
head(cut_quantile(dat, "AGE", n_groups = c(4, 3)))
#>   SUBJID AGE  CONC AGEQ4Q AGEQ4C   AGEQ4CC AGET3Q AGET3C AGET3CC
#> 1      1  36  0.00      1     Q1 [24,40.5)      1     T1 [24,43)
#> 2      1  36 33.51      1     Q1 [24,40.5)      1     T1 [24,43)
#> 3      1  36 10.91      1     Q1 [24,40.5)      1     T1 [24,43)
#> 4      2  42  3.97      2     Q2 [40.5,56)      1     T1 [24,43)
#> 5      2  42  0.81      2     Q2 [40.5,56)      1     T1 [24,43)
#> 6      2  42 22.76      2     Q2 [40.5,56)      1     T1 [24,43)

# With longitudinal data, quantiles use one row per subject
head(cut_quantile(dat, "AGE", n_groups = 4, id = "SUBJID"))
#>   SUBJID AGE  CONC AGEQ4Q AGEQ4C   AGEQ4CC
#> 1      1  36  0.00      1     Q1 [24,40.5)
#> 2      1  36 33.51      1     Q1 [24,40.5)
#> 3      1  36 10.91      1     Q1 [24,40.5)
#> 4      2  42  3.97      2     Q2 [40.5,56)
#> 5      2  42  0.81      2     Q2 [40.5,56)
#> 6      2  42 22.76      2     Q2 [40.5,56)

# Verbose output
result <- cut_quantile(dat, list(CONC = c(4, 3), AGE = 4), verbose = TRUE)
#> 
#> ============================================================
#> cut_quantile summary
#> ============================================================
#> # A tibble: 3 × 8
#>   cuts  var   n_groups n_total n_valid n_missing n_blq bins            
#>   <chr> <chr>    <int>   <int>   <int>     <int> <int> <list>          
#> 1 Q4    CONC         4      60      59         0     1 <tibble [5 × 2]>
#> 2 T3    CONC         3      60      59         0     1 <tibble [4 × 2]>
#> 3 Q4    AGE          4      60      60         0     0 <tibble [5 × 2]>
#> ------------------------------------------------------------
```
