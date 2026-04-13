# Generate a summary table of descriptive data for every individual in a dataset suitable for tabulation in a report.

Generate a summary table of descriptive data for every individual in a
dataset suitable for tabulation in a report.

## Usage

``` r
dgr_table(
  dat,
  fields,
  names = NULL,
  cutoff = 7,
  sig = 3,
  by = NULL,
  idvar = "ID",
  navars = c("-99", "-999"),
  mtype = "geomean"
)
```

## Arguments

- dat:

  An input data frame, with one row per unique individual.

- fields:

  A named character vector where names are column names and values are
  descriptive labels. For backward compatibility, can also be an unnamed
  character vector with separate `names` argument.

- names:

  Deprecated. Use a named vector for `fields` instead.

- cutoff:

  An integer defining the maximum number of unique values a variable
  should have to be considered categorical. Fields with more than this
  number of unique values are considered continuous for the purposes of
  the summary table (defaults to 7).

- sig:

  The number of significant digits summary values should have (defaults
  to 3).

- by:

  The field to use for grouping (a string). If not `NULL` (the default),
  the summary table will contain columns for each unique value of this
  field, as well as a column summarizing across all fields.

- idvar:

  The field in the dataset identifying each unique individual (defaults
  to "ID").

- navars:

  A vector containing values that are to be interpreted as missing
  (defaults to "-99" and "-999"). `NA` values are always considered to
  be missing.

- mtype:

  The type of mean to apply; `geomean`, the geometric mean (default) or
  `mean`, the arithmetic mean.

## Value

A data frame containing a summary of all the fields listed in `fields`,
for each individual in the dataset (the dataset should not contain
duplicated individuals), conditioned on the field in `by`. Continuous
values are summarized as median, mean, range and number of missing
values. Categorical values are summarized as count and relative
percentage.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
dat <- data.frame(
  ID = 1:10,
  AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
  SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F")
)

dgr_table(dat, c(AGE = "Age (years)", SEX = "Sex"))
#> # A tibble: 5 × 2
#>   Variable    Total                   
#>   <chr>       <chr>                   
#> 1 N           "10"                    
#> 2 Age (years) "46.5 (45.7)\n[29 ; 67]"
#> 3 Sex         ""                      
#> 4 - F         "5 (50%)"               
#> 5 - M         "5 (50%)"               
dgr_table(dat, c(AGE = "Age (years)"), by = "SEX")
#> # A tibble: 2 × 4
#>   Variable    F                      M                      Total               
#>   <chr>       <chr>                  <chr>                  <chr>               
#> 1 N           "5"                    "5"                    "10"                
#> 2 Age (years) "52 (48.4)\n[29 ; 67]" "43 (43.1)\n[33 ; 61]" "46.5 (45.7)\n[29 ;…
```
