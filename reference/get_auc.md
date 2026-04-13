# Calculate the area under the curve (AUC) for each subject over the time interval for dependent variables (`dv`) using the trapezoidal rule.

Calculate the area under the curve (AUC) for each subject over the time
interval for dependent variables (`dv`) using the trapezoidal rule.

## Usage

``` r
get_auc(data, time = "TIME", id = "ID", dv = "DV")
```

## Arguments

- data:

  A data frame.

- time:

  A string containing the name of the chronologically ordered time
  variable in `data`.

- id:

  A string containing the name of the ID column (defining subject level
  data) in `data`.

- dv:

  A string containing the name of the dependent variable column in
  `data`.

## Value

A data frame containing one AUC value for every subject as defined by
`id`.

Based on the `AUC` function originally written by Leonid Gibiansky in
package MIfuns 5.1, from Metrum Institute.

## References

<https://code.google.com/archive/p/mifuns/>

## Author

Leonid Gibiansky, <lgibiansky@quantpharm.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 AUCs <- get_auc(myAUCdata)
} # }
```
