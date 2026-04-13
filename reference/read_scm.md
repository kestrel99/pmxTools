# Read PsN SCM output into a format suitable for further use.

`read_scm` returns a summary of a Perl-speaks-NONMEM (PsN,
<https://uupharmacometrics.github.io/PsN/>) SCM (stepwise covariate
modeling) procedure. It depends on the presence of `scmlog.txt` and
`short_scmlog.txt` files in the specified directory.

## Usage

``` r
read_scm(dir, startPhase = "forward")
```

## Arguments

- dir:

  A PsN SCM folder (containing `scmlog.txt` and `short_scmlog.txt`).

- startPhase:

  Where to start collating the output; can be `"forward"` (the default)
  or `"backward"`.

## Value

A list of data frames, containing

- forward:

  all models evaluated during the forward inclusion step of covariate
  model building

- forwardSummary:

  the covariate relationships selected at each forward step

- forwardP:

  the P-value used for inclusion during the forward inclusion step

- backward:

  all models evaluated during the backward elimination step of covariate
  model building

- backwardSummary:

  the covariate relationships eliminated at each backward step

- backwardP:

  the P-value used for exclusion during the backward elimination step

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

Lindbom L, Ribbing J & Jonsson EN (2004). Perl-speaks-NONMEM (PsN) - A
Perl module for NONMEM related programming. Computer Methods and
Programs in Biomedicine, 75(2), 85-94.
[doi:10.1016/j.cmpb.2003.11.003](https://doi.org/10.1016/j.cmpb.2003.11.003)

Lindbom L, Pihlgren P & Jonsson N (2005). PsN-Toolkit - A collection of
computer intensive statistical methods for non-linear mixed effect
modeling using NONMEM. Computer Methods and Programs in Biomedicine,
79(3), 241-257.
[doi:10.1016/j.cmpb.2005.04.005](https://doi.org/10.1016/j.cmpb.2005.04.005)

Other NONMEM reading:
[`plot_scm()`](https://kestrel99.github.io/pmxTools/reference/plot_scm.md),
[`read_nm()`](https://kestrel99.github.io/pmxTools/reference/read_nm.md),
[`read_nm_all()`](https://kestrel99.github.io/pmxTools/reference/read_nm_all.md),
[`read_nm_multi_table()`](https://kestrel99.github.io/pmxTools/reference/read_nm_multi_table.md),
[`read_nmcov()`](https://kestrel99.github.io/pmxTools/reference/read_nmcov.md),
[`read_nmext()`](https://kestrel99.github.io/pmxTools/reference/read_nmext.md),
[`read_nmtables()`](https://kestrel99.github.io/pmxTools/reference/read_nmtables.md)

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
scm <- read_scm("E:/DrugX/ModelDevelopment/scm310")
} # }
```
