# Calculate C(t) for a 1-compartment linear model

Calculate C(t) for a 1-compartment linear model

## Usage

``` r
calc_sd_1cmt(t, dose, dur = NULL, tinf = NULL, ...)

calc_sd_1cmt_linear_bolus(t, dose, ...)

calc_sd_1cmt_linear_oral_1_lag(t, dose, ...)

calc_sd_1cmt_linear_infusion(t, dose, tinf, ...)

calc_sd_1cmt_linear_oral_0(t, dose, dur, ...)

calc_sd_1cmt_linear_oral_1(t, dose, ...)

calc_sd_1cmt_linear_oral_0_lag(t, dose, dur, ...)
```

## Arguments

- t:

  Time after dose (h)

- dose:

  Dose

- dur:

  Duration of zero-order absorption (h)

- tinf:

  Duration of infusion (h)

- ...:

  Passed to
  [`calc_derived_1cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)

## Value

Concentration of drug at requested time (`t`) after a single dose, given
provided set of parameters and variables.

## Functions

- `calc_sd_1cmt_linear_bolus()`: Calculate C(t) for a 1-compartment
  linear model after a single IV bolus dose

- `calc_sd_1cmt_linear_oral_1_lag()`: Calculate C(t) for a 1-compartment
  linear model with first-order absorption after a single oral dose,
  with lag time

- `calc_sd_1cmt_linear_infusion()`: Calculate C(t) for a 1-compartment
  linear model after a single IV infusion

- `calc_sd_1cmt_linear_oral_0()`: Calculate C(t) for a 1-compartment
  linear model with zero-order absorption after a single oral dose

- `calc_sd_1cmt_linear_oral_1()`: Calculate C(t) for a 1-compartment
  linear model with first-order absorption after a single oral dose

- `calc_sd_1cmt_linear_oral_0_lag()`: Calculate C(t) for a 1-compartment
  linear model with zero-order absorption after a single oral dose, with
  lag time

## References

Bertrand J & Mentre F (2008). Mathematical Expressions of the
Pharmacokinetic and Pharmacodynamic Models implemented in the Monolix
software.
<https://www.facm.ucl.ac.be/cooperation/Vietnam/WBI-Vietnam-October-2011/Modelling/Monolix32_PKPD_library.pdf>

Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics:
Concepts and Applications (4th). Lippincott Williams & Wilkins,
Philadelphia, 2010.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

Bill Denney, <wdenney@humanpredictions.com>

## Examples

``` r
Ct <- calc_sd_1cmt_linear_bolus(t=0:24, CL=6, V=25, dose=600)
Ct <- calc_sd_1cmt_linear_oral_1_lag(t=0:24, CL=6, V=25, ka=1.1, dose=600, tlag=2)
Ct <- calc_sd_1cmt_linear_infusion(t=0:24, CL=6, V=25, dose=600, tinf=1)
Ct <- calc_sd_1cmt_linear_oral_0(t=0:24, CL=6, V=25, dur=1.5, dose=600)
Ct <- calc_sd_1cmt_linear_oral_1(t=0:24, CL=6, V=25, ka=1.1, dose=600)
Ct <- calc_sd_1cmt_linear_oral_0_lag(t=0:24, CL=6, V=25, dur=1.5, dose=600, tlag=1.5)
```
