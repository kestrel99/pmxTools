# Calculate C(t) for a 3-compartment linear model

Calculate C(t) for a 3-compartment linear model

## Usage

``` r
calc_sd_3cmt(t, dose, dur = NULL, tinf = NULL, ...)

calc_sd_3cmt_linear_bolus(t, dose, ...)

calc_sd_3cmt_linear_oral_1_lag(t, dose, ...)

calc_sd_3cmt_linear_infusion(t, dose, tinf, ...)

calc_sd_3cmt_linear_oral_0(t, dose, dur, ...)

calc_sd_3cmt_linear_oral_0_lag(t, dose, dur, ...)

calc_sd_3cmt_linear_oral_1(t, dose, ...)
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
  [`calc_derived_3cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)

## Value

Concentration of drug at requested time (`t`) after a single dose, given
provided set of parameters and variables.

## Functions

- `calc_sd_3cmt_linear_bolus()`: Calculate C(t) for a 3-compartment
  linear model after a single IV bolus dose

- `calc_sd_3cmt_linear_oral_1_lag()`: Calculate C(t) for a 3-compartment
  linear model after a single oral dose

- `calc_sd_3cmt_linear_infusion()`: Calculate C(t) for a 3-compartment
  linear model after a single IV infusion

- `calc_sd_3cmt_linear_oral_0()`: Calculate C(t) for a 3-compartment
  linear model after a single dose, with zero-order absorption

- `calc_sd_3cmt_linear_oral_0_lag()`: Calculate C(t) for a 3-compartment
  linear model after a single dose, with zero-order absorption and a lag
  time

- `calc_sd_3cmt_linear_oral_1()`: Calculate C(t) for a 3-compartment
  linear model after a single oral dose

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
Ct <- calc_sd_3cmt_linear_bolus(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dose = 100)
Ct <- calc_sd_3cmt_linear_oral_1_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tlag = 1.5)
Ct <- calc_sd_3cmt_linear_infusion(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dose = 100, tinf=1)
Ct <- calc_sd_3cmt_linear_oral_0(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100)
Ct <- calc_sd_3cmt_linear_oral_0_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tlag=1.5)
Ct <- calc_sd_3cmt_linear_oral_1(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100)
```
