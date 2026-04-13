# Calculate C(t) for a 2-compartment linear model

Calculate C(t) for a 2-compartment linear model

## Usage

``` r
calc_sd_2cmt(t, dose, dur = NULL, tinf = NULL, ...)

calc_sd_2cmt_linear_bolus(t, dose, ...)

calc_sd_2cmt_linear_oral_1_lag(t, dose, ...)

calc_sd_2cmt_linear_infusion(t, dose, tinf, ...)

calc_sd_2cmt_linear_oral_0_lag(t, dose, dur, ...)

calc_sd_2cmt_linear_oral_1(t, dose, ...)

calc_sd_2cmt_linear_oral_0(t, dose, dur, ...)
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
  [`calc_derived_2cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)

## Value

Concentration of drug at requested time (`t`) after a single dose, given
provided set of parameters and variables.

## Functions

- `calc_sd_2cmt_linear_bolus()`: Calculate C(t) for a 2-compartment
  linear model after a single IV bolus dose

- `calc_sd_2cmt_linear_oral_1_lag()`: Calculate C(t) for a 2-compartment
  linear model after a single first-order oral dose with a lag time

- `calc_sd_2cmt_linear_infusion()`: Calculate C(t) for a 2-compartment
  linear model after a single infusion

- `calc_sd_2cmt_linear_oral_0_lag()`: Calculate C(t) for a 2-compartment
  linear model after a single zero-order oral dose, with lag time

- `calc_sd_2cmt_linear_oral_1()`: Calculate C(t) for a 2-compartment
  linear model after a single first-order oral dose

- `calc_sd_2cmt_linear_oral_0()`: Calculate C(t) for a 2-compartment
  linear model after a single zero-order oral dose

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
Ct <- calc_sd_2cmt_linear_bolus(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5, dose = 10)
Ct <- calc_sd_2cmt_linear_oral_1_lag(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
    dose = 1000, ka = 1, tlag = 2)
Ctrough <- calc_sd_2cmt_linear_infusion(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
    dose = 10, tinf = 1)
Ctrough <- calc_sd_2cmt_linear_oral_0_lag(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
    dose = 1000, dur = 1, tlag=2)
Ct <- calc_sd_2cmt_linear_oral_1(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
    dose = 1000, ka = 1)
Ct <- calc_sd_2cmt_linear_oral_0(t = 11.75, CL = 7.5, V1 = 20, V2 = 30, Q = 0.5,
    dose = 1000, dur = 1)
```
