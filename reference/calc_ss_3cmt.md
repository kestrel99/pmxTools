# Calculate C(t) for a 3-compartment linear model at steady-state

Calculate C(t) for a 3-compartment linear model at steady-state

## Usage

``` r
calc_ss_3cmt(tad, tau, dose, dur = NULL, tinf = NULL, ...)

calc_ss_3cmt_linear_bolus(tad, tau, dose, ...)

calc_ss_3cmt_linear_oral_1_lag(tad, tau, dose, ...)

calc_ss_3cmt_linear_infusion(tad, tau, dose, tinf, ...)

calc_ss_3cmt_linear_oral_0(tad, tau, dose, dur, ...)

calc_ss_3cmt_linear_oral_0_lag(tad, tau, dose, dur, ...)

calc_ss_3cmt_linear_oral_1(tad, tau, dose, ...)
```

## Arguments

- tad:

  Time after dose (h)

- tau:

  Dosing interval (h)

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

Concentration of drug at requested time (`t`) at steady-state, given
provided set of parameters and variables.

## Functions

- `calc_ss_3cmt_linear_bolus()`: Calculate C(t) for a 3-compartment
  linear model at steady state with IV bolus dosing

- `calc_ss_3cmt_linear_oral_1_lag()`: Calculate C(t) for a 3-compartment
  linear model at steady-state with first-order oral dosing with a lag
  time

- `calc_ss_3cmt_linear_infusion()`: Calculate C(t) for a 3-compartment
  linear model at steady state with IV infusions

- `calc_ss_3cmt_linear_oral_0()`: Calculate C(t) for a 3-compartment
  linear model at steady state, with zero-order absorption

- `calc_ss_3cmt_linear_oral_0_lag()`: Calculate C(t) for a 3-compartment
  linear model at steady state, with zero-order absorption and lag time

- `calc_ss_3cmt_linear_oral_1()`: Calculate C(t) for a 3-compartment
  linear model at steady-state with first-order oral dosing

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

## Examples

``` r
Ct <- calc_ss_3cmt_linear_bolus(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dose = 100, tau=24)
Ctrough <- calc_ss_3cmt_linear_oral_1_lag(t = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau=24, tlag = 1.5)
Ct <- calc_ss_3cmt_linear_infusion(tad = 11.75, CL = 2.5, V1 = 20, V2 = 50,
    V3 = 100, Q2 = 0.5, Q3 = 0.05, dose = 1000, tinf=1, tau=24)
Ct <- calc_ss_3cmt_linear_oral_0(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24)
Ct <- calc_ss_3cmt_linear_oral_0_lag(tad = 11.75, CL = 3.5, V1 = 20, V2 = 500,
    V3 = 200, Q2 = 0.5, Q3 = 0.05, dur = 1, dose = 100, tau = 24, tlag = 1.5)
Ct <- calc_ss_3cmt_linear_oral_1(tad = 11.75, CL = 3.5, V1 = 20,
    V2 = 500, V3 = 200, Q2 = 0.5, Q3 = 0.05, ka = 1, dose = 100, tau = 24)
```
