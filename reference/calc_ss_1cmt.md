# Calculate C(t) for a 1-compartment linear model at steady-state

Calculate C(t) for a 1-compartment linear model at steady-state

## Usage

``` r
calc_ss_1cmt(tad, tau, dose, dur = NULL, tinf = NULL, ...)

calc_ss_1cmt_linear_bolus(tad, tau, dose, ...)

calc_ss_1cmt_linear_infusion(tad, tau, dose, tinf, ...)

calc_ss_1cmt_linear_oral_0(tad, tau, dose, dur, ...)

calc_ss_1cmt_linear_oral_0_lag(tad, tau, dose, dur, ...)

calc_ss_1cmt_linear_oral_1_lag(tad, tau, dose, ...)

calc_ss_1cmt_linear_oral_1(tad, tau, dose, ...)
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
  [`calc_derived_1cpt()`](https://kestrel99.github.io/pmxTools/reference/calc_derived.md)

## Value

Concentration of drug at requested time (`t`) after a single dose, given
provided set of parameters and variables.

## Functions

- `calc_ss_1cmt_linear_bolus()`: Calculate C(t) for a 1-compartment
  linear model with IV bolus dosing at steady state

- `calc_ss_1cmt_linear_infusion()`: Calculate C(t) for a 1-compartment
  linear model with infusion at steady state

- `calc_ss_1cmt_linear_oral_0()`: Calculate C(t) for a 1-compartment
  linear model with zero-order oral absorption at steady state

- `calc_ss_1cmt_linear_oral_0_lag()`: Calculate C(t) for a 1-compartment
  linear model with zero-order oral absorption at steady state, with lag
  time

- `calc_ss_1cmt_linear_oral_1_lag()`: Calculate C(t) for a 1-compartment
  linear model with first-order oral absorption at steady state, with
  lag time

- `calc_ss_1cmt_linear_oral_1()`: Calculate C(t) for a 1-compartment
  linear model with first-order oral absorption at steady state

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
Ct <- calc_ss_1cmt_linear_bolus(t=0:24, CL=6, V=25, dose=600, tau=24)
Ct <- calc_ss_1cmt_linear_infusion(tad=0:36, CL=2, V=25, dose=600, tinf=1, tau=24)
Ct <- calc_ss_1cmt_linear_oral_0(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24)
Ct <- calc_ss_1cmt_linear_oral_0_lag(tad=0:36, CL=2, V=25, dose=600, dur=1, tau=24, tlag=1.5)
Ct <- calc_ss_1cmt_linear_oral_1_lag(tad=0:36, CL=2, V=25, dose=600,
    ka=0.25, tlag=0.75, tau=24)
Ct <- calc_ss_1cmt_linear_oral_1(tad=0:36, CL=2, V=25, dose=600, ka=0.25, tau=24)
```
