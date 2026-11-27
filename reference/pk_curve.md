# Provide concentration-time curves.

Concentrations after repeated dosing are obtained by superposition: each
dose given at time \\n \cdot ii\\ (\\n = 0, \ldots, addl\\) contributes
the single-dose concentration at the time elapsed since that dose, and
doses not yet given contribute nothing. Observation times need not be
regularly spaced or start at zero.

## Usage

``` r
pk_curve(
  t,
  model = "1cmt_oral",
  params = list(ka = 2.77, CL = 2.5, V = 25),
  dose = 600,
  ii = 24,
  addl = 0,
  ss = FALSE
)
```

## Arguments

- t:

  Observation time in h, specified as a vector of non-negative values
  relative to the first dose.

- model:

  The model to use. Must be one of "1cmt_bolus", "1cmt_infusion",
  "1cmt_oral", "2cmt_bolus", "2cmt_infusion", "2cmt_oral", "3cmt_bolus",
  "3cmt_infusion", "3cmt_oral". The default is "1cmt_oral".

- params:

  A named list containing parameter values for the selected model type.

- dose:

  Dose amount.

- ii:

  Interdose interval (or tau), in hours (default 24).

- addl:

  Number of additional doses (default 0).

- ss:

  Assume steady state concentration (default `FALSE`). When `TRUE`,
  steady state is assumed to have been reached with dosing interval `ii`
  before the dose at time 0 (as for `SS=1` in NONMEM), and the `addl`
  additional doses are given on top of it. After the last dose the curve
  declines as the drug washes out.

## Value

A data frame containing times (`t`) and concentrations (`cp`).

## Author

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
plot(pk_curve(t=seq(0,72,by=0.1), model="3cmt_oral", ii=12, addl=5,
  params=list(CL=2.5, V1=25, V2=2, V3=5, Q2=0.5, Q3=0.25, ka=1)), type="l")


# Steady state, followed by washout
plot(pk_curve(t=seq(0,72,by=0.1), model="1cmt_oral", ii=12, ss=TRUE,
  params=list(CL=2.5, V=25, ka=1)), type="l")
```
