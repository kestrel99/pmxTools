# Calculate derived pharmacokinetic parameters for a 1-, 2-, or 3-compartment linear model.

Calculate derived pharmacokinetic parameters for a 1-, 2-, or
3-compartment linear model.

## Usage

``` r
calc_derived(..., verbose = FALSE)

calc_derived_1cpt(
  CL,
  V = NULL,
  V1 = NULL,
  ka = NULL,
  dur = NULL,
  tlag = NULL,
  tinf = NULL,
  dose = NULL,
  tau = NULL,
  step = 0.1,
  type = "all",
  sigdig = 5
)

calc_derived_2cpt(
  CL,
  V1 = NULL,
  V2,
  Q2 = NULL,
  V = NULL,
  Q = NULL,
  dur = NULL,
  tinf = NULL,
  ka = NULL,
  tlag = NULL,
  dose = NULL,
  tau = NULL,
  step = 0.1,
  type = "all",
  sigdig = 5
)

calc_derived_3cpt(
  CL,
  V1 = NULL,
  V2,
  V3,
  Q2 = NULL,
  Q3,
  V = NULL,
  Q = NULL,
  ka = NULL,
  dur = NULL,
  tinf = NULL,
  tlag = NULL,
  dose = NULL,
  tau = NULL,
  step = 0.1,
  type = "all",
  sigdig = 5
)
```

## Arguments

- ...:

  Passed to the other `calc_derived_*()` functions.

- verbose:

  For `calc_derived()`, provide a message indicating the type of model
  detected.

- CL:

  Clearance (volume per time units, e.g. L/h)

- V1, V:

  Central volume of distribution (volume units, e.g. L). Values are
  synonyms; use only one.

- ka:

  Absorption rate (inverse time units, e.g. 1/h)

- dur:

  Duration of zero-order absorption (time units, e.g. h)

- tlag:

  Absorption lag time (time units, e.g. h)

- tinf:

  Duration of infusion (time units, e.g. h)

- dose:

  Dose (amount units, e.g. mg)

- tau:

  Duration of interdose interval (time units, e.g. h; defaults to 24)

- step:

  Time increment to use when estimating NCA parameters (time units, e.g.
  h; defaults to 0.1)

- type:

  Parameters to return. Default is `"all"`. If not `"all"`, this may be
  a vector from the names of the return value list.

- sigdig:

  Number of significant digits to be returned. Default is `5`.

- V2:

  First peripheral volume of distribution (volume units, e.g. L)

- Q2, Q:

  Intercompartmental clearance from central to first peripheral
  compartment (volume per time units, e.g. L/h). Values are synonyms;
  use only one.

- V3:

  Second peripheral volume of distribution (volume units, e.g. L)

- Q3:

  Intercompartmental clearance from central to second peripheral
  compartment (volume per time units, e.g. L/h)

## Value

Return a list of derived PK parameters for a 1-, 2-, or 3-compartment
linear model given provided clearances and volumes based on the `type`.
If a dose is provided, estimated non-compartmental analysis (NCA)
parameters will be provided as well, based on simulation of single-dose
and (if `tau` is specified) steady-state time courses.

- `Vss`: Volume of distribution at steady state, \\V\_{ss}\\ (volume
  units, e.g. L); 1-, 2-, and 3-compartment

- `k10`: First-order elimination rate, \\k\_{10}\\ (inverse time units,
  e.g. 1/h); 1-, 2-, and 3-compartment

- `k12`: First-order rate of transfer from central to first peripheral
  compartment, \\k\_{12}\\ (inverse time units, e.g. 1/h); 2- and
  3-compartment

- `k21`: First-order rate of transfer from first peripheral to central
  compartment, \\k\_{21}\\ (inverse time units, e.g. 1/h); 2- and
  3-compartment

- `k13`: First-order rate of transfer from central to second peripheral
  compartment, \\k\_{13}\\ (inverse time units, e.g. 1/h); 3-compartment

- `k31`: First-order rate of transfer from second peripheral to central
  compartment,\\k\_{31}\\ (inverse time units, e.g. 1/h); 3-compartment

- `thalf_alpha`: \\t\_{1/2,\alpha}\\ (time units, e.g. h); 1-, 2-, and
  3-compartment

- `thalf_beta`: \\t\_{1/2,\beta}\\ (time units, e.g. h); 2- and
  3-compartment

- `thalf_gamma`: \\t\_{1/2,\gamma}\\ (time units, e.g. h); 3-compartment

- `alpha`: \\\alpha\\; 1-, 2-, and 3-compartment

- `beta`: \\\beta\\; 2- and 3-compartment

- `gamma`: \\\beta\\; 3-compartment

- `trueA`: true A; 1-, 2-, and 3-compartment

- `trueB`: true B; 2- and 3-compartment

- `trueC`: true C; 3-compartment

- `fracA`: fractional A; 1-, 2-, and 3-compartment

- `fracB`: fractional B; 2- and 3-compartment

- `fracC`: fractional C; 3-compartment

- `AUCinf`: Area under the concentration-time curve to infinity (single
  dose)

- `AUCtau`: Area under the concentration-time curve over the dosing
  interval at steady state

- `Cmax`: Maximum concentration after a single dose

- `Cmaxss`: Maximum concentration over the dosing interval at steady
  state

- `Tmax`: Time after dose of maximum concentration

- `AUCinf_dose_normalized`: Dose-normalized area under the
  concentration-time curve to infinity (single dose)

- `AUCtau_dose_normalized`: Dose-normalized area under the
  concentration-time curve over the dosing interval at steady state

- `Cmax_dose_normalized`: Dose-normalized maximum concentration after a
  single dose

- `Cmaxss_dose_normalized`: Dose-normalized maximum concentration over
  the dosing interval at steady state

- `step`: Time increment used when estimating NCA parameters.

The input parameters with standardized names (`dose`, `V1`, `V2`, `V3`,
`CL`, `Q2`, and `Q3`) are also returned in the list, and if provided,
additional PK parameters of `ka`, `tlag`, `tinf` and `dur` are also
returned in the list. All inputs may be scalars or vectors.

## References

Shafer S. L. `CONVERT.XLS`

Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics:
Concepts and Applications (4th). Lippincott Williams & Wilkins,
Philadelphia, 2010.

## Author

Justin Wilkins, <justin.wilkins@occams.com>

Bill Denney, <wdenney@humanpredictions.com>

## Examples

``` r
params <- calc_derived(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
params <- calc_derived_1cpt(CL=16, V=25)
params <- calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
params <- calc_derived_3cpt(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
```
