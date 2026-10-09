# Design: pediatric weight/height simulation from NHANES

Date: 2026-10-08
Source: `C:\Users\justin\Downloads\nhanes.R` (stand-alone script)
Branch: `nhanes-peds`

## Goal

Turn the stand-alone NHANES script into pmxTools functionality that simulates
virtual pediatric populations with realistic, jointly distributed body weight
and height, for use as covariates in pharmacometric simulations.

The script simulated weight only, for ages 6-17, from data downloaded at run
time. The package version:

- ships a pre-processed NHANES reference dataset (works offline, testable on
  CRAN);
- covers ages 2-17;
- simulates weight **and** height jointly, preserving their correlation;
- separates simulation from validation (comparison table and plot).

## 1. Bundled reference data: `nhanes_peds`

Built by `data-raw/nhanes_peds.R` (excluded from the package build), which:

1. Downloads `DEMO_<suffix>.XPT` and `BMX_<suffix>.XPT` from
   `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<year>/DataFiles/` for:

   | CYCLE   | YEAR | SUFFIX |
   |---------|------|--------|
   | 2013-14 | 2013 | H      |
   | 2015-16 | 2015 | I      |
   | 2017-18 | 2017 | J      |
   | 2021-23 | 2021 | L      |

2. Joins on `SEQN`; keeps children aged 2-17 whole years at examination
   (`floor(RIDEXAGM / 12)`), sex 1 or 2, finite positive `WTMEC2YR`, and at
   least one of `BMXWT` / `BMXHT` measured. A missing weight or height is kept
   as `NA`; records are **not** dropped because one measure is missing. Only
   children with neither measure are excluded.
3. Saves `data/nhanes_peds.rda` (xz compression).

`haven` is used only by the data-raw script; it is not a package dependency.
If any cycle lacks `RIDEXAGM` or `WTMEC2YR`, stop and ask; do not substitute
another variable silently.

Columns (one row per child):

| Column       | Type      | Definition                                   |
|--------------|-----------|----------------------------------------------|
| `CYCLE`      | character | NHANES release, e.g. `"2017-18"`             |
| `SEQN`       | numeric   | NHANES respondent sequence number            |
| `AGE_MONTHS` | numeric   | Age in months at examination (`RIDEXAGM`)    |
| `AGE`        | integer   | Whole years, `floor(AGE_MONTHS / 12)`        |
| `SEX`        | character | `"Male"` or `"Female"` (`RIAGENDR`)          |
| `WT`         | numeric   | Body weight, kg (`BMXWT`); `NA` if not measured |
| `HT`         | numeric   | Standing height, cm (`BMXHT`); `NA` if not measured |
| `MEC_WT`     | numeric   | Two-year MEC exam weight (`WTMEC2YR`)        |

Probe of the source files (2026-10-08): all four cycles have `RIDEXAGM`,
`WTMEC2YR`, `BMXWT` and `BMXHT`; 11,323 children aged 2-17, of whom 11,071
have both weight and height.

Documented in `R/data.R`: source URLs, variable definitions, exclusions, the
number of children with weight or height missing, public-domain status, and a
statement that the equal cycle mixture used by the simulator is a modelling
choice, not an official pooled NHANES weight.

`DESCRIPTION` gains `LazyData: true`.

## 2. Simulator: `sample_nhanes()`

```r
sample_nhanes(
  n = 500,
  ages = 2:17,
  sex = c("Male", "Female"),
  vars = c("WT", "HT"),
  cycles = NULL,
  bandwidth_factor = 1,
  seed = NULL,
  data = nhanes_peds
)
```

- `n`: simulated children per age/sex stratum.
- `vars`: which measures to simulate: `"WT"`, `"HT"`, or both (default).
  Only children with **all** requested measures present are eligible as
  donors, so a child missing height is still used when only weight is
  requested. Both requested: joint 2-D kernel. One requested: 1-D kernel
  (identical to the script's method for `"WT"`).
- `ages`: whole years; must all be present in `data`.
- `sex`: subset of `"Male"`, `"Female"`.
- `cycles`: `NULL` uses every cycle in `data`; otherwise a subset of
  `unique(data$CYCLE)`.
- `bandwidth_factor`: multiplier on the kernel bandwidth. `0` gives plain
  weighted resampling; `>1` smooths more.
- `seed`: optional; if given, the RNG is seeded for the call and the caller's
  previous RNG state (`.Random.seed`, or its absence) is restored on exit.
- `data`: reference data with the `nhanes_peds` columns.

### Method (per AGE x SEX stratum)

Let `d = length(vars)` (1 or 2). Donors are the stratum's children with every
requested measure present.

1. Sampling probability of each donor:
   `p_i = MEC_WT_i / sum(MEC_WT in its cycle) * 1 / n_cycles`, then normalised
   to sum to 1 within the stratum. This is the script's equal cycle mixture:
   no cycle dominates because of a larger total MEC weight.
2. Work on `y = log(vars)` (one or two columns).
3. For each dimension `j`: weighted SD `s_j`, weighted IQR `q_j`, robust scale
   `sigma_j = min(s_j, q_j / 1.349)`, falling back to `s_j` and then to 0 when
   not finite or not positive (the script's fallbacks).
4. If `d = 2`: weighted correlation `r` of the two log variables (0 if either
   scale is 0); `R = [[1, r], [r, 1]]`. If `d = 1`: `R = 1`.
5. Kish effective sample size `n_eff = sum(p)^2 / sum(p^2)`.
6. Bandwidth matrix
   `H = (bandwidth_factor * 0.9)^2 * n_eff^(-2 / (d + 4)) * S %*% R %*% S`,
   with `S = diag(sigma)`. For `d = 1` this is exactly the script's
   `h = 0.9 * scale * n_eff^(-1/5)` (squared). For `d = 2` the exponent is
   Silverman's rule of thumb (constant `(4 / (d + 2))^(1 / (d + 4)) = 1`),
   keeping the script's 0.9 multiplier.
7. Draw `n` donor indices with `sample.int(prob = p, replace = TRUE)`, add
   `MASS::mvrnorm(n, rep(0, d), H)` noise to the donors' `y`, and
   exponentiate.

### Return value

A tibble ordered by AGE, SEX with columns `ID`, `AGE`, `SEXN` (1 = Male,
2 = Female), `SEX`, the requested `vars` (`WT` and/or `HT`), `SOURCE_CYCLE`,
`SOURCE_SEQN`.

Attributes:

- `"kernels"`: tibble per stratum with `AGE`, `SEX`, `N_NHANES` (eligible
  donors), `N_EFF`, `H_<var>` for each requested var (bandwidth SDs on the log
  scale, `sqrt(diag(H))`), and `RHO` (`r`) when both vars are requested.
- `"cycles"`: the cycles used, and `"vars"`: the vars simulated, so the
  comparison uses the same reference.

## 3. Validation helpers

### `compare_nhanes(sim, data = nhanes_peds)`

Builds the reference from the cycle-balanced MEC weights of `data`, restricted
to the cycles in `attr(sim, "cycles")`, the ages/sexes present in `sim`, and
the same donor eligibility as the simulation (all of `attr(sim, "vars")`
present). Returns a long tibble, one row per `AGE` x `SEX` x `VARIABLE`
(`VARIABLE` is each simulated var, plus `BMI = WT / (HT / 100)^2` when both
were simulated), with:

- `N`, `Mean`, `SD`, `P05`, `Q1`, `Median`, `Q3`, `P95`, each suffixed `_REF`
  (weighted) and `_SIM` (unweighted);
- `MEAN_RATIO`, `MEDIAN_RATIO`, `P05_RATIO`, `P95_RATIO` (SIM / REF);
- when both vars were simulated, `RHO_REF`, `RHO_SIM`: correlation of log WT
  and log HT per AGE x SEX (repeated on each variable's row).

Errors if `sim` lacks the `"cycles"` attribute or required columns.

### `plot_nhanes(comparison)`

ggplot of P05, Median and P95 against age; colour = source (NHANES vs
Simulated), linetype = percentile; `facet_grid(VARIABLE ~ SEX, scales =
"free_y")`; `theme_bw()`. Returns the ggplot object.

No CSV export in the package.

## 4. Errors and testing

### Input validation (`stop()`, as elsewhere in the package)

- `ages`, `sex` or `cycles` values absent from `data` (message names them);
- `vars` not a non-empty subset of `c("WT", "HT")`;
- `n` not a single positive whole number;
- `bandwidth_factor` not a single non-negative number;
- `data` missing required columns;
- any CYCLE x AGE x SEX stratum with zero eligible donors for the requested
  `vars` (message lists them).

Degenerate strata (zero spread in one dimension) get a zero bandwidth in that
dimension; `MASS::mvrnorm()` handles the singular `H`.

### Tests (testthat; offline, bundled data; no vdiffr snapshots)

- Shape: `n * length(ages) * length(sex)` rows; expected columns; `WT`, `HT`
  finite and positive.
- `bandwidth_factor = 0`: each simulated (WT, HT) equals its donor's values.
- `vars = "WT"`: only `WT` returned; donors include children missing height;
  matches the script's 1-D bandwidth `0.9 * scale * n_eff^(-1/5)`.
- `vars = c("WT", "HT")`: no donor has a missing WT or HT.
- `seed`: identical output for the same seed; `.Random.seed` unchanged by the
  call.
- Fidelity (fixed seed, large `n`): per-stratum medians of WT and HT within
  tolerance of the weighted reference medians; `RHO_SIM` within tolerance of
  `RHO_REF`.
- `cycles` subset: `SOURCE_CYCLE` limited to the subset; comparison uses it.
- Each validation error triggers.
- `compare_nhanes()`: expected columns; 3 rows per stratum.
- `plot_nhanes()`: returns a `ggplot`.
- `nhanes_peds`: expected columns; ages 2-17; no missing `MEC_WT`; every row
  has at least one of `WT` / `HT`; some rows have exactly one missing.

### Other

- NEWS bullet under `# pmxTools (development version)`.
- Runnable `@examples`.
- `.Rbuildignore`: `^data-raw$`, `^dev$`.

## Out of scope

- Downloading other cycles at run time (could be added later as a helper with
  `haven` in Suggests).
- Ages outside 2-17, other body measures, continuous-age simulation.
- Changing the existing `sample_omega()` / `sample_uncert()` seed handling
  (flagged separately: they call `set.seed()` and overwrite the user's RNG
  state).

## Revision 2026-10-09

- Functions renamed to `sample_nhanes()`, `compare_nhanes()`, `plot_nhanes()`
  (dataset keeps the name `nhanes_peds`).
- New argument `method = c("resample", "smooth")`, default `"resample"`:
  each simulated child takes one donor child's recorded values unchanged
  (weight and height stay paired; no noise is drawn). `"smooth"` is the
  kernel method described above. `bandwidth_factor` applies only to
  `"smooth"` and warns if supplied with `"resample"`. Kernel diagnostics
  report `H_WT`/`H_HT` = 0 when resampling; the result carries a `"method"`
  attribute.
- `AGE` is always returned as integer.
