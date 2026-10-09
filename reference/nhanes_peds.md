# NHANES pediatric body weight and height reference data

Body weight and standing height of children aged 2 to 17 years at the
examination, from four releases of the US National Health and Nutrition
Examination Survey (NHANES): 2013-2014, 2015-2016, 2017-2018 and August
2021-August 2023. Used as the reference population by
[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md).

## Usage

``` r
nhanes_peds
```

## Format

A data frame with one row per child and 8 columns:

- CYCLE:

  NHANES release, e.g. `"2017-18"`.

- SEQN:

  NHANES respondent sequence number (unique within a release).

- AGE_MONTHS:

  Age in months at the examination (`RIDEXAGM`).

- AGE:

  Age in whole years, `floor(AGE_MONTHS / 12)`.

- SEX:

  `"Male"` or `"Female"` (`RIAGENDR`).

- WT:

  Body weight in kg (`BMXWT`); `NA` if not measured.

- HT:

  Standing height in cm (`BMXHT`); `NA` if not measured.

- MEC_WT:

  Two-year MEC exam weight (`WTMEC2YR`).

## Source

National Center for Health Statistics, NHANES public data files
`DEMO_H`/`BMX_H`, `DEMO_I`/`BMX_I`, `DEMO_J`/`BMX_J` and
`DEMO_L`/`BMX_L`, <https://wwwn.cdc.gov/nchs/nhanes/>.

## Details

Children are included when they have a positive two-year mobile
examination center (MEC) exam weight and at least one of body weight or
standing height. A child with only one of the two measures is kept, with
the other set to `NA`;
[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md)
uses such children only when the missing measure is not requested.

The MEC weights are those of each individual release. They are not
combined into an official pooled weight;
[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md)
gives each release an equal share, which is a modelling choice.

NHANES data are produced by the US National Center for Health Statistics
and are in the public domain. The data were built with
`data-raw/nhanes_peds.R` in the package source.

## See also

[`sample_nhanes()`](https://kestrel99.github.io/pmxTools/reference/sample_nhanes.md)
