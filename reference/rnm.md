# Read NONMEM 7.2+ output into an R object.

Read NONMEM 7.2+ output into an R object.

## Usage

``` r
rnm(
  index,
  prefix = "run",
  pathNM,
  ndig = 3,
  ndigB = 3,
  ndigP = 1,
  Pci = 95,
  ext = ".lst",
  extmod = ".mod",
  Pvalues = TRUE,
  RawCI = FALSE,
  ...
)
```

## Arguments

- index:

  The NONMEM model index, i.e. the numeric part of the filename assuming
  it follows the convention 'run123.mod'.

- prefix:

  The NONMEM model prefix, assuming it follows the convention
  'run123.mod'. The default is `"run"`.

- pathNM:

  The path to the NONMEM output. This should not contain a trailing
  slash.

- ndig:

  Number of significant digits to use. The default is 3.

- ndigB:

  Number of significant digits to use. The default is 3.

- ndigP:

  Number of digits after the decimal point to use for percentages. The
  default is 1.

- Pci:

  Asymptotic confidence interval to apply when reporting parameter
  uncertainty. The default is 95.

- ext:

  NONMEM output file extension. The default is `".lst"`.

- extmod:

  NONMEM control stream file extension. The default is `"mod"`.

- Pvalues:

  Report P-values for parameters? The default is `TRUE`.

- RawCI:

  Report confidence intervals without estimate? The default is `FALSE`.

- ...:

  Additional arguments.

## Value

A list containing information extracted from the NONMEM output.

## Details

The output list is composed of the following objects:

- *Theta*: A data frame describing the structural (fixed-effect)
  parameters, containing parameter name, estimated value, standard error
  (SE), coefficient of variation (CV), lower and upper confidence limits
  (CIL and CIU, based on `Pci`), and P-value, calculated as
  `2\*(1-pnorm(abs(theta/theta.se)))`.

- *Eta*: A data frame describing the interindividual random-effects
  parameters, containing estimated value, standard error (SE),
  coefficient of variation (CV, calculated as `abs(100*(SE/OMEGA))`),
  coefficient of variation (EtaCV, calculated as `100*sqrt(OMEGA)`), and
  shrinkage.

- *Epsilon*: A data frame describing the residual random-effects
  parameters, containing estimated value, standard error (SE),
  coefficient of variation (CV, calculated as `abs(100*(SE/OMEGA))`),
  coefficient of variation (EtaCV, calculated as `100*sqrt(SIGMA)`), and
  shrinkage.

- *CorTheta*: A data frame containing the correlation matrix for fixed
  effects (`THETA`).

- *CorOmega*: A data frame containing the correlation matrix for
  interindividual random effects (`OMEGA`).

- *CorSigma*: A data frame containing the correlation matrix for
  residual random effects (`SIGMA`).

- *OmegaMatrix*: A data frame containing the `OMEGA` matrix.

- *SigmaMatrix*: A data frame containing the `SIGMA` matrix.

- *CovMatrixTheta*: A data frame containing the variance-covariance
  matrix for structural parameters (`THETA`).

- *CovMatrix*: A data frame containing the complete variance-covariance
  matrix.

- *OFV*: The objective function value.

- *ThetaString*: A data frame containing all relevant fixed-effects
  parameter information, suitable for use in a table of parameter
  estimates. Contains parameter name, estimate, standard error,
  coefficient of variation, combined estimate and asymptotic confidence
  interval, and P-value.

- *EtaString*: A data frame containing all relevant interindivudal
  random-effects parameter information, suitable for use in a table of
  parameter estimates. Contains parameter name, estimate (variance),
  standard error, coefficient of variation, percentage value (calculated
  as `100*sqrt(OMEGA)`), and shrinkage.

- *EpsString*: A data frame containing all relevant residual
  random-effects parameter information, suitable for use in a table of
  parameter estimates. Contains parameter name, estimate (variance),
  standard error, coefficient of variation, percentage value (calculated
  as `100*sqrt(SIGMA)`), and shrinkage.

- *RunTime*: Run time.

- *ConditionN*: Condition number.

## See also

NONMEM (<https://www.iconplc.com/solutions/technologies/nonmem>)

## Author

Rik Schoemaker, <rik.schoemaker@occams.com>

Justin Wilkins, <justin.wilkins@occams.com>

## Examples

``` r
if (FALSE) { # \dontrun{
nmOutput <- rnm("run315.lst")
} # }
```
