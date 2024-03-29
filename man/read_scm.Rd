% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_scm.R
\name{read_scm}
\alias{read_scm}
\title{Read PsN SCM output into a format suitable for further use.}
\usage{
read_scm(dir, startPhase = "forward")
}
\arguments{
\item{dir}{A PsN SCM folder (containing \code{scmlog.txt} and \code{short_scmlog.txt}).}

\item{startPhase}{Where to start collating the output; can be \code{"forward"} (the default) or \code{"backward"}.}
}
\value{
A list of data frames, containing
  \item{forward}{all models evaluated during the forward inclusion step of
  covariate model building}
  \item{forwardSummary}{the covariate relationships selected at each forward
  step}
  \item{forwardP}{the P-value used for inclusion during the forward inclusion step}
  \item{backward}{all models evaluated during the backward elimination step of covariate
  model building}
  \item{backwardSummary}{the covariate relationships eliminated at each backward step}
  \item{backwardP}{the P-value used for exclusion during the backward elimination step}
}
\description{
\code{read_scm} returns a summary of a Perl-speaks-NONMEM (PsN, \url{https://uupharmacometrics.github.io/PsN/}) SCM (stepwise covariate modeling)
procedure. It depends on the presence of \code{scmlog.txt} and \code{short_scmlog.txt} files in the
specified directory.
}
\examples{
\dontrun{
scm <- read_scm("E:/DrugX/ModelDevelopment/scm310")
}

}
\seealso{
NONMEM (\url{https://www.iconplc.com/innovation/nonmem/})

Lindbom L, Ribbing J & Jonsson EN (2004). Perl-speaks-NONMEM (PsN) - A Perl module for NONMEM related programming. Computer Methods and Programs in Biomedicine, 75(2), 85-94. \doi{10.1016/j.cmpb.2003.11.003}

Lindbom L, Pihlgren P & Jonsson N (2005). PsN-Toolkit - A collection of computer intensive statistical methods for non-linear mixed effect modeling using NONMEM. Computer Methods and Programs in Biomedicine, 79(3), 241-257. \doi{10.1016/j.cmpb.2005.04.005}

Other NONMEM reading: 
\code{\link{plot_scm}()},
\code{\link{read_nm_all}()},
\code{\link{read_nm_multi_table}()},
\code{\link{read_nmcov}()},
\code{\link{read_nmext}()},
\code{\link{read_nmtables}()},
\code{\link{read_nm}()}
}
\author{
Justin Wilkins, \email{justin.wilkins@occams.com}
}
\concept{NONMEM reading}
