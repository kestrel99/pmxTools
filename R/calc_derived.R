#' Calculate derived pharmacokinetic parameters for a 1-, 2-, or 3-compartment linear model.
#'
#' @param CL Clearance (L/h)
#' @param V1,V Central volume of distribution (L).  Values are synonyms; use only one.
#' @param V2 First peripheral volume of distribution (L)
#' @param V3 Second peripheral volume of distribution (L)
#' @param Q2,Q Intercompartmental clearance from central to first peripheral compartment (L/h).  Values are synonyms; use only one.
#' @param Q3 Intercompartmental clearance from central to second peripheral compartment (L/h)
#' @param type Type of half-life to return. Default is \code{"all"}; see details for other options that are specific to the number of compartments.
#' @param sigdig Number of significant digits to be returned. Default is \code{5}.
#' @param verbose For `calc_derived()`, provide a message indicating the type of model detected.
#' @param ... Passed to the other `calc_derived_*()` functions.
#' 
#' @return Return a list of derived PK parameters for a 3-compartment linear model given provided clearances and volumes based on the `type`.
#' \itemize{ 
#'   \item \code{Vss}: Volume of distribution at steady state, \eqn{V_{ss}} (L); 1-, 2-, and 3-compartment
#'   \item \code{k10}: First-order elimination rate, \eqn{k_{10}} (/h); 1-, 2-, and 3-compartment
#'   \item \code{k12}: First-order rate of transfer from central to first peripheral compartment, \eqn{k_{12}} (/h); 2- and 3-compartment
#'   \item \code{k21}: First-order rate of transfer from first peripheral to central compartment, \eqn{k_{21}} (/h); 2- and 3-compartment
#'   \item \code{k13}: First-order rate of transfer from central to second peripheral compartment, \eqn{k_{13}} (/h); 3-compartment
#'   \item \code{k31}: First-order rate of transfer from second peripheral to central compartment,\eqn{k_{31}} (/h); 3-compartment
#'   \item \code{thalf_alpha}: \eqn{t_{1/2,\alpha}} (h); 2- and 3-compartment
#'   \item \code{thalf_beta}: \eqn{t_{1/2,\beta}} (h); 2- and 3-compartment
#'   \item \code{thalf_gamma}: \eqn{t_{1/2,\gamma}} (h); 3-compartment
#'   \item \code{alpha}: \eqn{\alpha}; 1-, 2-, and 3-compartment
#'   \item \code{beta}: \eqn{\beta}; 2- and 3-compartment
#'   \item \code{gamma}: \eqn{\beta}; 3-compartment
#'   \item \code{trueA}: true A; 1-, 2-, and 3-compartment
#'   \item \code{trueB}: true B; 2- and 3-compartment
#'   \item \code{trueC}: true C; 3-compartment
#'   \item \code{fracA}: fractional A; 1-, 2-, and 3-compartment
#'   \item \code{fracB}: fractional B; 2- and 3-compartment
#'   \item \code{fracC}: fractional C; 3-compartment
#'  }
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @references Shafer S. L. \code{CONVERT.XLS}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams & Wilkins, Philadelphia, 2010.  
#'
#' @examples
#' params <- calc_derived(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
#' @export
calc_derived <- function(..., verbose=FALSE) {
  arg_names <- names(list(...))
  if (all(arg_names %in% names(formals(calc_derived_1cpt)))) {
    if (verbose) message("Detected 1-compartment model")
    calc_derived_1cpt(...)
  } else if (all(arg_names %in% names(formals(calc_derived_2cpt)))) {
    if (verbose) message("Detected 2-compartment model")
    calc_derived_2cpt(...)
  } else if (all(arg_names %in% names(formals(calc_derived_3cpt)))) {
    if (verbose) message("Detected 3-compartment model")
    calc_derived_3cpt(...)
  } else {
    unknown_params <- setdiff(arg_names, names(formals(calc_derived_3cpt)))
    stop(
      "Could not determine model type based on argument names.  Please check the following argument names: ",
      paste(unknown_params, collapse=", ")
    )
  }
}
