# A grid from 0 to `to` in increments of `step` that always ends exactly at
# `to`, so the whole dosing interval is searched, and that also
# contains the break points `at` lying within it.
seq_through <- function(to, step, at = NULL) {
  grid <- seq(0, to, by = step)
  sort(unique(c(grid, to, at[at >= 0 & at <= to])))
}

# Times within the dosing interval where a steady-state profile has a kink, and
# Cmaxss can lie off the grid: the end of an infusion, and the start and end of
# lagged or zero-order absorption (the lag shifts them by tlag modulo tau).
ss_break_points <- function(tau, tinf = NULL, tlag = NULL, dur = NULL) {
  c(tinf, c(tlag, sum(tlag, dur)) %% tau)
}

# Single-dose Cmax and Tmax from a concentration function of time, searched on
# a grid from 0 to `tend` that is doubled until the maximum lies inside it
# (rather than at its end, or nowhere when absorption has not yet begun).  The
# grid also includes the regimen's break points `at` (end of infusion, lag, end
# of zero-order absorption), where the peak can fall between grid steps.
find_peak <- function(conc_fun, tend, step, dose, at = NULL, max_doublings = 8) {
  if (dose == 0) return(list(Cmax = 0, Tmax = 0))
  for (i in 0:max_doublings) {
    t <- sort(unique(c(seq(0, tend, by = step), at)))
    C <- quiet_flipflop(conc_fun(t))
    peak <- which.max(C)
    if (peak < length(t) && C[peak] > 0) break
    tend <- tend * 2
  }
  if (peak == length(t) || C[peak] <= 0) {
    warning("No concentration peak was found within ", t[length(t)], " time units; Cmax and Tmax are unreliable.",
            call. = FALSE)
  }
  list(Cmax = C[peak], Tmax = t[peak])
}

# The concentration functions derive their rate constants with
# calc_derived_*(), which warns about flip-flop kinetics on every call.  The
# calling calc_derived_*() warns once itself, so silence the repeats.
muffle_flipflop <- function(w) {
  if (startsWith(conditionMessage(w), "Flip-flop")) invokeRestart("muffleWarning")
}
quiet_flipflop <- function(expr) withCallingHandlers(expr, warning = muffle_flipflop)

# Inputs that may be supplied as vectors, one element per parameter set.
derived_vector_args <- c("CL", "V", "V1", "V2", "V3", "Q", "Q2", "Q3", "ka", "dur", "tlag", "tinf", "dose", "tau")

#' Evaluate calc_derived_*() separately for each recycled parameter set
#'
#' The dose-dependent (NCA) outputs are computed from one simulated profile per
#' parameter set, so when `dose` is given with vector inputs each set is
#' evaluated on its own and the results are combined element-wise.  Messages
#' are emitted once rather than once per set.
#'
#' @param fun The calling `calc_derived_*()` function.
#' @param args Named list of the arguments supplied to `fun`.
#' @returns `NULL` when no row-wise evaluation is needed (scalar inputs, or no
#'   `dose`), otherwise the combined result.
#' @noRd
derived_rowwise <- function(fun, args) {
  vec <- Filter(Negate(is.null), args[intersect(names(args), derived_vector_args)])
  n_each <- lengths(vec)
  n <- max(c(1, n_each))
  bad <- n_each != 1 & n_each != n
  if (any(bad)) {
    stop("Parameter inputs must have length 1 or ", n, "; check: ",
         paste(names(vec)[bad], collapse = ", "), call. = FALSE)
  }
  if (n == 1 || is.null(args[["dose"]])) return(NULL)

  rows <- lapply(seq_len(n), function(i) {
    row_args <- args
    row_args[names(vec)] <- lapply(vec, function(x) x[min(i, length(x))])
    run <- function() do.call(fun, row_args)
    withCallingHandlers(
      if (i == 1) run() else suppressMessages(run()),
      warning = muffle_flipflop
    )
  })
  V1 <- if (is.null(vec[["V1"]])) vec[["V"]] else vec[["V1"]]
  maybe_warn_flipflop(ka = vec[["ka"]], k10 = vec[["CL"]] / V1)

  if (!is.list(rows[[1]])) return(unlist(rows))
  fields <- names(rows[[1]])
  out <- stats::setNames(lapply(fields, function(f) unlist(lapply(rows, `[[`, f))), fields)
  if ("step" %in% fields) out$step <- rows[[1]]$step
  out
}

#' Select the requested outputs of calc_derived_*()
#'
#' @param out Named list of all outputs; `NULL` entries are not available for
#'   the inputs given.
#' @param type `"all"`, one output name, or a vector of output names.
#' @returns The list without `NULL` entries for `"all"`, the value itself for
#'   one name, or a list of the named outputs for several names.
#' @noRd
select_derived <- function(out, type) {
  out <- out[!vapply(out, is.null, logical(1))]
  if (identical(type, "all")) return(out)
  if (!is.character(type) || length(type) == 0) {
    stop("`type` must be \"all\" or a character vector of output names.", call. = FALSE)
  }
  invalid <- setdiff(type, names(out))
  if (length(invalid) > 0) {
    stop("Invalid value for `type`: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (length(type) == 1) out[[type]] else out[type]
}

#' Calculate derived pharmacokinetic parameters for a 1-, 2-, or 3-compartment
#' linear model.
#'
#' @param CL Clearance (volume per time units, e.g. L/h)
#' @param V1,V Central volume of distribution (volume units, e.g. L).  Values
#'   are synonyms; use only one.
#' @param V2 First peripheral volume of distribution (volume units, e.g. L)
#' @param V3 Second peripheral volume of distribution (volume units, e.g. L)
#' @param Q2,Q Intercompartmental clearance from central to first peripheral
#'   compartment (volume per time units, e.g. L/h).  Values are synonyms; use
#'   only one.
#' @param Q3 Intercompartmental clearance from central to second peripheral
#'   compartment (volume per time units, e.g. L/h)
#' @param ka Absorption rate (inverse time units, e.g. 1/h)
#' @param dur Duration of zero-order absorption (time units, e.g. h)
#' @param tlag Absorption lag time (time units, e.g. h)
#' @param tinf Duration of infusion (time units, e.g. h)
#' @param dose Dose (amount units, e.g. mg)
#' @param tau Duration of interdose interval (time units, e.g. h; defaults to 24)
#' @param step Time increment of the grid searched for \code{Cmax}, \code{Tmax} and \code{Cmaxss} (time units, e.g. h; defaults to 0.1). The end of an infusion and the start and end of absorption are always included. \code{AUCtau} is calculated as \code{dose/CL}, which is exact at steady state for these linear models, so it does not depend on \code{step}.
#' @param type Parameters to return. Default is \code{"all"}.  If not
#'   \code{"all"}, this may be one name from the return value list, giving
#'   that value alone, or a vector of names, giving a list of those values.
#' @param sigdig Number of significant digits to be returned. Default is
#'   \code{5}.
#' @param verbose For \code{calc_derived()}, provide a message indicating the
#'   type of model detected.
#' @param ... Passed to the other \code{calc_derived_*()} functions.
#'
#' @return Return a list of derived PK parameters for a 1-, 2-, or 3-compartment
#'   linear model given provided clearances and volumes based on the
#'   \code{type}. If a dose is provided, estimated non-compartmental analysis (NCA) parameters will
#'   be provided as well, based on simulation of single-dose and (if `tau` is specified) steady-state time courses.
#' \itemize{
#'   \item \code{Vss}: Volume of distribution at steady state, \eqn{V_{ss}} (volume units, e.g. L); 1-, 2-, and 3-compartment
#'   \item \code{k10}: First-order elimination rate, \eqn{k_{10}} (inverse time units, e.g. 1/h); 1-, 2-, and 3-compartment
#'   \item \code{k12}: First-order rate of transfer from central to first peripheral compartment, \eqn{k_{12}} (inverse time units, e.g. 1/h); 2- and 3-compartment
#'   \item \code{k21}: First-order rate of transfer from first peripheral to central compartment, \eqn{k_{21}} (inverse time units, e.g. 1/h); 2- and 3-compartment
#'   \item \code{k13}: First-order rate of transfer from central to second peripheral compartment, \eqn{k_{13}} (inverse time units, e.g. 1/h); 3-compartment
#'   \item \code{k31}: First-order rate of transfer from second peripheral to central compartment,\eqn{k_{31}} (inverse time units, e.g. 1/h); 3-compartment
#'   \item \code{thalf_alpha}: \eqn{t_{1/2,\alpha}} (time units, e.g. h); 1-, 2-, and 3-compartment
#'   \item \code{thalf_beta}: \eqn{t_{1/2,\beta}} (time units, e.g. h); 2- and 3-compartment
#'   \item \code{thalf_gamma}: \eqn{t_{1/2,\gamma}} (time units, e.g. h); 3-compartment
#'   \item \code{alpha}: \eqn{\alpha}; 1-, 2-, and 3-compartment
#'   \item \code{beta}: \eqn{\beta}; 2- and 3-compartment
#'   \item \code{gamma}: \eqn{\beta}; 3-compartment
#'   \item \code{trueA}: true A; 1-, 2-, and 3-compartment
#'   \item \code{trueB}: true B; 2- and 3-compartment
#'   \item \code{trueC}: true C; 3-compartment
#'   \item \code{fracA}: fractional A; 1-, 2-, and 3-compartment
#'   \item \code{fracB}: fractional B; 2- and 3-compartment
#'   \item \code{fracC}: fractional C; 3-compartment
#'   \item \code{AUCinf}: Area under the concentration-time curve to infinity (single dose)
#'   \item \code{AUCtau}: Area under the concentration-time curve over the dosing interval at steady state
#'   \item \code{Cmax}: Maximum concentration after a single dose
#'   \item \code{Cmaxss}: Maximum concentration over the dosing interval at steady state
#'   \item \code{Tmax}: Time after dose of maximum concentration
#'   \item \code{AUCinf_dose_normalized}: Dose-normalized area under the concentration-time curve to infinity (single dose)
#'   \item \code{AUCtau_dose_normalized}: Dose-normalized area under the concentration-time curve over the dosing interval at steady state
#'   \item \code{Cmax_dose_normalized}: Dose-normalized maximum concentration after a single dose
#'   \item \code{Cmaxss_dose_normalized}: Dose-normalized maximum concentration over the dosing interval at steady state
#'   \item \code{step}: Time increment of the grid searched for \code{Cmax}, \code{Tmax} and \code{Cmaxss}.
#'  }
#'
#' The input parameters with standardized names (\code{dose}, \code{V1}, \code{V2},
#' \code{V3}, \code{CL}, \code{Q2}, and \code{Q3}) are also returned in the
#' list, and if provided, additional PK parameters of `ka`, `tlag`, `tinf` and `dur` are also
#' returned in the list.  All PK parameter inputs (including `dose` and `tau`)
#' may be scalars or vectors; vectors must all have the same length, and
#' scalars are recycled to it.  Each element of every returned parameter then
#' corresponds to one parameter set, with NCA parameters calculated separately
#' for each set.  `step`, `type` and `sigdig` must be single values.
#'
#' Clearances, volumes, `ka`, `dur`, `tinf`, `tau` and `step` must be positive;
#' `tlag` and `dose` must be non-negative.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @author Bill Denney, \email{wdenney@@humanpredictions.com}
#' @references Shafer S. L. \code{CONVERT.XLS}
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and
#'   Pharmacodynamics: Concepts and Applications (4th). Lippincott Williams &
#'   Wilkins, Philadelphia, 2010.
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

update_tend <- function(tend, guess_tmax) {
  mask_tend_tmax <- tend < guess_tmax
  if (any(mask_tend_tmax)) {
    tend[mask_tend_tmax] <-
      tend * 2^ceiling(log2(guess_tmax[mask_tend_tmax]/tend[mask_tend_tmax]))
  }
  tend
}

maybe_warn_flipflop <- function(ka, k10) {
  if(!is.null(ka)) {
    if(any(ka < k10)) {
      warning("Flip-flop kinetics detected. ka is lower than k10, and half-life may therefore not be similar to the half-life calculated by typical NCA.\n")
    }
  }
}

#' @rdname calc_derived
#' @examples
#' params <- calc_derived_1cpt(CL=16, V=25)
#' @export
calc_derived_1cpt <- function(CL, V=NULL, V1=NULL, ka=NULL, dur=NULL, tlag=NULL, tinf=NULL, dose=NULL, tau=NULL, step=0.1, type="all", sigdig=5) {
  rowwise <- derived_rowwise(calc_derived_1cpt, mget(names(match.call())[-1]))
  if (!is.null(rowwise)) return(rowwise)

  if (!xor(is.null(V), is.null(V1))) {
    stop("Exactly one of V or V1 may be provided since they are considered synonyms.")
  } else if (!is.null(V)) {
    # Use the preferred synonym
    V1 <- V
  }
  check_derived_inputs(
    positive = list(CL = CL, V1 = V1, ka = ka, dur = dur, tinf = tinf, tau = tau, step = step),
    nonnegative = list(tlag = tlag, dose = dose)
  )
  k10 <- CL/V1
  Vss <- V1
  trueA <- 1/V1
  thalf <- log(2)/k10
  alpha <- k10

  AUCinf <- NULL
  AUCtau  <- NULL
  Cmax   <- NULL
  Cmaxss <- NULL
  Tmax   <- NULL
  AUCinf_dose_normalized <- NULL
  AUCtau_dose_normalized <- NULL
  Cmax_dose_normalized   <- NULL
  Cmaxss_dose_normalized <- NULL

  guess_tmax <- 0
  if(!is.null(ka)) guess_tmax <- 1/ka
  if(!is.null(dur)) guess_tmax <- dur
  if(!is.null(tinf)) guess_tmax <- tinf
  if(!is.null(tlag)) guess_tmax <- guess_tmax + tlag

  tend <- update_tend(tend = 24, guess_tmax = guess_tmax)

  if(!is.null(dose)) {
    message("dose present. NCA parameters will be estimated.")
    if(!is.null(tau)) message("tau present. Steady state NCA parameters will be estimated.")
    AUCinf <- dose/CL
    AUCinf_dose_normalized <- AUCinf/dose

    # IV bolus
    if(is.null(dur) & is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("1-compartment bolus detected.")
      peak <- find_peak(function(t) calc_sd_1cmt_linear_bolus(t=t, dose=dose, CL=CL, V1=V1), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_1cmt_linear_bolus(tad = t, dose=dose, tau=tau, CL=CL, V1=V1)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }

    # IV infusion
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tinf)) {
      message("1-compartment infusion detected.")
      peak <- find_peak(function(t) calc_sd_1cmt_linear_infusion(t=t, dose=dose, CL=CL, V1=V1, tinf = tinf), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & !is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_1cmt_linear_infusion(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, tinf=tinf)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral first-order
    if(is.null(dur) & !is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("1-compartment first-order oral detected.")
      peak <- find_peak(function(t) calc_sd_1cmt_linear_oral_1(t=t, dose=dose, ka=ka, CL=CL, V1=V1), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & !is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_1cmt_linear_oral_1(tad = t, dose=dose, tau=tau, ka=ka, CL=CL, V1=V1)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral first-order lag
    if(is.null(dur) & !is.null(ka) & !is.null(tlag) & is.null(tinf)) {
      message("1-compartment first-order oral with lag time detected.")
      peak <- find_peak(function(t) calc_sd_1cmt_linear_oral_1_lag(t=t, dose=dose, ka=ka, CL=CL, V1=V1, tlag=tlag), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & !is.null(ka) & !is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_1cmt_linear_oral_1_lag(tad = t, dose=dose, tau=tau, ka=ka, CL=CL, V1=V1, tlag=tlag)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral zero-order
    if(!is.null(dur) & is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("1-compartment zero-order oral detected.")
      peak <- find_peak(function(t) calc_sd_1cmt_linear_oral_0(t=t, dose=dose, CL=CL, V1=V1, dur=dur), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(!is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_1cmt_linear_oral_0(tad = t, dose=dose, tau=tau, dur=dur, CL=CL, V1=V1)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral zero-order lag
    if(!is.null(dur) & is.null(ka) & !is.null(tlag) & is.null(tinf)) {
      message("1-compartment zero-order oral with lag time detected.")
      peak <- find_peak(function(t) calc_sd_1cmt_linear_oral_0_lag(t=t, dose=dose, CL=CL, V1=V1, tlag=tlag, dur=dur), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(!is.null(dur) & is.null(ka) & !is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_1cmt_linear_oral_0_lag(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, tlag=tlag, dur=dur)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
  }

  out <-
    list(
      k10=signif(CL/V1, sigdig),
      Vss=signif(Vss, sigdig),
      thalf=signif(thalf, sigdig),
      alpha=signif(alpha, sigdig),
      trueA=signif(1/V1, sigdig),
      fracA=1,
      tau=tau,
      AUCinf=AUCinf,
      AUCtau=AUCtau,
      Cmax=Cmax,
      Cmaxss=Cmaxss,
      Tmax=Tmax,
      AUCinf_dose_normalized=AUCinf_dose_normalized,
      AUCtau_dose_normalized=AUCtau_dose_normalized,
      Cmax_dose_normalized=Cmax_dose_normalized,
      Cmaxss_dose_normalized=Cmaxss_dose_normalized,
      step=step,
      # Include the macro parameters
      CL=CL,
      V1=V1,
      ka=ka,
      tlag=tlag,
      dur=dur,
      tinf=tinf,
      dose=dose
    )

  o <- select_derived(out, type)

  maybe_warn_flipflop(ka = ka, k10 = k10)
  o
}

#' @rdname calc_derived
#' @examples
#' params <- calc_derived_2cpt(CL=16, V1=25, V2=50, Q=0.5)
#' @export
calc_derived_2cpt <- function(CL, V1=NULL, V2, Q2=NULL, V=NULL, Q=NULL, dur=NULL, tinf=NULL, ka=NULL, tlag=NULL, dose=NULL, tau=NULL, step=0.1, type="all", sigdig=5) {
  rowwise <- derived_rowwise(calc_derived_2cpt, mget(names(match.call())[-1]))
  if (!is.null(rowwise)) return(rowwise)

  if (!xor(is.null(V), is.null(V1))) {
    stop("Exactly one of V or V1 may be provided since they are considered synonyms.")
  } else if (!is.null(V)) {
    # Use the preferred synonym
    V1 <- V
  }
  if (!xor(is.null(Q), is.null(Q2))) {
    stop("Exactly one of Q or Q2 may be provided since they are considered synonyms.")
  } else if (!is.null(Q)) {
    # Use the preferred synonym
    Q2 <- Q
  }
  check_derived_inputs(
    positive = list(CL = CL, V1 = V1, V2 = V2, Q2 = Q2, ka = ka, dur = dur, tinf = tinf, tau = tau, step = step),
    nonnegative = list(tlag = tlag, dose = dose)
  )

  ### variables
  k10   <- CL/V1
  k12   <- Q2/V1
  k21   <- Q2/V2

  a0 <- k10 * k21
  a1 <- -(k10 + k12 + k21)

  i1 <- (-a1 + sqrt(a1*a1-4*a0))/2
  i2 <- (-a1 - sqrt(a1*a1-4*a0))/2

  c1 <- (k21 - i1)/(i2 - i1)/V1
  c2 <- (k21 - i2)/(i1 - i2)/V1

  AUCinf <- NULL
  AUCtau  <- NULL
  Cmax   <- NULL
  Cmaxss <- NULL
  Tmax   <- NULL
  AUCinf_dose_normalized <- NULL
  AUCtau_dose_normalized <- NULL
  Cmax_dose_normalized   <- NULL
  Cmaxss_dose_normalized <- NULL

  guess_tmax <- 0
  if(!is.null(ka)) guess_tmax <- 1/ka
  if(!is.null(dur)) guess_tmax <- dur
  if(!is.null(tinf)) guess_tmax <- tinf
  if(!is.null(tlag)) guess_tmax <- guess_tmax + tlag

  tend <- update_tend(tend = 24, guess_tmax = guess_tmax)

  if(!is.null(dose)) {
    message("dose present. NCA parameters will be estimated.")
    if(!is.null(tau)) message("tau present. Steady state NCA parameters will be estimated.")
    AUCinf <- dose/CL
    AUCinf_dose_normalized <- AUCinf/dose

    # IV bolus
    if(is.null(dur) & is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("2-compartment bolus detected.")
      peak <- find_peak(function(t) calc_sd_2cmt_linear_bolus(t=t, dose=dose, CL=CL, V1=V1, V2=V2, Q=Q), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_2cmt_linear_bolus(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, V2=V2, Q=Q)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }

    # IV infusion
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tinf)) {
      message("2-compartment infusion detected.")
      peak <- find_peak(function(t) calc_sd_2cmt_linear_infusion(t=t, dose=dose, CL=CL, V1=V1, V2=V2, Q=Q, tinf = tinf), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & !is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_2cmt_linear_infusion(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, V2=V2, Q=Q, tinf=tinf)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral first-order
    if(is.null(dur) & !is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("2-compartment first-order oral detected.")
      peak <- find_peak(function(t) calc_sd_2cmt_linear_oral_1(t=t, dose=dose, ka=ka, CL=CL, V1=V1, V2=V2, Q=Q), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & !is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_2cmt_linear_oral_1(tad = t, dose=dose, tau=tau, ka=ka, CL=CL, V1=V1, V2=V2, Q=Q)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral first-order lag
    if(is.null(dur) & !is.null(ka) & !is.null(tlag) & is.null(tinf)) {
      message("2-compartment first-order oral with lag time detected.")
      peak <- find_peak(function(t) calc_sd_2cmt_linear_oral_1_lag(t=t, dose=dose, ka=ka, CL=CL, V1=V1, V2=V2, Q=Q, tlag=tlag), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & !is.null(ka) & !is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_2cmt_linear_oral_1_lag(tad = t, dose=dose, tau=tau, ka=ka, CL=CL, V1=V1, V2=V2, Q=Q, tlag=tlag)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral zero-order
    if(!is.null(dur) & is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("2-compartment zero-order oral detected.")
      peak <- find_peak(function(t) calc_sd_2cmt_linear_oral_0(t=t, dose=dose, CL=CL, V1=V1, V2=V2, Q=Q, dur=dur), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(!is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_2cmt_linear_oral_0(tad = t, dose=dose, tau=tau, dur=dur, CL=CL, V1=V1, V2=V2, Q=Q)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral zero-order lag
    if(!is.null(dur) & is.null(ka) & !is.null(tlag) & is.null(tinf)) {
      message("2-compartment zero-order oral with lag time detected.")
      peak <- find_peak(function(t) calc_sd_2cmt_linear_oral_0_lag(t=t, dose=dose, CL=CL, V1=V1, V2=V2, Q=Q, tlag=tlag, dur=dur), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(!is.null(dur) & is.null(ka) & !is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_2cmt_linear_oral_0_lag(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, V2=V2, Q=Q, tlag=tlag, dur=dur)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
  }

  out <-
    list(
      k10=signif(k10, sigdig),
      k12=signif(k12, sigdig),
      k21=signif(k21, sigdig),
      Vss=signif(V1 + V2, sigdig),
      thalf_alpha=signif(log(2)/i1, sigdig),
      thalf_beta=signif(log(2)/i2, sigdig),
      alpha=signif(i1, sigdig),
      beta=signif(i2, sigdig),
      trueA=signif(c1, sigdig),
      trueB=signif(c2, sigdig),
      fracA=signif(c1*V1, sigdig),
      fracB=signif(c2*V1, sigdig),
      tau=tau,
      AUCinf=AUCinf,
      AUCtau=AUCtau,
      Cmax=Cmax,
      Cmaxss=Cmaxss,
      Tmax=Tmax,
      AUCinf_dose_normalized=AUCinf_dose_normalized,
      AUCtau_dose_normalized=AUCtau_dose_normalized,
      Cmax_dose_normalized=Cmax_dose_normalized,
      Cmaxss_dose_normalized=Cmaxss_dose_normalized,
      step=step,
      # Include the macro parameters
      CL=CL,
      V1=V1,
      V2=V2,
      Q2=Q2,
      ka=ka,
      tlag=tlag,
      dur=dur,
      tinf=tinf,
      dose=dose
    )

  o <- select_derived(out, type)

  maybe_warn_flipflop(ka = ka, k10 = k10)

  o
}

#' @rdname calc_derived
#' @examples
#' params <- calc_derived_3cpt(CL=29.4, V1=23.4, V2=114, V3=4614, Q2=270, Q3=73)
#' @export
calc_derived_3cpt <- function(CL, V1=NULL, V2, V3, Q2=NULL, Q3, V=NULL, Q=NULL, ka=NULL, dur=NULL, tinf=NULL, tlag=NULL, dose=NULL, tau=NULL, step=0.1, type="all", sigdig=5) {
  rowwise <- derived_rowwise(calc_derived_3cpt, mget(names(match.call())[-1]))
  if (!is.null(rowwise)) return(rowwise)

  if (!xor(is.null(V), is.null(V1))) {
    stop("Exactly one of V or V1 may be provided since they are considered synonyms.")
  } else if (!is.null(V)) {
    # Use the preferred synonym
    V1 <- V
  }
  if (!xor(is.null(Q), is.null(Q2))) {
    stop("Exactly one of Q or Q2 may be provided since they are considered synonyms.")
  } else if (!is.null(Q)) {
    # Use the preferred synonym
    Q2 <- Q
  }
  check_derived_inputs(
    positive = list(CL = CL, V1 = V1, V2 = V2, V3 = V3, Q2 = Q2, Q3 = Q3, ka = ka, dur = dur, tinf = tinf,
                    tau = tau, step = step),
    nonnegative = list(tlag = tlag, dose = dose)
  )

  ### variables
  k10   <- CL/V1
  k12   <- Q2/V1
  k21   <- Q2/V2
  k13   <- Q3/V1
  k31   <- Q3/V3

  a0 <- k10 * k21 * k31
  a1 <- (k10 * k31) + (k21 * k31) + (k21 * k13) + (k10 * k21) + (k31 * k12)
  a2 <- k10 + k12 + k13 + k21 + k31

  p   <- a1 - (a2 * a2 / 3)
  q   <- (2 * a2 * a2 * a2 / 27) - (a1 * a2 / 3) + a0
  r1  <- sqrt(-(p * p * p)/27)
  phi <- acos((-q/2)/r1)/3
  r2  <- 2 * exp(log(r1)/3)

  root1 <- -(cos(phi) * r2 - a2/3)
  root2 <- -(cos(phi + 2 * pi/3) * r2 - a2/3)
  root3 <- -(cos(phi + 4 * pi/3) * r2 - a2/3)

  i1 <- pmax(root1, root2, root3)
  # There is no pmedian function
  i2 <- mapply(FUN=function(...) stats::median(c(...)), root1, root2, root3)
  i3 <- pmin(root1, root2, root3)

  c1 <- (k21 - i1) * (k31 - i1) / (i1 - i2) / (i1 - i3) / V1
  c2 <- (k21 - i2) * (k31 - i2) / (i2 - i1) / (i2 - i3) / V1
  c3 <- (k21 - i3) * (k31 - i3) / (i3 - i2) / (i3 - i1) / V1

  AUCinf <- NULL
  AUCtau  <- NULL
  Cmax   <- NULL
  Cmaxss <- NULL
  Tmax   <- NULL
  AUCinf_dose_normalized <- NULL
  AUCtau_dose_normalized <- NULL
  Cmax_dose_normalized   <- NULL
  Cmaxss_dose_normalized <- NULL

  guess_tmax <- 0
  if(!is.null(ka)) guess_tmax <- 1/ka
  if(!is.null(dur)) guess_tmax <- dur
  if(!is.null(tinf)) guess_tmax <- tinf
  if(!is.null(tlag)) guess_tmax <- guess_tmax + tlag

  tend <- update_tend(tend = 24, guess_tmax = guess_tmax)

  if(!is.null(dose)) {
    message("dose present. NCA parameters will be estimated.")
    if(!is.null(tau)) message("tau present. Steady state NCA parameters will be estimated.")
    AUCinf <- dose/CL
    AUCinf_dose_normalized <- AUCinf/dose

    # IV bolus
    if(is.null(dur) & is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("3-compartment bolus detected.")
      peak <- find_peak(function(t) calc_sd_3cmt_linear_bolus(t=t, dose=dose, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_3cmt_linear_bolus(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }

    # IV infusion
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tinf)) {
      message("3-compartment infusion detected.")
      peak <- find_peak(function(t) calc_sd_3cmt_linear_infusion(t=t, dose=dose, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, tinf = tinf), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & !is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_3cmt_linear_infusion(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, tinf=tinf)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral first-order
    if(is.null(dur) & !is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("3-compartment first-order oral detected.")
      peak <- find_peak(function(t) calc_sd_3cmt_linear_oral_1(t=t, dose=dose, ka=ka, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & !is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_3cmt_linear_oral_1(tad = t, dose=dose, tau=tau, ka=ka, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral first-order lag
    if(is.null(dur) & !is.null(ka) & !is.null(tlag) & is.null(tinf)) {
      message("3-compartment first-order oral with lag time detected.")
      peak <- find_peak(function(t) calc_sd_3cmt_linear_oral_1_lag(t=t, dose=dose, ka=ka, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, tlag=tlag), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(is.null(dur) & !is.null(ka) & !is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_3cmt_linear_oral_1_lag(tad = t, dose=dose, tau=tau, ka=ka, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, tlag=tlag)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral zero-order
    if(!is.null(dur) & is.null(ka) & is.null(tlag) & is.null(tinf)) {
      message("3-compartment zero-order oral detected.")
      peak <- find_peak(function(t) calc_sd_3cmt_linear_oral_0(t=t, dose=dose, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, dur=dur), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(!is.null(dur) & is.null(ka) & is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_3cmt_linear_oral_0(tad = t, dose=dose, tau=tau, dur=dur, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
    # oral zero-order lag
    if(!is.null(dur) & is.null(ka) & !is.null(tlag) & is.null(tinf)) {
      message("3-compartment zero-order oral with lag time detected.")
      peak <- find_peak(function(t) calc_sd_3cmt_linear_oral_0_lag(t=t, dose=dose, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, tlag=tlag, dur=dur), tend = tend, step = step,
                        dose = dose, at = c(tinf, tlag, sum(tlag, dur)))
      Cmax <- peak$Cmax
      Tmax <- peak$Tmax
      Cmax_dose_normalized <- Cmax/dose
    }
    if(!is.null(dur) & is.null(ka) & !is.null(tlag) & !is.null(tau) & is.null(tinf)) {
      t <- seq_through(tau, step, at = ss_break_points(tau, tinf, tlag, dur))
      ss_conc <- function(t) calc_ss_3cmt_linear_oral_0_lag(tad = t, dose=dose, tau=tau, CL=CL, V1=V1, V2=V2, V3=V3, Q2=Q2, Q3=Q3, tlag=tlag, dur=dur)
      C <- quiet_flipflop(ss_conc(t))
      Cmaxss <- max(C)
      AUCtau  <- dose/CL
      AUCtau_dose_normalized <- AUCtau/dose
      Cmaxss_dose_normalized <- Cmaxss/dose
    }
  }


  out <-
    list(
      k10=signif(k10, sigdig),
      k12=signif(k12, sigdig),
      k21=signif(k21, sigdig),
      k13=signif(k13, sigdig),
      k31=signif(k31, sigdig),

      Vss=signif(V1 + V2 + V3, sigdig),

      thalf_alpha=signif(log(2)/i1, sigdig),
      thalf_beta =signif(log(2)/i2, sigdig),
      thalf_gamma=signif(log(2)/i3, sigdig),

      alpha=signif(i1, sigdig),
      beta =signif(i2, sigdig),
      gamma=signif(i3, sigdig),

      trueA=signif(c1, sigdig),
      trueB=signif(c2, sigdig),
      trueC=signif(c3, sigdig),

      fracA=signif(c1*V1, sigdig),
      fracB=signif(c2*V1, sigdig),
      fracC=signif(c3*V1, sigdig),
      tau=tau,
      AUCinf=AUCinf,
      AUCtau=AUCtau,
      Cmax=Cmax,
      Cmaxss=Cmaxss,
      Tmax=Tmax,
      AUCinf_dose_normalized=AUCinf_dose_normalized,
      AUCtau_dose_normalized=AUCtau_dose_normalized,
      Cmax_dose_normalized=Cmax_dose_normalized,
      Cmaxss_dose_normalized=Cmaxss_dose_normalized,
      step=step,
      # Include the macro parameters
      CL=CL,
      V1=V1,
      V2=V2,
      V3=V3,
      Q2=Q2,
      Q3=Q3,
      ka=ka,
      tlag=tlag,
      dur=dur,
      tinf=tinf,
      dose=dose
    )
  o <- select_derived(out, type)

  maybe_warn_flipflop(ka = ka, k10 = k10)

  o
}
