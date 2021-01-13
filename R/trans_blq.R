#' Estimate the lower limit of quantification (LLOQ) from a vector
#'
#' Nonnegative values are considered to be above the LLOQ. \code{NA} values are
#' ignored.
#' 
#' @param x The numeric vector to use for estimation of the LLOQ
#' @return The lowest, nonzero value from \code{x}.  If all are \code{NA} or
#'   zero, 1 is returned, and a warning is issued.
#' @export
#' @examples
#' estimate_lloq(c(NA, 0, 2, 5))
#' @family BLQ Transformation
estimate_lloq <- function(x) {
  choices <- x[!is.na(x) & x > 0]
  if (length(choices)) {
    lloq <- min(choices)
  } else {
    warning("No samples above the lloq, using 1")
    lloq <- 1
  }
  lloq
}

#' A transform for ggplot2 with data that may be below the lower limit of
#' quantification
#' 
#' If the \code{lloq} is not provided, it will be estimated from the data as the
#' minimum value above zero.
#' 
#' @param lloq The value of the lower limit of quantification as a numeric
#'   scalar
#' @param x (only used if \code{lloq} is missing), the data for \code{lloq}
#'   estimation.
#' @param multiplier When data are \code{< lloq}, they are replaced by
#'   \code{lloq*multiplier} for display.
#' @param lloq_text The text to use on the axis to indicate values \code{<
#'   lloq}.  It will be automatically set to \code{paste0("<", lloq)} if
#'   missing.
#' @return A "trans" object based on the \code{scales} package for BLQ data.
#' @family BLQ Transformation
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(data=data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
#'   geom_point()
#' 
#' ggplot(data=data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
#'   geom_point() +
#'   scale_x_continuous(trans=blq_trans(lloq=3))
#' 
#' ggplot(data=data.frame(x=1:10, y=1:10), aes(x=x, y=y)) +
#'   geom_point() +
#'   scale_x_continuous(trans=blq_log10_trans(lloq=3))
#' }
#' @export
#' @importFrom scales trans_new breaks_extended
blq_trans <- function(lloq, x, multiplier=0.5, lloq_text) {
  if (missing(lloq)) {
    lloq <- estimate_lloq(x)
  }
  
  breaks_blq_linear <- breaks_blq_general(lloq=lloq, breakfun=scales::breaks_extended)
  scales::trans_new(
    name=paste0("blq-", lloq),
    transform=ftrans_blq_linear(lloq=lloq, multiplier=multiplier),
    inverse=itrans_blq_linear(lloq=lloq),
    breaks=breaks_blq_linear,
    #minor_breaks=breaks_blq_minor,
    format=label_blq(lloq=lloq, lloq_text=lloq_text)
  )
}

#' @describeIn blq_trans Log-scale transformation with BLQ
#' @param base The base for the logarithm
#' @export
#' @importFrom scales breaks_log trans_new
blq_log_trans <- function(lloq, x, multiplier=0.5, base=10, lloq_text) {
  if (missing(lloq)) {
    lloq <- estimate_lloq(x)
  }
  
  breaks_blq_log <- breaks_blq_general(lloq=lloq, breakfun=scales::breaks_log, trans=log)
  scales::trans_new(
    name=paste0("blq-log10-", lloq),
    transform=ftrans_blq_log(lloq=lloq, multiplier=multiplier, base=base),
    inverse=itrans_blq_log(lloq=lloq, base=base),
    breaks=breaks_blq_log,
    #minor_breaks=breaks_blq_minor,
    format=label_blq(lloq=lloq, lloq_text=lloq_text)
  )
}

#' Forward transformation for linear BLQ data
#' 
#' For ggplot2 scales.
#' 
#' @inheritParams blq_trans
#' @return A function of \code{x} that replaces \code{x < lloq} with
#'   \code{lloq*multiplier}
#' @family BLQ Transformation
#' @export
ftrans_blq_linear <- function(lloq, multiplier) {
  force(lloq)
  force(multiplier)
  function(x) {
    x[!is.na(x) & x < lloq] <- lloq*multiplier
    x
  }
}
#' Inverse transformation for linear BLQ data
#' 
#' For ggplot2 scales.
#' 
#' @inheritParams blq_trans
#' @return A function of \code{x} that replaces \code{x < lloq} with \code{lloq}
#' @family BLQ Transformation
#' @export
itrans_blq_linear <- function(lloq) {
  force(lloq)
  function(x) {
    x[!is.na(x) & x < lloq] <- lloq
    x
  }
}

#' @describeIn ftrans_blq_linear Log-scale transformation
#' @inheritParams blq_trans
#' @export
ftrans_blq_log <- function(lloq, multiplier, base=10) {
  force(lloq)
  force(multiplier)
  force(base)
  function(x) {
    mask_blq <- !is.na(x) & x < lloq
    mask_above_lloq <- !is.na(x) & x >= lloq
    x[mask_blq] <- log(lloq*multiplier, base)
    x[mask_above_lloq] <- log(x[mask_above_lloq], base)
    x
  }
}
#' @describeIn itrans_blq_linear Log-scale inverse transform
#' @inheritParams blq_trans
#' @export
itrans_blq_log <- function(lloq, base) {
  force(lloq)
  force(base)
  function(x) {
    mask_blq <- !is.na(x) & (base^x) < lloq
    mask_above_lloq <- !is.na(x) & (base^x) >= lloq
    x[mask_blq] <- lloq
    x[mask_above_lloq] <- base^x[mask_above_lloq]
    x
  }
}

#' Generate breaks for measurements below the limit of quantification
#' 
#' Breaks that are \code{< lloq} are removed.  If the lowest break is removed if
#' it is too close to the lloq.
#' 
#' For ggplot2 scales.  This is not usually used directly.  See
#' \code{blq_trans()} and \code{blq_log10_trans()} for the functions that are
#' more commonly used.
#' 
#' @inheritParams blq_trans
#' @param breakfun The function used for normal scale breaks if the \code{lloq}
#'   were not present.
#' @param trans A parameter translation function (typically either
#'   \code{identity} for linear scale or \code{log} for log scale).
#' @param ... passed as \code{breakfun(n=n, ...)}
#' @return A function for calculating breaks with arguments \code{x} and
#'   \code{n}
#' @family BLQ Transformation
#' @examples
#' breaks_blq_general(lloq=3, breakfun=scales::breaks_extended)(1:100, n=5)
#' @export
breaks_blq_general <- function(lloq, breakfun, trans=identity, ...) {
  force(lloq)
  force(breakfun)
  force(trans)
  function(x, n=5) {
    if (all(is.na(x))) {
      # If all data are NA, then the normal way to create all-NA breaks should be used
      breakfun(n=n, ...)(x=x, n=n)
    } else if (sum(!is.na(x) & x >= lloq) == 0) {
      # All data are NA or BQL
      lloq
    } else {
      x_mask_blq <- !is.na(x) & x < lloq
      x_loq_censored <- x
      x_loq_censored[x_mask_blq] <- lloq
      # Use a standard function to generate normal breaks
      breaks_orig <- breakfun(n=n, ...)(x=x_loq_censored, n=n)
      breaks_has_blq <- any(breaks_orig < lloq*(1+sqrt(.Machine$double.eps)))
      # Drop breaks below the LLOQ
      breaks_final <- breaks_orig[breaks_orig >= lloq]
      if (length(breaks_final) == 0) {
        breaks_final <- lloq
      } else if (any(x_mask_blq) | breaks_has_blq) {
        # include the LLOQ in the breaks if any of the data or the suggested
        # breaks are BLQ.
        
        # If the lloq is added, it shouldn't be too close to the first break
        spacing <- diff(trans(breaks_orig))[1]
        # Drop breaks that are too close together
        if ((min(trans(breaks_final)) - (0.5*spacing)) < trans(lloq)) {
          breaks_final <- breaks_final[breaks_final > min(breaks_final)]
        }
        breaks_final <- c(lloq, breaks_final)
      }
      breaks_final
    }
  }
}

#' Label axes with censoring labels for BLQ
#' 
#' For ggplot2 scales.
#' 
#' @inheritParams blq_trans
#' @return A function of \code{x} which returns the formatted values.
#' @family BLQ Transformation
#' @export
label_blq <- function(lloq, lloq_text) {
  force(lloq)
  if (missing(lloq_text)) {
    lloq_text <- paste0("<", lloq)
  } else {
    force(lloq_text)
  }
  function(x, ...) {
    ret <- format(x, ..., trim = TRUE, justify = "left")
    ret[is.na(x)] <- NA
    # Protect from approximations of LLOQ when formatting by using sqrt(...)
    mask_lloq <- !is.na(x) & (x < (lloq*(1+sqrt(.Machine$double.eps))))
    if (any(mask_lloq)) {
      ret[mask_lloq] <- lloq_text
    }
    ret
  }
}
