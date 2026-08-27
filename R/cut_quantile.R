#' Create quantile-based bins for continuous variables
#'
#' @param dat A data frame containing the variables to bin.
#' @param var Variable(s) to bin: single name, character vector, or named list.
#'   Named list allows different quantile cuts per variable, e.g.
#'   `list(AGE = c(4, 3), WT = 4)` creates both quartiles and tertiles for AGE.
#' @param n_groups Number of quantile groups (2-5). Can be:
#'   - Single value applied to all variables
#'   - Vector applied in order
#'   - Named list for per-variable settings
#' @param missing_codes Values to treat as missing (defaults to c(-99, -999)).
#' @param blq_label Label for BLQ/zero values (defaults to "BLQ").
#' @param unit Unit for display, appended to interval labels. Can be single
#'   value or named list per variable.
#' @param id Optional subject ID column. If provided, quantile calculation
#'   uses 1 row per subject. Stops if a subject has conflicting values.
#' @param verbose If TRUE, prints summary and returns list with data and summary.
#'
#' @return
#' If `verbose = FALSE` (default): returns modified data frame invisibly.
#' If `verbose = TRUE`: returns a list with:
#'   - `data`: modified data frame
#'   - `summary`: tibble with cut details per variable
#'   - `skipped`: tibble of any skipped cuts (due to zero-range bins)
#'
#' Output columns added (for var = "CONC", n_groups = 4):
#'   - `CONCQ4Q`: numeric factor (1, 2, 3, 4)
#'   - `CONCQ4C`: character ("Q1", "Q2", "Q3", "Q4", "BLQ")
#'   - `CONCQ4CC`: continuous factor with intervals
#'
#' @examples
#' \dontrun{
#' # Single variable, quartiles
#' dat <- cut_quantile(dat, "AGE", n_groups = 4)
#'
#' # Multiple cuts on same variable
#' dat <- cut_quantile(dat, "AGE", n_groups = c(4, 3))
#'
#' # With longitudinal data
#' dat <- cut_quantile(dat, list(CONC = 4, AGE = 4), id = "SUBJID")
#'
#' # Verbose output
#' result <- cut_quantile(dat, list(CONC = c(4, 3), AGE = 4), verbose = TRUE)
#' }
#'
#' @import ggplot2
#' @importFrom dplyr mutate select group_by summarize ungroup filter pull left_join slice
#' @importFrom dplyr n_distinct bind_rows
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr starts_with
#' @importFrom purrr map_lgl
#' @export
cut_quantile <- function(dat,
                         var,
                         n_groups = 4,
                         missing_codes = c(-99, -999),
                         blq_label = "BLQ",
                         unit = NULL,
                         id = NULL,
                         verbose = FALSE) {

  cuts <- parse_cuts_input(var, n_groups)

  if (!is.null(id)) {
    dat <- check_and_prepare_id(dat, cuts$variable, id, missing_codes)
  }

  if (!is.null(id)) {
    dat_subj <- dat %>%
      group_by(.data[[id]]) %>%
      slice(1) %>%
      ungroup()
    n_id <- nrow(dat_subj)
  } else {
    dat_subj <- dat
    n_id <- NULL
  }

  summaries <- vector("list", nrow(cuts))
  skipped_list <- vector("list", nrow(cuts))
  cut_cols <- character(0)

  for (i in seq_len(nrow(cuts))) {
    v <- cuts$variable[i]
    ng <- cuts$n_groups[i]
    u <- cuts_unit(var, unit, v)
    mc <- cuts_missing_codes(var, missing_codes, v)

    result <- cut_single(
      dat_subj, v, ng, mc, u, blq_label,
      verbose = FALSE, n_id = n_id
    )

    dat_subj <- result$data

    if (!is.null(result$summary)) {
      summaries[[i]] <- result$summary
    }
    if (!is.null(result$skipped)) {
      skipped_list[[i]] <- result$skipped
    }
    if (is.null(result$skipped)) {
      cut_cols <- c(cut_cols, result$new_cols)
    }
  }

  if (!is.null(id)) {
    dat <- map_cuts_to_original(dat, dat_subj, id, cuts)
  } else {
    dat <- dat_subj
  }

  if (verbose) {
    summary_df <- bind_rows(summaries)
    skipped_df <- bind_rows(skipped_list)

    print_verbose_output(summary_df, skipped_df, n_id)

    return(invisible(list(
      data = dat,
      summary = summary_df,
      skipped = skipped_df
    )))
  }

  invisible(dat)
}

parse_cuts_input <- function(var, n_groups) {
  if (is.character(var) && length(var) == 1 && is.null(names(var))) {
    if (length(n_groups) > 1) {
      tibble::tibble(variable = rep(var, length(n_groups)), n_groups = as.integer(n_groups))
    } else {
      tibble::tibble(variable = var, n_groups = as.integer(n_groups))
    }
  } else if (is.character(var) && length(var) > 1 && is.null(names(var))) {
    if (length(n_groups) == 1) {
      tibble::tibble(variable = var, n_groups = rep(as.integer(n_groups), length(var)))
    } else {
      tibble::tibble(variable = var, n_groups = as.integer(n_groups))
    }
  } else if (is.list(var) && !is.null(names(var))) {
    var <- var[names(var) != ""]
    tibble::tibble(
      variable = rep(names(var), times = lengths(var)),
      n_groups = as.integer(unlist(var))
    )
  } else {
    stop("'var' must be a variable name, character vector, or named list")
  }
}

check_and_prepare_id <- function(dat, vars, id, missing_codes) {
  for (v in unique(vars)) {
    mc <- cuts_missing_codes(NULL, missing_codes, v)

    id_check <- dat %>%
      group_by(.data[[id]]) %>%
      summarize(
        n_vals = n_distinct(.data[[v]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_vals > 1)

    if (nrow(id_check) > 0) {
      stop(
        "Variable '", v, "' has conflicting values across records for ",
        nrow(id_check), " subject(s). ",
        "Use aggregate function or provide pre-processed data.\n",
        "Affected IDs: ", paste(head(id_check[[id]], 5), collapse = ", "),
        if (nrow(id_check) > 5) "..." else ""
      )
    }
  }
  dat
}

map_cuts_to_original <- function(dat, dat_subj, id, cuts) {
  subj_cols_to_join <- c(id)
  for (i in seq_len(nrow(cuts))) {
    v <- cuts$variable[i]
    ng <- cuts$n_groups[i]
    prefix <- get_prefix(ng)
    q_col <- paste0(v, prefix, ng, "Q")
    char_col <- paste0(v, prefix, ng, "C")
    subj_cols_to_join <- c(subj_cols_to_join, q_col, char_col)
  }
  subj_cols_to_join <- unique(subj_cols_to_join)

  subj_cuts <- dat_subj[, colnames(dat_subj) %in% subj_cols_to_join, drop = FALSE]

  dat <- left_join(dat, subj_cuts, by = id)

  dat
}

cuts_unit <- function(var, unit, v) {
  if (is.null(unit)) {
    return(NULL)
  }
  if (is.character(unit) && length(unit) == 1) {
    return(unit)
  }
  if (is.list(unit) && !is.null(names(unit))) {
    return(unit[[v]])
  }
  if (is.character(unit) && !is.null(names(unit))) {
    return(unit[v])
  }
  NULL
}

cuts_missing_codes <- function(var, missing_codes, v) {
  if (is.list(missing_codes) && !is.null(names(missing_codes))) {
    return(missing_codes[[v]] %||% missing_codes[[1]])
  }
  missing_codes
}

get_prefix <- function(n) {
  prefixes <- c("2" = "H", "3" = "T", "4" = "Q", "5" = "Q")
  prefixes[as.character(n)] %||% "Q"
}

cut_single <- function(dat,
                       var,
                       n_groups,
                       missing_codes,
                       unit,
                       blq_label,
                       verbose = FALSE,
                       n_id = NULL) {

  prefix <- get_prefix(n_groups)
  unit_suffix <- if (!is.null(unit) && nzchar(unit)) paste0(" ", unit) else ""

  q_col <- paste0(var, prefix, n_groups, "Q")
  c_col <- paste0(var, prefix, n_groups, "CC")
  char_col <- paste0(var, prefix, n_groups, "C")

  clean_vals <- dat[[var]]
  clean_vals[clean_vals %in% missing_codes] <- NA

  n_missing <- sum(is.na(clean_vals) | clean_vals %in% missing_codes)
  n_blq <- sum(!is.na(clean_vals) & !clean_vals %in% missing_codes & clean_vals <= 0, na.rm = TRUE)
  n_total <- nrow(dat)
  n_valid <- sum(!is.na(clean_vals) & clean_vals > 0 & !clean_vals %in% missing_codes, na.rm = TRUE)

  skipped <- NULL
  qtiles <- NULL

  if (n_valid == 0) {
    skipped <- tibble::tibble(
      var = var,
      n_groups = n_groups,
      reason = "No valid (non-BLQ, non-missing) values"
    )
  } else {
    valid_vals <- clean_vals[!is.na(clean_vals) & clean_vals > 0 & !clean_vals %in% missing_codes]
    qtiles <- stats::quantile(valid_vals, probs = seq(0, 1, length.out = n_groups + 1))

    zero_range <- which(diff(qtiles) == 0)
    if (length(zero_range) > 0) {
      skipped <- tibble::tibble(
        var = var,
        n_groups = n_groups,
        reason = paste0(
          "Zero-range bins at positions: ", paste(zero_range, collapse = ", "),
          " (likely due to ties in data)"
        )
      )
    }
  }

  if (!is.null(skipped)) {
    warning("Skipping cut '", var, "' with ", n_groups, " groups: ", skipped$reason)
    return(list(
      data = dat,
      summary = NULL,
      skipped = skipped,
      new_cols = character(0)
    ))
  }

  clean_vals <- ifelse(
    is.na(clean_vals) | clean_vals <= 0 | clean_vals %in% missing_codes,
    NA,
    clean_vals
  )

  binned <- cut(
    clean_vals,
    breaks = qtiles,
    labels = 1:n_groups,
    include.lowest = TRUE,
    right = FALSE
  )

  binned_char <- ifelse(
    is.na(binned),
    blq_label,
    paste0(prefix, n_groups, binned, unit_suffix)
  )

  dat[[q_col]] <- binned
  dat[[char_col]] <- binned_char

  bin_counts <- table(binned, useNA = "always")
  bin_df <- tibble(
    bin = names(bin_counts),
    n = as.integer(bin_counts)
  )
  bin_df$bin <- ifelse(is.na(bin_df$bin), blq_label, bin_df$bin)

  summary <- NULL
  if (verbose) {
    summary <- tibble(
      var = var,
      n_groups = n_groups,
      n_total = n_total,
      n_valid = n_valid,
      n_missing = n_missing,
      n_blq = n_blq,
      n_id = n_id,
      bins = list(bin_df)
    )
  }

  list(
    data = dat,
    summary = summary,
    skipped = NULL,
    new_cols = c(q_col, char_col)
  )
}

print_verbose_output <- function(summary_df, skipped_df, n_id) {
  message("\n", strrep("=", 60))
  message("cut_quantile summary")
  message(strrep("=", 60))

  if (nrow(summary_df) > 0) {
    summary_df <- summary_df %>%
      mutate(cuts = paste0(get_prefix(n_groups), n_groups), .before = 1) %>%
      select(-n_id)
    print(summary_df, n = Inf, width = Inf)
  } else {
    message("No cuts applied.")
  }

  if (nrow(skipped_df) > 0) {
    message(strrep("-", 60))
    message("Skipped cuts (zero-range bins):")
    print(skipped_df, n = Inf)
  }

  if (!is.null(n_id)) {
    message(strrep("-", 60))
    message("Note: n_id = ", n_id, " unique subjects used for quantile calculation")
  }

  message(strrep("-", 60))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
