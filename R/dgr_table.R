#' @noRd
parse_fields <- function(fields, names = NULL) {
  if (!is.null(names)) {
    if (!is.null(names(fields))) {
      warning("Both 'fields' and 'names' provided. Using 'fields' names.")
      return(list(fields = setNames(unname(fields), names(fields)), labels = names))
    }
    if (length(fields) != length(names)) {
      stop("'fields' and 'names' must have the same length")
    }
    warning("The 'names' argument is deprecated. Use a named vector: dgr_table(dat, c(FIELD = 'Label'))")
    return(list(fields = setNames(fields, names), labels = names))
  }
  
  if (is.null(names(fields))) {
    stop("When 'fields' is unnamed, 'names' must be provided")
  }
  
  list(fields = fields, labels = unname(fields))
}

detect_field_type <- function(x, cutoff = 7) {
  n_unique <- length(unique(x[!is.na(x)]))
  is_numeric <- is.numeric(x)
  list(
    is_continuous = n_unique > cutoff && is_numeric,
    n_unique = n_unique
  )
}

convert_navars <- function(x, navars) {
  if (is.numeric(x) && !is.numeric(navars)) {
    converted <- suppressWarnings(as.numeric(navars))
    if (any(is.na(converted))) {
      warning("Could not convert some 'navars' values to numeric")
      navars
    } else {
      converted
    }
  } else {
    navars
  }
}

summarize_continuous <- function(dat, field, label, cutoff, sig, by, navars, mtype) {
  x <- dat[[field]]
  navars_conv <- convert_navars(x, navars)
  x[x %in% navars_conv] <- NA
  
  type_info <- detect_field_type(x, cutoff)
  if (!type_info$is_continuous) {
    return(summarize_categorical(dat, field, label, cutoff, sig, by, navars, mtype))
  }
  
  n_miss <- sum(is.na(x))
  
  if (!is.null(by)) {
    by_vals <- unique(dat[[by]])
    grp_stats <- dat %>%
      group_by(.data[[by]]) %>%
      summarize(
        gm = if (mtype == "geomean") {
          gm(.data[[field]], na.rm = TRUE)
        } else {
          mean(.data[[field]], na.rm = TRUE)
        },
        md = median(.data[[field]], na.rm = TRUE),
        mn = min(.data[[field]], na.rm = TRUE),
        mx = max(.data[[field]], na.rm = TRUE),
        n_miss = sum(is.na(.data[[field]])),
        .groups = "drop"
      )
    
    total_stats <- tibble(
      gm = if (mtype == "geomean") gm(x, na.rm = TRUE) else mean(x, na.rm = TRUE),
      md = median(x, na.rm = TRUE),
      mn = min(x, na.rm = TRUE),
      mx = max(x, na.rm = TRUE),
      n_miss = n_miss
    )
    total_stats[[by]] <- "Total"
    total_stats <- total_stats[, c(by, "gm", "md", "mn", "mx", "n_miss")]
    
    grp_stats <- bind_rows(grp_stats, total_stats)
    
    by_levels <- c(sort(by_vals), "Total")
    grp_stats[[by]] <- factor(grp_stats[[by]], levels = by_levels)
    grp_stats <- grp_stats[order(grp_stats[[by]]), ]
    
    result <- map(by_levels, ~{
      row <- grp_stats[grp_stats[[by]] == .x, ]
      if (nrow(row) == 0) {
        row <- tibble(md = NA, gm = NA, mn = NA, mx = NA, n_miss = NA)
      }
      format_cell(row$md, row$gm, row$mn, row$mx, row$n_miss, sig)
    })
    names(result) <- by_levels
    
  } else {
    result <- list(Total = format_cell(
      median(x, na.rm = TRUE),
      if (mtype == "geomean") gm(x, na.rm = TRUE) else mean(x, na.rm = TRUE),
      min(x, na.rm = TRUE),
      max(x, na.rm = TRUE),
      n_miss,
      sig
    ))
  }
  
  var_row <- c(list(Variable = label), result)
  as_tibble_row(c(Variable = label, result), .name_repair = "unique")
}

summarize_categorical <- function(dat, field, label, cutoff, sig, by, navars, mtype) {
  x <- dat[[field]]
  navars_conv <- convert_navars(x, navars)
  x[x %in% navars_conv] <- NA
  x <- factor(x, exclude = NULL)
  
  lvls <- levels(x)
  counts <- table(x, useNA = "ifany")
  
  format_cat_value <- function(cnt, total) {
    pct <- if (total > 0) cnt / total * 100 else 0
    paste0(cnt, " (", signif(pct, sig), "%)")
  }
  
  if (!is.null(by)) {
    by_vals <- sort(unique(dat[[by]]))
    by_levels <- c(by_vals, "Total")
    
    rows <- list()
    rows[[1]] <- c(Variable = label, setNames(rep("", length(by_levels)), by_levels))
    
    for (j in seq_along(lvls)) {
      lvl <- lvls[j]
      lvl_val <- if (is.na(lvl)) "Missing" else paste0("- ", lvl)
      row_vals <- map_chr(by_levels, ~{
        if (.x == "Total") {
          grp_x <- x
        } else {
          grp_x <- x[dat[[by]] == .x]
        }
        grp_counts <- table(grp_x, useNA = "ifany")
        grp_total <- sum(grp_counts)
        cnt <- grp_counts[as.character(lvl)]
        if (is.na(cnt)) cnt <- 0
        format_cat_value(cnt, grp_total)
      })
      rows[[j + 1]] <- c(Variable = lvl_val, setNames(row_vals, by_levels))
      names(rows)[j + 1] <- NULL
    }
    
    out <- tibble(Variable = map_chr(rows, 1))
    for (col in by_levels) {
      out[[col]] <- map_chr(rows, ~.x[[col]])
    }
    
  } else {
    total <- sum(counts)
    
    rows <- list()
    rows[[1]] <- c(Variable = label, Total = "")
    
    for (j in seq_along(lvls)) {
      lvl <- lvls[j]
      lvl_val <- if (is.na(lvl)) "Missing" else paste0("- ", lvl)
      rows[[j + 1]] <- c(Variable = lvl_val, Total = format_cat_value(counts[j], total))
    }
    
    out <- tibble(
      Variable = map_chr(rows, 1),
      Total = map_chr(rows, 2)
    )
  }
  
  out
}

format_cell <- function(md, gm, mn, mx, n_miss, sig) {
  paste0(
    signif(md, sig), " (", signif(gm, sig), ")\n",
    "[", signif(mn, sig), " ; ", signif(mx, sig), "]",
    if (n_miss > 0) paste0(" {", n_miss, "}") else ""
  )
}

as_tibble_row <- function(x, .name_repair = "unique") {
  n_cols <- length(x)
  n_rows <- max(lengths(x))
  
  mat <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
  colnames(mat) <- names(x)
  
  for (i in seq_along(x)) {
    val <- x[[i]]
    if (length(val) == 1) {
      mat[, i] <- val
    } else {
      mat[seq_along(val), i] <- val
    }
  }
  
  tibble::as_tibble(mat, .name_repair = .name_repair)
}

#' Generate a summary table of descriptive data for every individual in a dataset suitable for tabulation in a report.
#'
#' @param dat An input data frame, with one row per unique individual.
#' @param fields A named character vector where names are column names and values are
#'   descriptive labels. For backward compatibility, can also be an unnamed character
#'   vector with separate `names` argument.
#' @param names Deprecated. Use a named vector for `fields` instead.
#' @param cutoff An integer defining the maximum number of unique values a variable should have to be considered categorical. Fields with more than this number of unique values are considered continuous for the purposes of the summary table (defaults to 7).
#' @param sig The number of significant digits summary values should have (defaults to 3).
#' @param by The field to use for grouping (a string). If not \code{NULL} (the default), the summary table will contain columns for each unique value of this field, as well as a column summarizing across all fields.
#' @param idvar The field in the dataset identifying each unique individual (defaults to "ID").
#' @param navars A vector containing values that are to be interpreted as missing (defaults to "-99" and "-999"). `NA` values are always considered to be missing.
#' @param mtype The type of mean to apply; `geomean`, the geometric mean (default) or `mean`, the arithmetic mean.
#'
#' @return A data frame containing a summary of all the fields listed in \code{fields}, for each individual in the dataset (the dataset should not contain duplicated individuals), conditioned on the field in \code{by}. Continuous values are summarized as median, mean, range and number of missing values. Categorical values are summarized as count and relative percentage.
#'
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#'
#' @examples
#' dat <- data.frame(
#'   ID = 1:10,
#'   AGE = c(45, 52, 38, 61, 29, 55, 43, 67, 33, 48),
#'   SEX = c("M", "F", "M", "M", "F", "F", "M", "F", "M", "F")
#' )
#'
#' dgr_table(dat, c(AGE = "Age (years)", SEX = "Sex"))
#' dgr_table(dat, c(AGE = "Age (years)"), by = "SEX")
#'
#' @export
#' @importFrom patchwork area
#' @importFrom dplyr group_by summarize bind_rows mutate select rename relocate
#' @importFrom purrr map map_chr pmap map_dfr map_depth
#' @importFrom stringr str_replace_all
#' @importFrom tibble as_tibble tibble
dgr_table <- function(dat, fields, names = NULL, cutoff = 7, sig = 3,
                      by = NULL, idvar = "ID", navars = c("-99", "-999"),
                      mtype = "geomean") {
  
  field_info <- parse_fields(fields, names)
  field_cols <- names(field_info$fields)
  field_labels <- field_info$labels
  
  n_groups <- if (!is.null(by)) length(unique(dat[[by]])) + 1 else 1
  
  n_row <- tibble(Variable = "N")
  if (!is.null(by)) {
    by_vals <- sort(unique(dat[[by]]))
    by_levels <- c(by_vals, "Total")
    n_by <- tapply(dat[[idvar]], dat[[by]], length)
    n_total <- nrow(dat)
    n_row <- tibble(Variable = "N", !!!setNames(c(as.list(n_by), list(n_total)), by_levels))
  } else {
    n_row$Total <- as.character(nrow(dat))
  }
  
  n_row[] <- lapply(n_row, as.character)
  
  summaries <- pmap(list(field_cols, field_labels), ~{
    summarize_continuous(dat, ..1, ..2, cutoff, sig, by, navars, mtype)
  })
  
  summaries <- map(summaries, ~{
    .x[] <- lapply(.x, as.character)
    .x
  })
  
  out <- bind_rows(c(list(n_row), summaries))
  
  if (!is.null(by)) {
    col_names <- c("Variable", sort(unique(dat[[by]])), "Total")
  } else {
    col_names <- c("Variable", "Total")
  }
  
  out <- out[, col_names, drop = FALSE]
  
  out
}
