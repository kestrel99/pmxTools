# Build data/nhanes_peds.rda from the public NHANES files.
#
# Run from the package root:  Rscript data-raw/nhanes_peds.R
# Needs internet access and the haven and usethis packages.
#
# Children aged 2-17 whole years at the examination, from four NHANES
# releases. Children are kept when at least one of body weight or standing
# height was measured; a missing measure is kept as NA.

cycles <- data.frame(
  CYCLE = c("2013-14", "2015-16", "2017-18", "2021-23"),
  YEAR = c(2013, 2015, 2017, 2021),
  SUFFIX = c("H", "I", "J", "L")
)

read_nhanes_cycle <- function(cycle, year, suffix) {
  url <- paste0(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", year, "/DataFiles/"
  )
  demo <- haven::read_xpt(paste0(url, "DEMO_", suffix, ".XPT"))
  bmx <- haven::read_xpt(paste0(url, "BMX_", suffix, ".XPT"))

  demo_vars <- c("SEQN", "RIAGENDR", "RIDEXAGM", "WTMEC2YR")
  bmx_vars <- c("SEQN", "BMXWT", "BMXHT")
  missing_vars <- c(setdiff(demo_vars, names(demo)), setdiff(bmx_vars, names(bmx)))
  if (length(missing_vars) > 0) {
    stop("Cycle ", cycle, " lacks: ", paste(missing_vars, collapse = ", "))
  }

  d <- merge(demo[demo_vars], bmx[bmx_vars], by = "SEQN")
  positive <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x) | x <= 0] <- NA_real_
    x
  }

  out <- data.frame(
    CYCLE = cycle,
    SEQN = as.numeric(d$SEQN),
    AGE_MONTHS = as.numeric(d$RIDEXAGM),
    AGE = as.integer(floor(as.numeric(d$RIDEXAGM) / 12)),
    SEX = c("Male", "Female")[match(as.integer(d$RIAGENDR), 1:2)],
    WT = positive(d$BMXWT),
    HT = positive(d$BMXHT),
    MEC_WT = positive(d$WTMEC2YR),
    stringsAsFactors = FALSE
  )

  out[
    !is.na(out$AGE) & out$AGE %in% 2:17 &
      !is.na(out$SEX) &
      !is.na(out$MEC_WT) &
      (!is.na(out$WT) | !is.na(out$HT)),
  ]
}

nhanes_peds <- do.call(
  rbind,
  Map(read_nhanes_cycle, cycles$CYCLE, cycles$YEAR, cycles$SUFFIX)
)
nhanes_peds <- nhanes_peds[order(nhanes_peds$CYCLE, nhanes_peds$SEQN), ]
rownames(nhanes_peds) <- NULL

message("Rows: ", nrow(nhanes_peds))
message("Missing WT only: ", sum(is.na(nhanes_peds$WT)))
message("Missing HT only: ", sum(is.na(nhanes_peds$HT)))

usethis::use_data(nhanes_peds, overwrite = TRUE, compress = "xz")
