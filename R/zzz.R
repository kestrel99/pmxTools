# .onLoad <- function(libname, pkgname) {
#     createPmxClasses()
# }

utils::globalVariables(c(
  "ID", "TIME", "EVENT",  # datamap
  "n_vals",               # cut_quantile::check_and_prepare_id
  "n_groups"              # cut_quantile::print_verbose_output
))
