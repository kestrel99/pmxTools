#' Extract problem and estimation information from a NONMEM output object.
#'
#' @param x A NONMEM output object generated using \code{\link{read_nm}}.
#' @param sigdig Specifies the number of significant digits to be provided (default=6).
#' @param est.step Specifies which estimation step to return parameters from (default is the last).
#' 
#' @seealso NONMEM (\url{https://www.iconplc.com/solutions/technologies/nonmem})
#' 
#' @examples
#' \dontrun{
#'  nmOutput <- read_nm("run315.xml")
#'  probInfo <- get_probinfo(nmOutput)
#' }
#'
#' @export

get_probinfo <- function(x, sigdig=6, est.step=NULL){
  
  ## estimation options as a dataframe
  if(is.null(est.step)) {
    no_steps <- sum(stringr::str_count(names(x$nonmem$problem), "estimation"))
  } else {
    no_steps <- est.step
  }
  
  ind_est  <- match("estimation", names(x$nonmem$problem))-1+no_steps
  
  ## problem options as a dataframe - 61 columns with information
  prob_options <- data.frame(attributes(x$nonmem$problem$problem_options), stringsAsFactors = F)
  
  ## Get the most important information only
  ## Get the 'important' problem info only and rename columns to end user friendly names
  prob_info <- data.frame(nRecords = ifelse(sum(match(names(prob_options),"data_nrec"), na.rm = T) == 1, as.numeric(prob_options$data_nrec), NA),
                          nOBS = ifelse(sum(match(names(prob_options),"data_nobs"), na.rm = T) == 1, as.numeric(prob_options$data_nobs), NA),
                          nIND = ifelse(sum(match(names(prob_options),"data_nind"), na.rm = T) == 1, as.numeric(prob_options$data_nind), NA),
                          nTHETA = ifelse(sum(match(names(prob_options),"nthetat"), na.rm = T) == 1, as.numeric(prob_options$nthetat), NA),
                          nOMEGA = ifelse(sum(match(names(prob_options),"omega_blockdim"), na.rm = T) == 1, as.numeric(prob_options$omega_blockdim), NA),
                          nSIGMA = ifelse(sum(match(names(prob_options),"sigma_diagdim"), na.rm = T) == 1, as.numeric(prob_options$sigma_diagdim), NA),
                          COV_STEP_OMITTED = ifelse(prob_options$cov_omitted == 'no', 'NO', 'YES'),
                          COV_MATRIX = ifelse(sum(match(names(prob_options),"cov_matrix"), na.rm = T) == 1, prob_options$cov_matrix, NA))
  
  
  ## Info in parallel_est, table_series, estimation_options, estimation_information,
  ## monitor, termination_* (except termination_status), covariance*, covariance_*, invcovariance,
  ## correlation, eigenvalues is not used
  est_info <- data.frame(final_objective_function = signif(as.numeric(unlist(x$nonmem$problem[[ind_est]]$final_objective_function)), sigdig),
                         estimation_method = unlist(x$nonmem$problem[[ind_est]]$estimation_method),
                         estimation_title = unlist(x$nonmem$problem[[ind_est]]$estimation_title),
                         estimation_elapsed_time = unlist(x$nonmem$problem[[ind_est]]$estimation_elapsed_time),
                         termination_status = ifelse(as.numeric(unlist(x$nonmem$problem[[ind_est]]$termination_status))== 0,"SUCCESSFUL","TERMINATED"),
                         condition_number=NA)
  
  if(!is.null(x$nonmem$problem[[ind_est]]$eigenvalues)) {
    est_info$condition_number = ifelse(min(as.numeric(unlist(x$nonmem$problem[[ind_est]]$eigenvalues))) != 0, 
                                       signif(abs(max(as.numeric(unlist(x$nonmem$problem[[ind_est]]$eigenvalues))))/abs(min(as.numeric(unlist(x$nonmem$problem[[ind_est]]$eigenvalues)))), sigdig),
                                       Inf)   ## ratio of absolute max and min eigenvalues
  }
  
 
  ## license information
  licInfo <- unlist(x$nonmem$license_information)
  
  ## program information
  progInfo <- unlist(x$nonmem$program_information)
  
  ## Problem Title
  probTitle <- unlist(x$nonmem$problem$problem_title)
  
  ## Problem Start time, stop time and elapsed time
  probTimes <- data.frame(startTime = unlist(x$stop_datetime),
                          stopTime = unlist(x$stop_datetime),
                          elapsedTime = signif(as.numeric(unlist(x$total_cputime)), sigdig))
  
  ## nmtran - NOT USED
  # nmtran <- unlist(x$nmtran)
  
  ## control_stream - NOT USED
  # control_stream <- unlist(x$control_stream)
  
  out <- list()
  
  out$prob_title <- probTitle
  out$prog_info <- progInfo
  out$lic_info <- licInfo
  out$prob_times <- probTimes
  out$prob_info <- prob_info
  out$est_info <- est_info
  
  out
  
}