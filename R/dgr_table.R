#' Generate a summary table of descriptive data for every individual in a dataset suitable for tabulation in a report.
#'
#' @param dat An input data frame, with one row per unique individual.
#' @param fields A vector of strings containing the names of the fields to be included in the summary table.
#' @param names A vector of strings containing descriptive names for the fields to be included in the summary table.
#' @param cutoff An integer defining the maximum number of unique values a variable should have to be considered categorical. Fields with more than this number of unique values are considered continuous for the purposes of the summary table (defaults to 7).
#' @param sig The number of significant digits summary values should have (defaults to 3).
#' @param by The field to use for grouping (a string). If not \code{NULL} (the default), the summary table will contain columns for each unique value of this field, as well as a column summarizing across all fields.
#' @param idvar The field in the dataset identifying each unique individual (defaults to "ID").
#' @param navars A vector containing values that are to be interpreted as missing (defaults to "-99" and "-999"). `NA` values are always considered to be missing.
#' 
#' @return A data frame containing a summary of all the fields listed in \code{fields}, for each individual in the dataset (the dataset should not contain duplicated individuals), conditioned on the field in \code{by}. Continuous values are summarized as median, mean, range and number of missing values. Categorical values are summarized as count and relative percentage. 
#' 
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#'  count_na(c(0,5,7,NA,3,3,NA))
#' }
#'
#' @export
#' @importFrom dplyr mutate_if 
#' @importFrom patchwork area 

# generate summary table
dgr_table <- function(dat, fields, names, cutoff=7, sig=3, by=NULL, idvar="ID", navars=c("-99","-999")) {
  
  # get proportions (for summary tables)
  ptable <- function(x, sig=3) signif(as.numeric(prop.table(table(x)))*100, sig)
  
  # summarise variable column (for summary tables)
  summ_field <- function(dat, field, name, cutoff=7, sig=3, by=NULL) {
    dat[[field]][dat[[field]] %in% navars] <- NA
    
    if(length(unique(dat[[field]])) > cutoff & is.numeric(dat[[field]])) {
      cont <- T
    } else {
      cont <- F
    }
    if(cont) {
      if(!is.null(by)) {
        gm <- c(tapply(dat[[field]], dat[[by]], mean, na.rm=T), mean(dat[[field]], na.rm=T, neg.rm=T))
        md <- c(tapply(dat[[field]], dat[[by]], median, na.rm=T), median(dat[[field]], na.rm=T))
        mx <- c(tapply(dat[[field]], dat[[by]], max, na.rm=T), max(dat[[field]], na.rm=T))
        mn <- c(tapply(dat[[field]], dat[[by]], min, na.rm=T), min(dat[[field]], na.rm=T))
        miss <- c(tapply(dat[[field]], dat[[by]], count_na), length(dat[[field]][is.na(dat[[field]])]))
        row <- t(as.matrix(c(name, paste(signif(md, sig), " (", signif(gm, sig), ")\n", 
                                         "[", signif(mn, sig), " ; ", signif(mx, sig), "] ", 
                                         "{", miss, "}", sep=""))))
      } else {
        row <- t(as.matrix(c(name, paste(signif(median(dat[[field]], na.rm=T), sig), 
                                         " (", signif(mean(dat[[field]], na.rm=T), sig), ")\n", 
                                         "[", signif(min(dat[[field]], na.rm=T), sig), " ; ", signif(max(dat[[field]], na.rm=T), sig), "] ", 
                                         "{", length(dat[[field]][is.na(dat[[field]])]), "}", sep=""))))
      }
    } else {
      dat[[field]] <- factor(dat[[field]])
      if(!is.null(by)) {
        dat$Summary <- dat[[by]]
        dat1 <- dat
        dat1$Summary <- "zzzAll"
        dat1 <- rbind(dat, dat1)
        
        tab  <- table(dat1[[field]], dat1$Summary)
        ptab <- cbind(matrix(unlist(tapply(dat[[field]], dat[[by]], ptable)), ncol=length(unique(dat[[by]]))), ptable(dat[[field]]))
        dimnames(ptab) <- dimnames(tab)
        
        mat <- matrix(paste(tab , " (", signif(ptab,3), "%)", sep=""), nrow=nrow(tab))
        mat <- cbind(paste(" - ", dimnames(tab)[[1]], sep=""), mat)
        row <- rbind(c(name,rep("", times=ncol(mat)-1)), mat)
      } else {
        tab  <- table(dat[[field]], useNA = "ifany")
        ptab <- tab/nrow(dat)*100
        
        mat <- matrix(paste(tab , " (", signif(ptab,3), "%)", sep=""), nrow=length(tab))
        mat <- cbind(paste(" - ", dimnames(tab)[[1]], sep=""), mat)
        row <- rbind(c(name,""), mat)
      }
    }
    row <- dplyr::mutate_if(tibble::as_tibble(row),
                     is.character,
                     stringr::str_replace_all, pattern = " \\{0\\}", replacement = "")
    row <- dplyr::mutate_if(tibble::as_tibble(row),
                     is.character,
                     stringr::str_replace_all, pattern = " - NA", replacement = " - Missing")
    row  
  }
  
  ## main function
  if(!is.null(by)) {
    out <- matrix(nrow=1, c("N", as.numeric(tapply(dat[[idvar]], dat[[by]], length)), nrow(dat)))
  } else {
    out <- matrix(nrow=1, c("N", nrow(dat)))
  }
  for (i in 1:length(fields)) {
    out <- rbind(out, summ_field(dat=dat, fields[i], names[i], cutoff=cutoff, sig=sig, by=by))
  }
  
  out <- data.frame(out, stringsAsFactors = T)
  if(!is.null(by)) {
    names(out) <- c("Variable",dimnames(tapply(dat[[idvar]], dat[[by]], length))[[1]],"Total")
  } else {
    names(out) <- c("Variable","Total")
  }
  out
}
