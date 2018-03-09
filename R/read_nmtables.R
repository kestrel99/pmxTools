#' Reads NONMEM output tables.
#'
#' @param tableFiles NONMEM table files to be read.
#' @param runNo Run number.
#' @param tabSuffix Table file suffix.
#' @param tableNames List of root table names, using the Xpose naming convention as the default.
#' @param quiet Flag for displaying intermediate output (defaults to \code{FALSE}).
#' @param ... Additional arguments.
#'
#' @return A data frame.
#'
#' @note Adapted from Xpose 4 (\url{https://CRAN.R-project.org/package=xpose4}).
#' @seealso Jonsson EN, Karlsson MO. Xpose--an S-PLUS based population pharmacokinetic/pharmacodynamic model building aid for NONMEM. Comput Methods Programs Biomed. 1999 Jan;58(1):51-64, \url{http://xpose4.sourceforge.net}
#' @author Justin Wilkins, Niclas Jonsson, Andrew Hooker
#'
#' @examples
#' \dontrun{
#' tables <- read_nmtables(315)
#' }
#'
#' @import utils
#' @export

read_nmtables <-
  function(tableFiles = NULL,
           runNo      = NULL,
           tabSuffix  = "",
           tableNames = c("sdtab","mutab","patab","catab",
                          "cotab","mytab","extra","xptab"),
           quiet = FALSE,
           ...) {

    ### based on read.nm.tables function from Xpose 4
    ### written by Andrew C. Hooker, Justin J. Wilkins, Mats O. Karlsson and E. Niclas Jonsson

    if (is.null(tableFiles)){
      if(is.null(runNo)) {
        stop("Run number must be specified if no table files are provided.\n")
        return(NULL)
      }
      tabFiles <- sapply(tableNames, paste, runNo, tabSuffix, sep="")
    } else {
      tabFiles <- tableFiles
    }

    ## Read in the table files
    totab      <- NULL
    totnam     <- NULL
    seen.files <- NULL
    filedim    <- NULL

    for(i in 1:length(tabFiles)) {
      filename <- tabFiles[i]
      if(!file.exists(filename)) {
        next
      } else {
        cat(paste("    Reading",filename,"\n"))

        ## Check which type of separator we have in our tables
        header.line = scan(file=filename, nlines=1, skip=1, what="character", sep="\n", quiet=T)
        sep.char = ""
        if(length(grep(",",header.line))!=0) sep.char = ","

        ## Check if we have unequal number of fields in the file
        ## used for multiple simulations
        fields.per.line      <- count.fields(filename)
        fields.in.first.line <- fields.per.line[1]
        fields.in.rest       <- fields.per.line[-1]
        if((length(unique(fields.in.rest))!=1) ||
           (all(fields.in.first.line==fields.in.rest))){
          if(!quiet) {
            cat(paste(filename," conatins varying numbers of fields.\n",sep=""))
            cat("This may be due to multiple TABLE and header rows \n")
            cat("caused by running multiple simulations in NONMEM (NSIM > 1).\n")
            cat("Will attempt to remove these rows. Please be patient...\n")
          }
          tmp   <- readLines(filename, n = -1)
          inds  <- grep("TABLE",tmp)
          if (length(inds)!=1){
            inds  <- inds[c(2:length(inds))]
            inds2 <- inds+1
            tempfile<- paste(filename,".xptmp",sep="")
            write.table(tmp[-c(inds,inds2)],file=tempfile,
                        row.names=FALSE,quote=FALSE)
            assign(paste("n.",filename,sep=""),read.table(tempfile,skip=2,header=T,sep=sep.char))
            unlink(tempfile)
          } else {
            assign(paste("n.",filename,sep=""),read.table(filename,skip=1,header=T,sep=sep.char))
          }
        } else {
          assign(paste("n.",filename,sep=""),read.table(filename,skip=1,header=T,sep=sep.char))
        }

        ## Remember the files seen
        ##if(is.null(seen.files)) {
        ##  seen.files <- paste("n.",filename,sep="")
        ##} else {
        seen.files <- c(seen.files,paste("n.",filename,sep=""))
        ##}
      }
    }

    ## Check if we found any table files

    if(any(is.null(seen.files))) {
      #if(tab.suffix!=sim.suffix) {
      cat("There don't seem to be any table files matching this row number (",
          runNo, ")!\n")
      return(NULL)
    }

    ## Check if the files have the same length
    for(nfile in seen.files) {
      if(is.null(filedim)) {
        filedim <- nrow(get(nfile))
      } else {
        filedim <- c(filedim,nrow(get(nfile)))
      }
    }

    file.df <- data.frame(seen.files=seen.files,filedim=filedim)
    lngths  <- sort(unique(file.df$filedim))

    if(length(lngths) !=1) {
      cat("\nThe table files associated with this run number (",runNo,
          ") appear\n")
      cat("to have different lengths.\n")
      cat("Please check your output, it is likely files have been modified or the $TABLE step has failed.\n")
      return(NULL)
    }


    ## Add the tables to totab and replicate the shorter ones to match
    ## the size of the longest one
    maxlngth <- max(file.df$filedim)

    ##singdef <-
    ##  c("id","idlab","idv","dv","pred","ipred","iwres","wres","res")

    for(ii in 1:nrow(file.df)) {
      filnam <- as.character(file.df[ii,"seen.files"])
      new.df <- get(filnam)
      sz     <- file.df[ii,"filedim"]
      rl     <- maxlngth/sz

      if(any(is.null(totab))) {
        totab <- new.df
      } else {
        totab <- cbind(totab,new.df)
      }

      totnam <- c(totnam,names(new.df))

      ## store parameters & covariates for Data.R & SData.R

      if(!is.na(pmatch("n.patab", filnam))){
        write(names(new.df), file=".patab.names.tmp")
      } else {
        if(!is.na(pmatch("n.catab", filnam))){
          write(names(new.df), file=".catab.names.tmp")
        } else {
          if(!is.na(pmatch("n.cotab", filnam))){
            write(names(new.df), file=".cotab.names.tmp")
          } else {
            if(!is.na(pmatch("n.sdtab", filnam))){
              write(names(new.df), file=".sdtab.names.tmp")
            }
          }
        }
      }


    }

    # cat(totnam, "\n")

    ## Delete all duplicates

    totab <- totab[, !duplicated(totnam)]
    return(totab)
  }

