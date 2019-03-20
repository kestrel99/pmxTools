#' Plot PsN SCM results.
#'
#' \code{plot_scm} returns a step-plot summary of a Perl-speaks-NONMEM (PsN, \url{https://uupharmacometrics.github.io/PsN/}) SCM (stepwise covariate modeling)
#' procedure. It depends on the presence of an \code{scmlog.txt} file in the
#' specified directory, and is inspired by the plot code provided with PsN 4.8.1. 
#'
#' @param dir A PsN SCM folder (containing \code{scmlog.txt}).
#' @param phase SCM phase. Can be \code{"both"} (the default), \code{"forward"} or \code{"backward"}.
#' @param lineCol  Line color. Default is '#902C10'.
#' @param lineType  Line type. Default is '1'.
#' @param lineSize  Line color. Default is '1'.
#' @param pointCol Point color. Default is '#902C10'.
#' @param pointShape Point shape. Default is '16'.
#' @param pointSize Point size. Default is '3'.
#' @param textCol Point color. Default is 'black'.
#' @param textSize Point color. Default is '6'.
#'
#' @return A \code{\link{ggplot2}} plot object.
#' 
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' 
#' @seealso Lindbom L, Ribbing J & Jonsson EN (2004). Perl-speaks-NONMEM (PsN) - A Perl module for NONMEM related programming. Computer Methods and Programs in Biomedicine, 75(2), 85-94. \url{https://doi.org/10.1016/j.cmpb.2003.11.003}
#' @seealso Lindbom L, Pihlgren P & Jonsson N (2005). PsN-Toolkit - A collection of computer intensive statistical methods for non-linear mixed effect modeling using NONMEM. Computer Methods and Programs in Biomedicine, 79(3), 241-257. \url{https://doi.org/10.1016/j.cmpb.2005.04.005}
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' 
#' @examples
#' \dontrun{
#' scm <- plot_scm("E:/DrugX/ModelDevelopment/scm310")
#' }
#'
#' @import utils ggplot2 ggrepel gridExtra
#' @export

plot_scm <- function(dir, phase="both", lineCol="#902C10", lineType=1, lineSize=1, 
                     pointCol="#902C10", pointShape=16, pointSize=3,
                     textCol="black", textSize=5) {

  s <- read_scm(dir = dir)
  
  ## forward
  
  s$forwardSummary$Step <- 1:nrow(s$forwardSummary)
  
  fdata <- s$forwardSummary[1:nrow(s$forwardSummary)-1,]
  nfsteps <- nrow(s$forwardSummary) -1
  
  sdata <- data.frame(step = c(0, rep(seq(1, max(s$forwardSummary$Step)-1, by=1), each=2), max(s$forwardSummary$Step)),
                      OFV  = rep(s$forwardSummary$BaseOFV, each=2))
  
  pdata <- data.frame(step = c(0:nfsteps)+0.5,
                      OFV = s$forwardSummary$BaseOFV,
                      dOFV = c(0, s$forwardSummary$dOFV[1:nfsteps]))
  pdata$dOFV_txt <- sprintf("%.3f", pdata$dOFV)
  
  fwdplot <- ggplot(sdata, aes_string("step", "OFV")) + 
    geom_line(col=lineCol, size=lineSize, linetype=lineType) +
    geom_point(data=pdata, aes_string("step", "OFV"), col=pointCol, size=pointSize) +
    geom_text_repel(data=pdata, aes_string("step", "OFV", label="dOFV_txt"), col=textCol, size=textSize) +
    scale_x_continuous("Included covariate", limits=c(0, nfsteps+1), breaks=(c(0:nfsteps)+0.5),
                       labels=c("None", as.character(s$forwardSummary$BestModel[1:nfsteps]))) +
    labs(title=paste("Forward inclusion step (p<", s$forwardP, ")", sep="")) +
    theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  ## backward
  
  a1 <- split(s$backward, s$backward$Step)
  a2 <- lapply(a1, function(x) {
    if(nrow(x) > 0) {
      out <- x[x$Drop==max(x$Drop),]
      out
    }
  })
  bdata <- do.call(rbind, a2)
  bdata <- bdata[bdata[["Insignificant"]]==1,]
  
  sdata <- data.frame(step = c(0, rep(seq(1, nrow(bdata), by=1), each=2), nrow(bdata)+1),
                      OFV  = c(rep(bdata$BaseOFV, each=2), rep(bdata$NewOFV[nrow(bdata)], times=2)))
  
  pdata <- data.frame(step = c(0:nrow(bdata))+0.5,
                      OFV = c(bdata$BaseOFV, bdata$NewOFV[nrow(bdata)]),
                      dOFV = c(0, bdata$Drop))
  pdata$dOFV_txt <- sprintf("%.3f", pdata$dOFV)
  
  bwdplot <- ggplot(sdata, aes_string("step", "OFV")) + 
    geom_line(col=lineCol, size=lineSize, linetype=lineType) +
    geom_point(data=pdata, aes_string("step", "OFV"), col=pointCol, size=pointSize) +
    geom_text_repel(data=pdata, aes_string("step", "OFV", label="dOFV_txt"), col=textCol, size=textSize, 
                    box.padding = 0.5, min.segment.length = 5) +
    scale_x_continuous("Excluded covariate", limits=c(0, nrow(bdata)+1), breaks=(c(0:nrow(bdata))+0.5),
                       labels=c("None", as.character(bdata$Model))) +
    labs(title=paste("Backward elimination step (p<", s$backwardP, ")", sep=""))+
    theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  ## plot
  
  if(phase=="forward") {
    print(fwdplot)
  }
  
  if(phase=="backward") {
    print(bwdplot)
  }
  
  if(phase=="both") {
    grid.arrange(fwdplot + coord_cartesian(ylim=c(min(s$forwardSummary$BaseOFV), max(s$forwardSummary$BaseOFV))), 
                 bwdplot + coord_cartesian(ylim=c(min(s$forwardSummary$BaseOFV), max(s$forwardSummary$BaseOFV))), nrow=1)
  }
}

