#' Visualize PsN SCM output.
#'
#' \code{plot_scm} returns a visualization of a Perl-speaks-NONMEM (PsN, \url{https://uupharmacometrics.github.io/PsN/}) SCM (stepwise covariate modeling)
#' procedure. It depends on the presence of \code{scmlog.txt} and \code{short_scmlog.txt} files in the
#' specified directory.
#'
#' @param dir A PsN SCM folder (containing \code{scmlog.txt} and \code{short_scmlog.txt}).
#' @param startPhase Where to start collating the output; can be \code{"forward"} (the default) or \code{"backward"}.
#' @param defCol Default node outline color.
#' @param defFillCol Default node fill color.
#' @param defFontCol Default node font color.
#' @param defWidth Default node outline width.
#' @param fwdSuccessCol Node outline color for a model fit matching the forward inclusion criterion.
#' @param fwdSuccessFillCol Node fill color for a model fit matching the forward inclusion criterion.
#' @param fwdSuccessFontCol Node font color for a model fit matching the forward inclusion criterion.
#' @param fwdFailCol Node outline color for a model fit not matching the forward inclusion criterion.
#' @param fwdFailFillCol Node fill color for a model fit not matching the forward inclusion criterion.
#' @param fwdFailFontCol Node font color for a model fit not matching the forward inclusion criterion.
#' @param bwdSuccessCol Node outline color for a model fit matching the backward elimination criterion.
#' @param bwdSuccessFillCol Node fill color for a model fit matching the backward elimination criterion.
#' @param bwdSuccessFontCol Node font color for a model fit matching the backward elimination criterion.
#' @param bwdFailCol Node outline color for a model fit not matching the backward elimination criterion.
#' @param bwdFailFillCol Node fill color for a model fit not matching the backward elimination criterion.
#' @param bwdFailFontCol Node font color for a model fit not matching the backward elimination criterion.
#' @param fullFwdCol Node outline color for the full forward model (i.e. the final model before the backward elimination procedure in SCM).
#' @param fullFwdFillCol Node fill color for the full forward model (i.e. the final model before the backward elimination procedure in SCM).
#' @param fullFwdFontCol Node font color for the full forward model (i.e. the final model before the backward elimination procedure in SCM).
#' @param fullFwdWidth Node outline width for the full forward model (i.e. the final model before the backward elimination procedure in SCM).
#' @param finalCol Node outline color for the final reduced model (i.e. the final model reached after the backward elimination procedure in SCM).
#' @param finalFillCol Node fill color for the final reduced model (i.e. the final model reached after the backward elimination procedure in SCM).
#' @param finalFontCol Node font color for the final reduced model (i.e. the final model reached after the backward elimination procedure in SCM).
#' @param finalWidth Node outline width for the final reduced model (i.e. the final model reached after the backward elimination procedure in SCM).
#' @param nodeStyle Node style. A string containing a comma-separated list of options (which include "filled", "striped", "wedged", "diagonals" and "rounded"). See the GraphViz documentation for further details.
#' @param nodeShape Node shape. Options include "box" (the default), "oval", "diamond", "egg", "plaintext", "point", "square", "triangle" and many more. See the GraphViz documentation for further details.
#' @param fontname Font for nodes. Options depend heavily on the local system - see the GraphViz documentation for further details.
#' @param rankdir Direction of graph layout. Possible values are "TB" (the default), "LR", "BT", "RL", corresponding to directed graphs drawn from top to bottom, from left to right, from bottom to top, and from right to left, respectively. 
#' @param layout Graph layout. Possible values are "dot" (the default), "neato", "twopi", and "circo". Note that of these, "dot" is the easiest to interpret and the others may produce odd results. 
#' @param lookupDF A data frame containing a lookup table for node labels. By default, {plot_scm} will use the PSN model names. If a lookup table containing the fields `Model` and `Alias` is provided, model names in `Model` will be replaced in the output plots by mtaching labels in `Alias`.   
#' @param ... Additional parameters passed to the underlying \code{\link[data.tree]{SetNodeStyle}} and \code{\link[data.tree]{SetEdgeStyle}} functions, which in turn rely on \code{\link[DiagrammeR]{DiagrammeR}}.  
#' 
#' @return A \code{grViz} object.
#' 
#' @details This function parses PsN SCM output and displays it as a GraphViz graph (effectively, an HTML widget). It is built on \code{\link[data.tree]{plot.Node}} - please refer to doucmentation for this function for a more detailed overview of what is possible (a lot). For more specific details, see \url{http://rich-iannone.github.io/DiagrammeR/docs.html}. 
#' 
#' @seealso NONMEM (\url{http://www.iconplc.com/innovation/nonmem/})
#' @seealso GraphViz (\url{https://graphviz.org/Documentation.php})
#' 
#' @seealso Lindbom L, Ribbing J & Jonsson EN (2004). Perl-speaks-NONMEM (PsN) - A Perl module for NONMEM related programming. Computer Methods and Programs in Biomedicine, 75(2), 85-94. \url{https://doi.org/10.1016/j.cmpb.2003.11.003}
#' @seealso Lindbom L, Pihlgren P & Jonsson N (2005). PsN-Toolkit - A collection of computer intensive statistical methods for non-linear mixed effect modeling using NONMEM. Computer Methods and Programs in Biomedicine, 79(3), 241-257. \url{https://doi.org/10.1016/j.cmpb.2005.04.005}
#' @author Justin Wilkins, \email{justin.wilkins@@occams.com}
#' @family NONMEM reading
#' @examples
#' \dontrun{
#' scm <- plot_scm("E:/DrugX/ModelDevelopment/scm310")
#' }
#'
#' @export
#' @import data.tree DiagrammeR
#' @importFrom stringr str_extract str_split

plot_scm <- function(dir, startPhase="forward", 
                     fwdSuccessCol = "#66C2A5", fwdFailCol="black",
                     bwdSuccessCol = "#FC8D62", bwdFailCol="black", defCol="black",
                     fwdSuccessFillCol = "#B3E2CD", fwdFailFillCol="white",
                     bwdSuccessFillCol = "#FDCDAC", bwdFailFillCol="white", defFillCol="white",
                     fwdSuccessFontCol = "black", fwdFailFontCol="black", 
                     bwdSuccessFontCol = "black", bwdFailFontCol="black",
                     defFontCol="black",
                     fullFwdCol = "#8DA0CB", finalCol = "#E78AC3", 
                     fullFwdFillCol = "#CBD5E8", finalFillCol = "#F4CAE4", 
                     fullFwdFontCol = "black", finalFontCol = "black", 
                     fullFwdWidth = "2px", finalWidth = "2px", defWidth = "1px",
                     nodeStyle="filled,rounded", nodeShape="box", 
                     fontname="helvetica", 
                     rankdir="TB", 
                     layout="dot",
                     lookupDF=NULL,
                     ...
                     ) {
  
  if(!is.null(lookupDF)) {
    if(!("data.frame" %in% class(lookupDF))) {
      stop("lookupDF must be a data frame with columns `Model` and `Alias`.")
    }
    if(!("Model" %in% names(lookupDF))) {
      stop("lookupDF must be a data frame with columns `Model` and `Alias`.")
    }
    if(!("Alias" %in% names(lookupDF))) {
      stop("lookupDF must be a data frame with columns `Model` and `Alias`.")
    }
  } 
  
  Step <- NULL 
  scm <- read_scm(dir, startPhase = startPhase)
  
  frows <- nrow(scm$forward)
  brows <- nrow(scm$backward)
  
  if(is.null(frows)) frows <- 0
  if(is.null(brows)) brows <- 0
  
  if((frows>50) | (brows>50)) {
    warning(paste("This SCM has ", frows, " forward records and ", brows, " backward elimination records.\nThis plot will have many nodes and may not be easy to read.\n",sep=""))
  }
  
  #scmTree <- Node$new("Base model")
  
  ### forward
  
  if("forward" %in% names(scm)) {
    scm$forward$Parent <- "Base model"
    scm$forward$ModelOld <- scm$forward$Model
    scm$forward$Model <- paste(scm$forward$Model, as.character(scm$forward$Step), sep="_F")
    scm$forwardSummary$BestModelOld <- scm$forwardSummary$BestModel
    scm$forwardSummary$BestModel <- paste(scm$forwardSummary$BestModel, 1:nrow(scm$forwardSummary), sep="_F")
    
    for(f in 1:nrow(scm$forward)) {
      if(scm$forward$Step[f]>1) {
        scm$forward$Parent[f] <- scm$forwardSummary$BestModel[as.numeric(as.character(scm$forward$Step[f]))-1]
      }
    }
    
    ## lookup list
    
    if(!is.null(lookupDF)) {
      lookupDF$ModelOld <- lookupDF$Model
      scm$forward <- merge(scm$forward, lookupDF[,c("ModelOld","Alias")], all.x=T)
      scm$forward$Alias[is.na(scm$forward$Alias)] <- scm$forward$ModelOld[is.na(scm$forward$Alias)]
    } else {
      scm$forward$Alias <- scm$forward$ModelOld
    }
    
    scm$forward$pathString <- paste("Base model", scm$forward$Model, sep="/")
    
    for (n in 2:max(as.numeric(as.character(scm$forward$Step)))) {
      parent <- "Base model"
      for(i in 2:n) {
        parent <- paste(parent, unique(scm$forward$Parent[as.numeric(as.character(scm$forward$Step))==i]), sep="/")
      }
      scm$forward$pathString[as.numeric(as.character(scm$forward$Step))==n] <- paste(parent, unique(scm$forward$Model[as.numeric(as.character(scm$forward$Step))==n]), sep="/")
    }
    
    scm$forward$FullFwd <- 0
    scm$forward$FullFwd[scm$forward$NewOFV==min(scm$forward$NewOFV[scm$forward$Significant==1])] <- 1
    
    scm$forward$Final <- 0
    
    scm_fwd <- as.Node(scm$forward)
    scm_fwd$Drop <- 0
    scm_fwd$Goal <- 0
  
  }
    ### backward
    
    if("backward" %in% names(scm)) {
      scm$backward$Parent <- "Full forward model"
      
      if(as.numeric(as.character(scm$backward$Step))[1]>1) {
        scm$backward$Step <- as.numeric(as.character(scm$backward$Step)) - 1 
      }
      
      scm$backward$ModelOld <- scm$backward$Model
      scm$backward$Model <- paste(scm$backward$Model, as.numeric(as.character(scm$backward$Step)), sep="_B")
      scm$backwardSummary$BestModelOld <- scm$backwardSummary$BestModel
      scm$backwardSummary$BestModel <- paste(scm$backwardSummary$BestModel, 1:nrow(scm$backwardSummary), sep="_B")
      
      ## lookup list
      
      if(!is.null(lookupDF)) {
        lookupDF$ModelOld <- lookupDF$Model
        scm$backward <- merge(scm$backward, lookupDF[,c("ModelOld","Alias")], all.x=T)
        scm$backward$Alias[is.na(scm$backward$Alias)] <- scm$backward$ModelOld[is.na(scm$backward$Alias)]
      } else {
        scm$backward$Alias <- scm$backward$ModelOld
      }
      
      parent <- "Full forward model"
      lowcov <- ""
      for (f in 1:max(as.numeric(as.character(scm$backward$Step)))) {
        if (f ==1) {
          scm$backward$pathString[as.numeric(as.character(scm$backward$Step))==f] <- paste(parent, scm$backward$Model[as.numeric(as.character(scm$backward$Step))==f], sep="/")
        }
        st <- subset(scm$backward, as.numeric(as.character(Step))==f)
        lowcov <- paste(lowcov, st$Model[abs(st$Drop) == min(abs(st$Drop))], sep="/")
        if(f>1) {
          scm$backward$pathString[as.numeric(as.character(scm$backward$Step))==f] <- paste("Full forward model", lowcov, sep="/")
        }
      }
      
      scm$backward$FullFwd <- 0
      
      scm$backward$Final <- 0
      scm$backward$Final[scm$backward$NewOFV==max(scm$backward$NewOFV[scm$backward$Insignificant==1])] <- 1
      
      scm_bwd <- as.Node(scm$backward)
      scm_bwd$Drop <- 0
      scm_bwd$Goal <- 0
      
    }
  
  if(exists("scm_bwd") & exists("scm_fwd")) {
    scm_b <- data.frame(Step=c(scm$forward$Step,scm$backward$Step),
                        Model=c(scm$forward$Model,scm$backward$Model),
                        BaseOFV=c(scm$forward$BaseOFV,scm$backward$BaseOFV),
                        NewOFV=c(scm$forward$NewOFV,scm$backward$NewOFV),
                        Drop=c(scm$forward$Drop,scm$backward$Drop),
                        Goal=c(scm$forward$Goal,scm$backward$Goal),
                        PVal=c(scm$forward$PVal,scm$backward$PVal),
                        ModelOld=c(scm$forward$ModelOld,scm$backward$ModelOld),
                        Parent=c(scm$forward$Parent,scm$backward$Parent),
                        pathString=c(scm$forward$pathString,scm$backward$pathString),
                        Phase=c(scm$forward$Phase,scm$backward$Phase),
                        Significant=c(scm$forward$Significant,scm$backward$Insignificant),
                        FullFwd=c(scm$forward$FullFwd,scm$backward$FullFwd),
                        Final=c(scm$forward$Final,scm$backward$Final),
                        Alias=c(scm$forward$Alias,scm$backward$Alias))

    rootBackward <- scm_b$pathString[scm_b$FullFwd==1]
    scm_b$pathString[scm_b$Phase=="backward"] <- gsub("Full forward model", replacement = rootBackward, 
                                                      x = scm_b$pathString[scm_b$Phase=="backward"])
    
    scm_p <- as.Node(scm_b)
    scm_p$Drop = 0
    scm_p$Goal = 0
  }
  
  if(exists("scm_bwd") & !(exists("scm_fwd"))) {
    scm_p <- scm_bwd  
  }
  
  if(exists("scm_fwd") & !(exists("scm_bwd"))) {
    scm_p <- scm_fwd  
  }
  
  
  SetGraphStyle(scm_p, rankdir=rankdir, layout=layout)
  
  SetNodeStyle(scm_p, style = nodeStyle, shape = nodeShape, 
               color = function(node) {
                 col <- defCol
                 
                 if(node$isRoot) col <- defCol
                 
                 if(!node$isRoot) {
                   col <- defCol
                   if((node$Drop > node$Goal) & (node$Phase=="forward")) {
                     col <- fwdSuccessCol
                   }
                   if((node$Drop <= node$Goal) & (node$Phase=="forward")) {
                     col <- fwdFailCol
                   }
                   if((node$Drop > node$Goal) & (node$Phase=="backward")) {
                     col <- bwdSuccessCol
                   }
                   if((node$Drop <= node$Goal) & (node$Phase=="backward")) {
                     col <- bwdFailCol
                   }
                   if(node$FullFwd==1) col <- fullFwdCol
                   if(node$Final==1) col <- finalCol
                 }
                 
                 col
               },
               fillcolor = function(node) {
                 fill <- defFillCol
                 
                 if(node$isRoot) fill <- defFillCol
                 
                 if(!node$isRoot) {
                   fill <- defCol
                   if((node$Drop > node$Goal) & (node$Phase=="forward")) {
                     fill <- fwdSuccessFillCol
                   }
                   if((node$Drop <= node$Goal) & (node$Phase=="forward")) {
                     fill <- fwdFailFillCol
                   }
                   if((node$Drop > node$Goal) & (node$Phase=="backward")) {
                     fill <- bwdSuccessFillCol
                   }
                   if((node$Drop <= node$Goal) & (node$Phase=="backward")) {
                     fill <- bwdFailFillCol
                   }
                   if(node$FullFwd==1) fill <- fullFwdFillCol
                   if(node$Final==1) fill <- finalFillCol
                 }
                 
                 fill
               }, 
               fontcolor = function(node) {
                 fontcol <- defFontCol
                 
                 if(node$isRoot) fontcol <- defFontCol
                 
                 if(!node$isRoot) {
                   fontcol <- defCol
                   if((node$Drop > node$Goal) & (node$Phase=="forward")) {
                     fontcol <- fwdSuccessFontCol
                   }
                   if((node$Drop <= node$Goal) & (node$Phase=="forward")) {
                     fontcol <- fwdFailFontCol
                   }
                   if((node$Drop > node$Goal) & (node$Phase=="backward")) {
                     fontcol <- bwdSuccessFontCol
                   }
                   if((node$Drop <= node$Goal) & (node$Phase=="backward")) {
                     fontcol <- bwdFailFontCol
                   }
                   
                   if(node$FullFwd==1) fontcol <- fullFwdFontCol
                   if(node$Final==1) fontcol <- finalFontCol
                 }
                 
                 fontcol
               }, 
               fontname = fontname, tooltip = GetDefaultTooltip, 
               penwidth = function(node) {
                 pw <- defWidth
                 if(!node$isRoot) {
                   if(node$FullFwd==1) pw <- fullFwdWidth
                   if(node$Final==1) pw <- finalWidth
                 }
                 pw
               },
               label = function(node) {
                 if(!node$isRoot) {
                   if(node$Phase=="backward") {
                     paste(paste("Drop", node$Alias), paste0("+", sprintf("%.3f", abs(node$Drop))), sep="\n")
                   } else {
                     paste(paste("Add", node$Alias), paste0("-", sprintf("%.3f", abs(node$Drop))), sep="\n")
                   }
                 }
               }, ...)
  
  SetEdgeStyle(scm_p, ...)
  
  plot(scm_p)
  
}