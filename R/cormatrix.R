#' A Function to Create a Correlation Matrix
#'
#' This function allows you to print the output of data frames and most R functions in an APA formatted table.
#' @param ivs List of variables from which to compute all pairwise correlations. Input must be a list (i.e., c("var1","var2")). Defaults to NULL.
#' @param data A data frame containing the observed variables to correlate. Defaults to NULL.
#' @param method Type of correlation coefficient to computed. Options include "pearson","kendall", or "spearman". Defaults to "pearson".  
#' @param digits Number of digits you want to round numbers too. Leave blank if you don't want to round numbers or if the Data Frame contains strings. Defaults to NULL.
#' @param output the format desired for output. Options include "html" or "r". Defaults to "r". 
#' @param columnnames Option to determine whether variable names are used as column headers (TRUE) or numbers are used instead (FALSE). Defaults to TRUE.
#' @param CI Option to include confidence intervals in the correlation matrix table. Defaults to FALSE. 
#' @param CIlevel If confidence intervals are requested, this sets the desired confidence interval. Defaults to .95
#' @param displayDF Option to include the degrees of freedom used in each test in the correlation matrix table. Defults to FALSE.
#' @param displayP Option to include the statistical significance (p-value) of each correlation in the matrix. Null hypothesis that there is no relationship between the variables (r=0). Defaults to FALSE.
#' @param showSig Option to include symbols beside correlation coefficient denoting the level of significance of each correlation. Levels include p < .1, p < .05, p < .1, p < .001. Defaults to FALSE.
#' @keywords Correlation
#' @export
#' @import broom
#' @import knitr
#' @examples
#' cormatrix(ivs = c("mpg","wt","hp"), data=mtcars, digits=3, output="html", columnnames = FALSE, CI=TRUE, CIlevel=.95, displayDF=TRUE, displayP=TRUE, showSig=TRUE)

cormatrix <- function(ivs=NULL, 
                      data="NULL",
                      method="pearson", 
                      output="r", 
                      means=FALSE, 
                      columnnames=TRUE, 
                      CI=FALSE, 
                      CIlevel=.95, 
                      displayDF=FALSE, 
                      displayP=FALSE, 
                      showSig=TRUE, 
                      digits=3){
  
  if (is.null(data)) stop('No dataframe provided')
  if (is.null(ivs)) stop('No variables provided')
  if (length(ivs) < 2) stop('Not enough variables provided')
  
  for (varName in ivs){
    data[[varName]] <- as.numeric(data[[varName]])
  }
  
  
  ##Define matrices used##
  cor.matrix <- matrix(nrow=length(ivs),ncol=length(ivs))
  p <-  matrix(nrow=length(ivs),ncol=length(ivs))
  n <-  matrix(nrow=length(ivs),ncol=length(ivs))
  conf <-  matrix(nrow=length(ivs),ncol=length(ivs))
  
  ##Conduct cor.test for each combination of variables in IVS and save output in specific matrices##
  for (i in 1:length(ivs)){
    for (w in 1:length(ivs)){
      cor.matrix[i,w] <- cor.test(data[[as.character(ivs[i])]],data[[as.character(ivs[w])]], method=method)$estimate
      p[i,w] <- round(cor.test(data[[as.character(ivs[i])]],data[[as.character(ivs[w])]], method=method)$p.value,digits)
      if (method=="pearson"){
        n[i,w] <- cor.test(data[[as.character(ivs[i])]],data[[as.character(ivs[w])]], method=method)$parameter
        conf[i,w] <- paste("[",paste(round(cor.test(data[[as.character(ivs[i])]],data[[as.character(ivs[w])]], method=method,conf.level = CIlevel)$conf.int[1:2],digits),collapse=", "),"]",sep="")
      }
      if(i==w){
        
        cor.matrix[i,w] <- "-"
        p[i,w] <- ""
        n[i,w] <- ""
        conf[i,w] <- ""
      
      }else if (p[i,w]<=.1 && p[i,w]>.05 && i!=w && showSig==TRUE){
        est <- as.character(round(as.numeric(cor.matrix[i,w]),digits))
        cor.matrix[i,w] <- paste(est,"†", sep="")
      }else if(p[i,w]<=.05 && p[i,w]>.01 && i!=w && showSig==TRUE){
        cor.matrix[i,w] <- paste(round(as.numeric(cor.matrix[i,w]),digits),"*",sep="")
      }else if(p[i,w]<=.01 && p[i,w]>.001 && i!=w && showSig==TRUE){
        cor.matrix[i,w] <- paste(round(as.numeric(cor.matrix[i,w]),digits),"**",sep="")
      }else if(p[i,w]<=.001 && i!=w && showSig==TRUE){
        cor.matrix[i,w] <- paste(round(as.numeric(cor.matrix[i,w]),digits),"***",sep="")
      }else{
        cor.matrix[i,w] <- paste(round(as.numeric(cor.matrix[i,w]),digits),"", sep="")
      }  
    }
  }
  
  ##Store method name for display##
  method.name <- cor.test(data[[as.character(ivs[i])]],data[[as.character(ivs[w])]], method=method)$method
  
  ##Add rows and set rownames based on output requested CI=TRUE/FALSE, displayDF=TRUE/FALSE, displayP=TRUE/FALSE##
  rnames <- ""
  cor.matrix2 <- matrix(nrow=0,ncol=length(ivs))
  extrarow <-0
  for (i in 1:length(ivs)){
    cor.matrix2 <- rbind(cor.matrix2,cor.matrix[i,])
    rnames <- cbind(rnames,ivs[i])
    if(CI==TRUE && method=="pearson"){
      extrarow <- extrarow + 1
      cor.matrix2 <- rbind(cor.matrix2,conf[,i])
      rnames <- cbind(rnames,"")
    }
    if(displayDF==TRUE && method=="pearson"){
      extrarow <- extrarow + 1
      cor.matrix2 <- rbind(cor.matrix2,n[,i])
      rnames <- cbind(rnames,"")
    }
    if (displayP==TRUE){
      extrarow <- extrarow + 1
      cor.matrix2 <- rbind(cor.matrix2,p[,i])
      rnames <- cbind(rnames,"")
    }
    rnames2 <- rep("",((length(ivs)*(extrarow+1))+1))
  }
  
  ##Set column names to IVs or numbers according to columnnames=TRUE/FALSE##
  if(columnnames==TRUE){
    
    cnames <- ivs
    colnames(cor.matrix2) <- cnames
    cor.matrix2 <- rbind(as.character(cnames),cor.matrix2)
    cor.matrix2 <- cbind(as.character(rnames),cor.matrix2)
    
  }else{
    
    cnames=1:length(ivs)
    for (i in 1:length(ivs)){
      rnames <- replace(rnames,rnames==as.character(ivs[i]), paste(as.character(i),". ",as.character(ivs[i]), sep=""))
    }
    cor.matrix2 <- rbind(as.character(cnames),cor.matrix2)
    cor.matrix2 <- cbind(as.character(rnames),cor.matrix2)
    
  }
  
  ##Print output in r or HTML##
  if(output=="r"){
    
    cnames <- rep("",(length(ivs)+1))
    cat("\nCorrelation Matrix for Variables: ", paste(as.character(ivs),collapse=", "), "\n\n")
    prmatrix(cor.matrix2, right=TRUE, quote=FALSE,rowlab=rnames2, collab=cnames)
    if(CI==TRUE && method=="pearson"){
      cilevelprint <- paste("Confidence Level = ",as.character(CIlevel),sep="")
    }else{
      cilevelprint <- ""
    }
    if(showSig==TRUE){
      sigprint <- "\n† p < .1, * p < .05, ** p < .01, *** p < .001"
    }else{
      sigprint <- ""
    }
    cat("note: ", method.name, ". ", cilevelprint,". ", sigprint,sep="")
    
  }else if(output=="html"){
    
    if(CI==TRUE && method=="pearson"){
      cilevelprint <- paste("Confidence Level= ",as.character(CIlevel),".",sep="")
    }else{
      cilevelprint <- ""
    }
    if(showSig==TRUE){
      sigprint <- "<br/>† <i>p</i> < .1, &#42; <i>p</i> < .05, &#42;&#42; <i>p</i> < .01, &#42;&#42;&#42; <i>p</i> < .001."
    }else{
      sigprint <- ""
    }
   
    cnames <- cor.matrix2[1,]
  
    note <- paste(method.name, ".&nbsp;", cilevelprint, sigprint,sep="")
    title <- paste("Correlation Matrix for Variables: ", paste(as.character(ivs),collapse=", "))
    APAhtmlTable(cor.matrix2[-1,], digits=digits, genNote=note,tableTitle=title,header=cnames)
  }
}