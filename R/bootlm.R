library(boot)
library(broom)
library(lm.beta)
library(plyr)
library(knitr)

bootReg <- function(formula,data,indices){
  d <- data[indices,]
  fit <- lm(formula,data=d)
  return(coef(fit))
}

semicorR <- function(x=NULL, y=NULL, cov=NULL, data=NULL){
  lm1.formula <- paste(x,"~",paste(cov,collapse="+"))
  yvar <- as.character(y)
  resX <- residuals(lm(as.formula(lm1.formula), data=data, na.action=na.exclude))
  s = as.matrix(cbind(data[y],resX))
  return(cor(s[,1], s[,2], use="complete.obs"))
}

round_df <- function(df=NULL, digits=3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

bootlm <- function(model = NULL, 
                   output = "r", 
                   data = NULL, 
                   boot = 1000, 
                   boottype="perc", 
                   digits=3, 
                   conf=0.95,
                   semicor=FALSE) {
  
  if (is.null(data)) stop('No dataframe provided')
  #if (is.null(y)) stop('No dependent variable provided')
  if (is.null(model)) stop('No independent variables provided')
  
  ##Bootstrapped Linear Model
  lm.model <- as.formula(model)
  bootResults <- boot(statistic=bootReg, formula=lm.model, data=data, R=boot)
  # if (semicor==TRUE) {
  #   bootResults.partcor <- boot(statistic=bootCor, formula=lm.model, data=data, R=boot)
  # }
  
  lm.model.object <- lm(lm.model,data = data)
  lm.summary <- summary(lm.model.object)
  totSummary.tbl <- tidy(lm.summary)
  
  lm.D9.beta <- lm.beta(lm.model.object)
  beta <- coef(lm.D9.beta, standardized = TRUE)
  
  vars <- all.vars(formula(lm.model.object))
  ci <- matrix(nrow =nrow(totSummary.tbl), ncol=5)
  partcor.ci <- matrix(nrow =nrow(totSummary.tbl), ncol=5)
  partcor <- vector()
  boot.se <- vector()
  boot.bias <- vector()
  boot.p <- vector()
  null.dist <- vector()
  #partcor[1,1] <- NA
  #colnames(partcor) <- "Part"
  for(i in 1:nrow(totSummary.tbl)){
    if (boottype=="perc"){
      ci[i,] <- boot.ci(bootResults, type=boottype, index=i, conf=conf)$percent
      if(semicor==TRUE){
        #print(boot.ci(bootResults.partcor, type=boottype, index=i, conf=conf)$percent)
        #partcor.ci[i,] <- boot.ci(bootResults.partcor, type=boottype, index=i, conf=conf)$percent
      }
    }else if(boottype=="bca"){
      ci[i,] <- boot.ci(bootResults, type=boottype, index=i, conf=conf)$bca
      if(semicor==TRUE){
        #partcor.ci[i,] <- boot.ci(bootResults.partcor, type=boottype, index=i, conf=conf)$percent
      }
    }
    boot.se[i] <- sqrt(var(bootResults$t[,i]))
    boot.bias[i] <- mean(bootResults$t[,i]) - totSummary.tbl$estimate[i]
    
    # The bootstrap p-value can then be approximated by
    null.dist <- bootResults$t[,i] - mean(bootResults$t[,i])
    bootP <- (1+sum(abs(null.dist) > abs(bootResults$t0[i])))/(1+bootResults$R)
    boot.p[i] <- bootP
    
    
    if(semicor==TRUE){  
      if(i<nrow(totSummary.tbl)){
        w <- i+1
        #partcor[w] <- semicorR(x=vars[w],y=vars[1],cov=vars[-c(1,w)], data=data)
      }
    }
  }
  
  
  #Remove extra information
  ci<-ci[,-c(1,2,3)] #delete columns 6, 7 and 8
  #partcor.ci<-partcor.ci[,-c(1,2,3)] #delete columns 6, 7 and 8
  
  
  #add bias and boot se and then reorder columns and then add estimate confididence interval
  totSummary.tbl <- cbind(totSummary.tbl,boot.bias,boot.se,beta)
  totSummary.tbl <- totSummary.tbl[c("term","estimate","beta","std.error","boot.bias","boot.se","statistic","p.value")]
  totSummary.tbl <- cbind(totSummary.tbl,ci)
  #totSummary.tbl <- cbind(totSummary.tbl,ci,boot.p)
  
  #Rename columns
  totSummary.tbl <- plyr::rename(totSummary.tbl, c("term"=" ","beta"="Beta","std.error"="SE","boot.bias"="Bias","boot.se"="Boot SE","statistic"="t","p.value"="p","1"="Boot LLCI","2"="Boot ULCI"))
  #totSummary.tbl <- plyr::rename(totSummary.tbl, c("term"=" ","beta"="Beta","std.error"="SE","boot.bias"="Bias","boot.se"="Boot SE","statistic"="t","p.value"="p","1"="Boot LLCI","2"="Boot ULCI","boot.p"="Bootstrap p"))
  
  # if(semicor==TRUE){
  #   #add , semi-partial correlation and it's confidence interval 
  #   #partcor.ci[1,1] <- NA
  #   #partcor.ci[1,2] <- NA
  #   totSummary.tbl <- cbind(totSummary.tbl,partcor)
  #   #totSummary.tbl <- cbind(totSummary.tbl,partcor.ci)
  #   
  #   #totSummary.tbl <- rename(totSummary.tbl, c("partcor"="semipartial r","1"= "semipartial LLCI","2"="semipartial ULCI"))
  #   totSummary.tbl <- plyr::rename(totSummary.tbl, c("partcor"="semipartial r"))
  # }
  
  
  
  if (boottype=="perc"){
    print.boottype <- "Percentile"
  }else if (boottype=="bca"){
    print.boottype <- "BCa"
  }
  
  if(output=="r"){
    
    #round results to 4 digits
    totSummary.tbl <- totSummary.tbl[-1]
    totSummary.tbl <- round_df(totSummary.tbl,digits)
    
    #Print Header
    cat(paste("Bootstrapped Regression Model \n \n Model = ", as.character(model), "\n Outcome = ", vars[1],"\n IVs = ",paste(vars[2:length(vars)], collapse = ", "),"\n Observations = ",nrow(data),"\n Bootstrap Type = ", print.boottype,"\n Bootstrap Samples = ",boot,"\n Confidence Interval = ", as.character(conf),  "\n \n",sep=""))
    
    #Print output
    print(totSummary.tbl, na.print = "")
    
    #Print Model Parameters
    f.sig = pf(lm.summary$fstatistic[1], lm.summary$fstatistic[2],lm.summary$fstatistic[3], lower.tail = FALSE)
    model.params <- paste("F(",lm.summary$fstatistic[2],",",lm.summary$fstatistic[3],") = ",round(lm.summary$fstatistic[1],3), ", p= ",round(f.sig,3) , ", R^2: ",round(lm.summary$r.squared,3), ", Adjusted R^2: ",round(lm.summary$adj.r.squared,3),sep="")
    cat("\n",model.params,"\n\n")
  }else if(output=="html"){
    
    totSummary.tbl <- round_df(totSummary.tbl,digits)
    
    #Model Parameters
    f.sig = pf(lm.summary$fstatistic[1], lm.summary$fstatistic[2],lm.summary$fstatistic[3], lower.tail = FALSE)
    model.params <- paste("F(",lm.summary$fstatistic[2],",",lm.summary$fstatistic[3],") = ",round(lm.summary$fstatistic[1],3), ", p= ",round(f.sig,3) , ", R^2: ",round(lm.summary$r.squared,3), ", Adjusted R^2: ",round(lm.summary$adj.r.squared,3),sep="")
    
    
    #Print Header
    header <- paste("<p><b>Bootstrapped Regression Model</b> <br /><br /> <b>Model:</b> ", as.character(model), "<br /> <b>Outcome =</b> ", vars[1],"<br /> <b>IVs =</b> ",paste(vars[2:length(vars)], collapse = ", "),"<br /> <b>Observations=</b> ",nrow(data),"<br /> <b>Bootstrap Type:</b> ", print.boottype,"<br /> <b>Bootstrap Samples=</b> ",boot,"<br /> <b>Confidence Interval =</b> ", conf, "<br /></p>",sep="")
    cat(header)
    
    #Print output
    APAhtmlTable(totSummary.tbl,round = FALSE, tidy=FALSE,tableTitle="Bootstrapped Regression results",genNote=model.params)
  }
}

APAhtmlTable <- function(df=NULL, 
                         round=TRUE, 
                         tidy=TRUE, 
                         digits=3, 
                         genNote=NULL, 
                         font="Times New Roman", 
                         fontSize=12, tableNum=NULL, 
                         tableTitle=NULL, 
                         headerBold=FALSE, 
                         headerAlign="center", 
                         leftColAlign="left", 
                         contentAlign="center", 
                         header=NULL){
  
  # for (df1 in df){
  #   #df1 <- tidy(df1)
  #   print(df1)
  # }
  #Set font options
  fontWeight <- ifelse(headerBold==TRUE,"font-weight:bold;","font-weight:normal;")
  if ((!headerAlign=="left") & (!headerAlign=="right") & (!headerAlign=="center")) stop("Invalid headerAlign set. Options include 'left', 'center', 'right'.")
  if ((!leftColAlign=="left") & (!leftColAlign=="right") & (!leftColAlign=="center")) stop("Invalid leftColAlign set. Options include 'left', 'center', 'right'.")
  if ((!contentAlign=="left") & (!contentAlign=="right") & (!contentAlign=="center")) stop("Invalid contentAlign set. Options include 'left', 'center', 'right'.")
  
  #Check header length
  #if((!is.null(header)) && (!length(header)==(ncol(df)+1))) stop("Number of header columns provided does not match number of columns in dataframe")
  
  #Convert object to dataFrame
  if(tidy==TRUE){
    df <- tidy(df)
  }
  
  #Round digits in data frame
  if(round==TRUE){
    df <- round_df(df,digits)
  }
  
  #Convert dataframe to characters and replace NAs with blanks
  df[, ] <- lapply(df[, ], as.character)
  df[is.na(df)] <- ""
  
  #Populate table header with column names
  if (is.null(header)){
    toprow.names <- as.list(colnames(df))
    toprow <- paste(as.character(toprow.names),collapse = paste("</th><th style='border-bottom:solid 1px #000000;border-top:solid 1px #000000;padding-right: 10px;text-align: ", headerAlign,";", as.character(fontWeight),"'>",sep=""),sep="")
  }else{
    toprow.names <- as.list(header)
    toprow <- paste(as.character(toprow.names),collapse = paste("</th><th style='border-bottom:solid 1px #000000;border-top:solid 1px #000000;padding-right: 10px;text-align: ", headerAlign,";", as.character(fontWeight),"'>",sep=""),sep="")
  }
  
  #Add start table tag
  tbl.summary <- paste("<br/><table style='border-collapse: collapse;font-family: ", font ,";font-size:", fontSize ,"pt'>",sep="")
  
  #Add Table Number if not Null
  if(!is.null(tableNum)){
    tbl.summary <- paste(tbl.summary, "<tr><td colspan=", ncol(df) ," style='line-height: 2;'>Table", as.character(tableNum),"</td></tr>",sep="")
  }
  
  #Add Table Title if not Null
  if(!is.null(tableTitle)){
    tbl.summary <- paste(tbl.summary, "<tr><td colspan='", ncol(df) ,"', style='font-style:italic;line-height: 2;'>", as.character(tableTitle),"</td></tr>",sep="")
  }
  
  #ADD header row
  tbl.summary <- paste(tbl.summary, "<tr><th style='border-bottom:solid 1px #000000;border-top:solid 1px #000000;padding-right: 10px;text-align: ", leftColAlign,";", fontWeight,"'>",toprow,"</th></tr>")
  
  #Add content rows
  for (i in 1:nrow(df)){
    
    cell.style <- ifelse(i==nrow(df),paste("padding-right: 10px;text-align: ", contentAlign,";border-bottom:solid 1px #000000;"),paste("padding-right: 10px;text-align: ", contentAlign,";"))
    firstcell.style <- ifelse(i==nrow(df),paste("padding-right: 10px;border-bottom:solid 1px #000000;text-align: ", leftColAlign,";"),paste("padding-right: 10px;text-align: ", leftColAlign,";"))
    
    
    cell.content <- paste(df[i,],collapse = paste("</td><td style='", cell.style,"'>",sep=""),sep="")
    tbl.summary <- paste(tbl.summary,"<tr><td style='", firstcell.style,"'>", as.character(cell.content),"</td></tr>",sep="")
  }
  
  #Add general not if not Null
  if(!is.null(genNote)){
    tbl.summary <- paste(tbl.summary,"<tr><td colspan='", ncol(df) ,"'><i>Note.</i>&nbsp;", genNote ,"</td></tr>",sep="")
  }
  
  #Add closing table tag
  tbl.summary <- paste(tbl.summary,"</table><br />",sep="")
  
  #Print Final table
  #cat(tbl.summary)
  asis_output(tbl.summary)
  
}