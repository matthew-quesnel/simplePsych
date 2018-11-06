#' Linear Regression and Related Functions
#'
#' This function allows you to print the output of data frames and most R functions in an APA formatted table.
#' @param ivs List of variables from which to compute all pairwise correlations. Input must be a list (i.e., c("var1","var2")). Defaults to NULL.
#' @param data A data frame containing the observed variables to correlate. Defaults to NULL.
#' @param method Type of correlation coefficient to be computed. Options include "pearson","kendall", or "spearman". Defaults to "pearson".
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
#' @import lm.beta
#' @import boot
#' @import plyr
#' @import dplyr
#' @import car
#' @import ggplot2
#' @import jmvcore
#' @import jtools
#' @examples
#' cormatrix(ivs = c("mpg","wt","hp"), data=mtcars, digits=3, output="html", columnnames = FALSE, CI=TRUE, CIlevel=.95, displayDF=TRUE, displayP=TRUE, showSig=TRUE)


regression <- function(y=NULL,
                       model=NULL,
                       simple.slopes=NULL,
                       conf=.95,
                       digits=3,
                       output="r",
                       assumptions=TRUE,
                       data=NULL,
                       boot=NULL,
                       boottype="perc",
                       bootci=0.95,
                       plot=FALSE){

  if (is.null(data)) stop('No dataframe provided')
  #if (is.null(y)) stop('No dependent variable provided')
  if (is.null(model)) stop('No model provided')

  lm.model.object2 <- list()
  modelnum <- 0
  modelname <- list()
  for (formula.item in model){
    modelnum <- modelnum + 1
    modelname <- c(modelname,paste("Model",modelnum))
    lm.model <- as.formula(formula.item)
    lm.obj <- lm(lm.model,data = data)



     vars <- all.vars(lm.model)
     #data[,vars] <- lapply(data[,vars], as.numeric)
     #data[vars] <- jmvcore::toNumeric(data[vars])
     if (any(lapply(data[,vars],is.factor)==TRUE)) stop('No factors allowed')
    # data[vars] <- as.numeric(data[vars])

    lm.model.object2 <- c(lm.model.object2,list(lm.obj))
  }

  if(length(model)>1){
    modelcompare <- do.call(anova, lm.model.object2)
    class(modelcompare)
    rownames(modelcompare) <- modelname
    #modelcompare <- broom::tidy(modelcompare)
    modelcompare <- as.data.frame(modelcompare)
    modelcompare <- plyr::rename(modelcompare, c("Res.Df"= "Residual Df","RSS"="Residual Sum Sq","Pr(>F)"="p"))

    ##Round p-values and remove "NAs"
    modelcompare <- round_df(modelcompare,digits)
    modelcompare$p <- format.pval(modelcompare$p, eps = .001, digits = digits)
    modelcompare[is.na(modelcompare)] <- ""
    modelcompare[modelcompare=="NA"] <- ""

    if (output=="r"){
      cat("Model Comparison \n", sep = "\n")
      print(modelcompare)
    }else if (output=="html"){
      cat(APAhtmlTable(modelcompare, round = FALSE, tidy=FALSE, digits=digits, tableTitle="Model Comparison"))
    }
  }

  modelnum <- 0
  for (formula.item in model){
    modelnum <- modelnum + 1
    ##Linear Model
    lm.model <- as.formula(formula.item)
    lm.model.object <- lm(lm.model,data = data)
    conf.int <- confint(lm.model.object,level = conf)
    lm.summary <- summary(lm.model.object)
    totSummary.tbl <- tidy(lm.summary)


    #Standardized Beta
    lm.D9.beta <- lm.beta(lm.model.object)
    beta <- coef(lm.D9.beta, standardized = TRUE)

    #Partial Correlations

    #create vars
    vars <- all.vars(formula(lm.model.object))
    if (is.factor(data[,vars])) stop('No factors allowed')


    #Combine and Round everything
    colnames(conf.int) <- c("LLCI","ULCI")
    totSummary.tbl <- cbind(totSummary.tbl,beta,conf.int)
    totSummary.tbl <- rename(totSummary.tbl, c("std.error"="SE","statistic"="t","p.value"="p","beta"="β"))

    if (nrow(totSummary.tbl)>2){
      partcor <- vector()
      partcor[1] <- NA
      varitems <- list()
      d <- matrix(nrow=nrow(data))
      cname <- list()
      int.term <- FALSE
      for(i in 1:nrow(totSummary.tbl)){
        if(i<nrow(totSummary.tbl)){

          w <- i+1


          if(grepl(":", totSummary.tbl$term[w])==TRUE){
            vars3 <- strsplit(totSummary.tbl$term[w],split=":")
            int.term <- TRUE

            #multiply all interaction terms
            int1 <- ""
            int1 <- paste(vars3[[1]],collapse=".")
            int <- 0
            for(i in 1:length(vars3[[1]])){
              if(i==1){
                int <- data[vars3[[1]][i]]
              }else{
                int <- int*data[vars3[[1]][i]]
              }
            }

            additem <- int
            cname <- c(cname,int1)
          }else{
            varitems <- totSummary.tbl$term[w]
            #print(data[varitems])
            additem <- data[varitems]
            cname <- c(cname,totSummary.tbl$term[w])
          }
          if(w==2){
            d <- cbind(additem)
          }else{
            d <- cbind(d,additem)
          }

        }
      }
      additem <- as.matrix(data[vars[1]])
      d <- cbind(additem,d)
      cname <- c(vars[1],cname)
      colnames(d) <- cname
      d <- as.data.frame(d)
      for(w in 2:nrow(totSummary.tbl)){
        partcor[w] <- semicorR(x=as.character(cname[w]),y=as.character(cname[1]),cov=as.character(cname[-c(1,w)]), data=d)
      }
      totSummary.tbl <- cbind(totSummary.tbl, partcor)


    }else{
      r <- vector()
      r[1] <- NA
      r[2] <- cor(data[,vars[1]],data[,vars[2]])
      totSummary.tbl <- cbind(totSummary.tbl,r)

    }

    ##Handle Bootstrapped Regressions

    if (!is.null(boot)){

      bootResults <- boot(statistic=bootReg, formula=lm.model, data=data, R=boot)

      ci <- matrix(nrow =nrow(totSummary.tbl), ncol=5)
      boot.se <- vector()
      boot.bias <- vector()
      boot.p <- vector()
      null.dist <- vector()

      for(i in 1:nrow(totSummary.tbl)){
        if (boottype=="perc"){
          ci[i,] <- boot::boot.ci(bootResults, type=boottype, index=i, conf=conf)$percent
        }else if(boottype=="bca"){
          ci[i,] <- boot::boot.ci(bootResults, type=boottype, index=i, conf=conf)$bca
        }
        boot.se[i] <- sqrt(var(bootResults$t[,i]))
        boot.bias[i] <- mean(bootResults$t[,i]) - totSummary.tbl$estimate[i]

        # The bootstrap p-value can then be approximated by
        null.dist <- bootResults$t[,i] - mean(bootResults$t[,i])
        bootP <- (1+sum(abs(null.dist) > abs(bootResults$t0[i])))/(1+bootResults$R)
        boot.p[i] <- bootP

      }
      #Remove extra information
      ci<-ci[,-c(1,2,3)] #delete columns 6, 7 and 8

      #Add Boostrap columns
      totSummary.tbl <- cbind(totSummary.tbl,boot.bias,boot.se,ci)

      # Round dataframe
      totSummary.tbl <- round_df(totSummary.tbl,digits)

      totSummary.tbl$p <- format.pval(totSummary.tbl$p, eps = .001, digits = digits)

      #Rename and Reorder
      if (nrow(totSummary.tbl)>2){
        totSummary.tbl <- totSummary.tbl[c("term","estimate","β","SE","boot.bias","boot.se","t","p","partcor","1","2")]
        totSummary.tbl <- plyr::rename(totSummary.tbl, c("term"="Term","estimate"="Coefficient","boot.bias"="Bias","boot.se"="Boot SE","1"="Boot LLCI","2"="Boot ULCI","partcor"="Semi-partial r"))
        #totSummary.tbl["Semi-partial r"][is.na(totSummary.tbl["Semi-partial r"])] <- ""
      }else{
        totSummary.tbl <- totSummary.tbl[c("term","estimate","β","SE","boot.bias","boot.se","t","p","r","1","2")]
        totSummary.tbl <- plyr::rename(totSummary.tbl, c("boot.bias"="Bias","boot.se"="Boot SE","1"="Boot LLCI","2"="Boot ULCI"))
      }
    }else{
      #Round dataframe
      totSummary.tbl <- round_df(totSummary.tbl,digits)

      totSummary.tbl$p <- format.pval(totSummary.tbl$p, eps = .001, digits = digits)

      if (nrow(totSummary.tbl)>2){
        totSummary.tbl <- totSummary.tbl[c("term","estimate","LLCI","ULCI","β","SE","t","p","partcor")]
        totSummary.tbl <- plyr::rename(totSummary.tbl, c("term"="Term","estimate"="Coefficient","partcor"="Semi-partial r"))
        #totSummary.tbl["Semi-partial r"][is.na(totSummary.tbl["Semi-partial r"])] <- ""
      }else{
        totSummary.tbl <- totSummary.tbl[c("term","estimate","LLCI","ULCI","β","SE","t","p","r")]
        totSummary.tbl <- plyr::rename(totSummary.tbl, c("term"="Term","estimate"="Coefficient"))
      }
    }




    #Calculate simple slopes
    slopemeans <- NULL
    slopetest <- NULL
    #slopesT <- FALSE
    # atval <- list() #vector("list",length(simple.slopes))
    # #atval <- "list("
    # lsformula <- paste("~",paste(simple.slopes,collapse=":"))

   #print(simpleTest(model=lm.model.object,simple.slopes=simple.slopes, data=data, digits=digits))


      # for(var in simple.slopes){
    #
    #     unique_x <- length(unique(data[[var]]))
    #
    #     if(unique_x==2){
    #       #Add parameters to list for at= command
    #       atval[[var]] <- c(0,1)
    #       #jtools::sim_slopes(lm.model.object, pred = simple.slopes[var], modx = var[2], modx.values=unique(data[[var]]), johnson_neyman = FALSE, digits=digits)
    #     }else{
    #       plusSD <- 0
    #       minusSD <- 1
    #       #print(sd(data[[var]]))
    #       #plusSD <- mean(data[,var])+sd(data[[var]])
    #       #minusSD <- mean(data[[,var])-sd(data[[var]])
    #       plusSD <- mean(data[[var]])+sd(data[[var]])
    #       minusSD <- mean(data[[var]])-sd(data[[var]])
    #       slopeval <- list(round(minusSD,digits),round(plusSD,digits))
    #
    #
    #       #Add parameters to list for at= command
    #       atval[[var]] <- c(as.numeric(slopeval[1]),as.numeric(slopeval[2]))
    #       #jtools::sim_slopes(lm.model.object, pred = var[1], modx = mindset, modx.values=c(0,1), johnson_neyman = FALSE, digits=digits)
    #     }
    #
    # }
    #
    # slopesT <- is.contained(simple.slopes, vars)
    # if(!is.null(simple.slopes) && slopesT==TRUE && int.term==TRUE){
    #   #Estimate y at values of var
    #   slopemeans <- lsmeans(lm.model.object,as.formula(lsformula), at=atval)
    #   slopetest <- summary(pairs(slopemeans), adjust="none")
    #   slopemeans <- summary(slopemeans)
    #   slopemeans <- plyr::rename(slopemeans, c("lsmean"="Predicted Mean","lower.CL"="LLCI","upper.CL"="ULCI"))
    #
    #   #Tests of Slopes
    #   slopetest <- round_df(slopetest,digits)
    #   slopetest$p.value <- format.pval(slopetest$p.value, eps = .001, digits = digits)
    #   slopetest <- plyr::rename(slopetest, c("contrast"="Contrast","estimate"="Estimate","t.ratio"="t","p.value"="p"))
    #
    #   #TRY to test slopes with continuous
    #
    #     data[[simple.slopes[2]]] <- as.factor(data[[simple.slopes[2]]])
    #
    #   slopemeans <- lsmeans(lm.model.object,as.formula(lsformula))
    #   slopetest <- summary(pairs(slopemeans), adjust="none")
    #   slopemeans <- summary(slopemeans)
    #   slopemeans <- plyr::rename(slopemeans, c("lsmean"="Predicted Mean","lower.CL"="LLCI","upper.CL"="ULCI"))
    #
    #   #Tests of Slopes
    #   slopetest <- round_df(slopetest,digits)
    #   slopetest$p.value <- format.pval(slopetest$p.value, eps = .001, digits = digits)
    #   slopetest <- plyr::rename(slopetest, c("contrast"="Contrast","estimate"="Estimate","t.ratio"="t","p.value"="p"))
    #   }

    #Print output
    if (output=="r"){
      options(width = 200)
      cat(paste("\n",modelname[modelnum],"\n"), sep = "\n")
      #Print Header
      cat(paste("OLS Regression \n \n Model = ", as.character(formula.item), "\n Outcome = ", vars[1],"\n Predictors = ",paste(vars[2:length(vars)], collapse = ", "),"\n Observations = ",nrow(data),"\n"),sep="")
      if(!is.null(boot)){
        if (boottype=="perc"){
          print.boottype <- "Percentile"
        }else if (boottype=="bca"){
          print.boottype <- "BCa"
        }
        cat(" Bootstrap Type = ", print.boottype,"\n Bootstrap Samples = ", boot,"\n Bootstrap Confidence Interval = ", as.character(bootci),  "\n \n",sep="")
      }
      ## Assumptions ##
      if(assumptions==TRUE){
        cat("\nCheck Assumptions of Model\n", sep="")
        fit <- lm.model.object
        #####---------- Plots of the residuals ------------######
        par(mfrow=c(1,2),mar=c(2, 4, 2, 1) + 0.1)
        plot(lm.model.object)

        model.diag.metrics <- augment(fit)
        leverage<- model.diag.metrics %>%
          filter(.hat>2*(length(fit$terms)+1)/(fit$df.residual+length(fit$terms)+1))

        influence <- model.diag.metrics %>%
          filter(.cooksd>(4/fit$df.residual))

        cook.extreme <- which(model.diag.metrics$.cooksd>(4/fit$df.residual))
        par(mfrow=c(1,1))
        plot(model.diag.metrics$.cooksd, ylab="Cook's D Values")
        abline(h=(4/fit$df.residual), lty=2)

        if(length(influence>0)){
          text(cook.extreme, influence$.cooksd, labels=influence$.rownames, pos=4)
        }

        cat("\nHigh Leverage Observations (hat > 2(k+1)/n)\n")
        print(leverage, row.names = FALSE)

        cat("\nHigh Influence Observations (Cook's D > 4/(n-k-1))\n")
        print(influence, row.names = FALSE)

        #Print VIF
        if (nrow(totSummary.tbl)>2){
          VIFactor <- car::vif(lm.model.object)
          VIFactor <- broom::tidy(VIFactor)
          Tolerance <- 1/VIFactor$x
          VIFactor <- plyr::rename(VIFactor, c("names"="Terms","x"="VIF"))
          VIFactor <- cbind(VIFactor,Tolerance)
          cat("\nVariance Inflation Factors\n")
          print(VIFactor, row.names = FALSE)
        }

        #Print Durban Watson
        cat("\nDurban-Watson Test\n")
        dbw.test <- car::durbinWatsonTest(lm.model.object)
        dbw.test <- as.vector(dbw.test)
        dbw.test2 <- cbind(dbw.test$r,dbw.test$dw,dbw.test$p)
        colnames(dbw.test2) <- c("Autocorrelation","D-W Statistic","p")
        dbw.test2 <- as.data.frame(dbw.test2)
        print(dbw.test2, row.names = FALSE)
      }
      cat("\nRegression results", sep="\n")
      print(totSummary.tbl, row.names = FALSE)
      f.sig = pf(lm.summary$fstatistic[1], lm.summary$fstatistic[2],lm.summary$fstatistic[3], lower.tail = FALSE)
      cat("\nF(",lm.summary$fstatistic[2],",",lm.summary$fstatistic[3],") = ",round(lm.summary$fstatistic[1],digits), ", p = ", round(f.sig,digits) ,", R^2 = ", round(lm.summary$r.squared, digits), ", Adjusted R^2 = ", round(lm.summary$adj.r.squared, digits),"\n", sep="")

      ## Check if simple.slope includes an interaction term and if interaction term is in the model
      slopesT <- list()
      int.term <- list()
      for(simple.slopes2 in simple.slopes){
         ##Check if simple.slope includes an interaction term
         int.term[simple.slopes2] <- grepl("*", simple.slopes2, fixed = TRUE)
          simple.slopes2 <- sort(unlist(strsplit(simple.slopes2, "*", fixed = TRUE)))
          for (vars in totSummary.tbl[["Term"]]){
            RegTerms <- sort(unlist(strsplit(vars, ":", fixed = TRUE)))
            test_slope <- identical(as.list(simple.slopes2),as.list(RegTerms))
            if(test_slope==TRUE){
              break
          }
        }
        if(test_slope==TRUE){
          slopesT <- c(slopesT,TRUE)
        }else{
          slopesT <- c(slopesT,FALSE)
        }
     }
      #Print simple slopes if requested, and if interaction term is in request and in the model
      if(!is.null(simple.slopes) && all(slopesT==TRUE) && all(int.term==TRUE)){
        cat("\n")
        simpleTest(formula=formula.item,simple.slopesV=simple.slopes,data=data,digits=digits, output=output)
      }else if(!all(int.term==TRUE)){
        cat("\nSimple Slopes Analysis\n")
        cat("\nWarning: One or more simple.slopes provided does not contain an interaction term. Must be in form 'x1*x2'\n")
      }else if(!all(slopesT==TRUE)){
        cat("\nSimple Slopes Analysis\n")
        cat("\nWarning: One or more interaction terms in simple.slopes was not found in the model: ",model," \n")
      }
    }else if (output=="html"){
      cat(paste("<h3>", modelname[modelnum],"</h3>"))

      #Print Header
      cat(paste("OLS Regression <br/><br/> Model = ", as.character(formula.item), "<br/> Outcome = ", vars[1],"<br/> Predictors = ",paste(vars[2:length(vars)], collapse = ", "),"<br/> Observations = ",nrow(data),"<br/>"),sep="")
      if(!is.null(boot)){
        if (boottype=="perc"){
          print.boottype <- "Percentile"
        }else if (boottype=="bca"){
          print.boottype <- "BCa"
        }
        cat(" Bootstrap Type = ", print.boottype,"<br/> Bootstrap Samples = ",boot,"<br/> Bootstrap Confidence Interval = ", as.character(bootci),  "<br/><br/>",sep="")
      }

      ## Assumptions ##
      if(assumptions==TRUE){
        cat("<h5>Check Assumptions of Model</h5><br/>", sep="")
        fit <- lm.model.object
        #####---------- Plots of the residuals ------------######
        par(mfrow=c(1,2),mar=c(2, 4, 2, 1) + 0.1)
        plot(lm.model.object)

        model.diag.metrics <- augment(fit)
        leverage<- model.diag.metrics %>%
          filter(.hat>2*(length(fit$terms)+1)/(fit$df.residual+length(fit$terms)+1))

        influence <- model.diag.metrics %>%
          filter(.cooksd>(4/fit$df.residual))

        cook.extreme <- which(model.diag.metrics$.cooksd>(4/fit$df.residual))

        par(mfrow=c(1,1))
        plot(model.diag.metrics$.cooksd, ylab="Cook's D Values")
        abline(h=(4/fit$df.residual), lty=2)

        if(length(influence>0)){
          text(cook.extreme, influence$.cooksd, labels=influence$.rownames, pos=4)
        }

        #Print Leverage
        cat(APAhtmlTable(leverage,round = TRUE, digits=3, tidy=FALSE,tableTitle="High Leverage Observations (hat > 2(k+1)/n)"))

        #Print Leverage
        cat(APAhtmlTable(influence,round = TRUE, digits=3, tidy=FALSE,tableTitle="High Influence Observations (Cook's D > 4/(n-k-1))"))

        if (nrow(totSummary.tbl)>2){
          VIFactor <- car::vif(lm.model.object)
          #Print VIF
          VIFactor <- broom::tidy(VIFactor)
          Tolerance <- 1/VIFactor$x
          VIFactor <- cbind(VIFactor,Tolerance)
          cat(APAhtmlTable(VIFactor,round = TRUE, digits=3, tidy=FALSE,tableTitle="Variance Inflation Factors",header = c("Terms","VIF","Tolerance")))
        }

        #Print Durban Watson
        dbw.test <- car::durbinWatsonTest(lm.model.object)
        dbw.test <- as.vector(dbw.test)
        dbw.test2 <- cbind(dbw.test$r,dbw.test$dw,dbw.test$p)
        colnames(dbw.test2) <- c("Autocorrelation","D-W Statistic","p")
        cat(APAhtmlTable(dbw.test2,round = TRUE, digits=3, tidy=TRUE,tableTitle="Durban-Watson Test"))
      }
      cat("<br/><h5>Regression Model Output</h5>")
      #Model Parameters
      f.sig = pf(lm.summary$fstatistic[1], lm.summary$fstatistic[2],lm.summary$fstatistic[3], lower.tail = FALSE)
      model.params <- paste("F(",lm.summary$fstatistic[2],",",lm.summary$fstatistic[3],") = ",round(lm.summary$fstatistic[1],digits), ", p = ",round(f.sig,digits) , ", R^2 = ",round(lm.summary$r.squared,digits), ", Adjusted R^2 = ",round(lm.summary$adj.r.squared,digits),sep="")
      #Print APA style Table
      cat(APAhtmlTable(totSummary.tbl,round = FALSE, tidy=FALSE,tableTitle="Regression results",genNote=model.params))

      ##Check if simple.slope includes an interaction term
      int.term <- ifelse((grepl("*", simple.slopes)==TRUE),TRUE,FALSE)
      ## Check if interaction term is in the model
      simple.slopes2 <- sort(unlist(strsplit(simple.slopes, "*", fixed = TRUE)))
      for (vars in totSummary.tbl[["Term"]]){
        RegTerms <- sort(unlist(strsplit(vars, ":", fixed = TRUE)))
        slopesT <- identical(as.list(simple.slopes2),as.list(RegTerms))
        if(slopesT==TRUE){
          break
        }
      }

      #Print simple slopes if requested
      if(!is.null(simple.slopes) && slopesT==TRUE && int.term==TRUE){
        cat("<h5>Simple Slopes Analyses</h5>")
        cat(simpleTest(formula=formula.item,simple.slopesV=simple.slopes,data=data,digits=digits, output=output))
      }
    }
    if(plot==TRUE){
      if(length(vars)>2 && slopesT==TRUE && int.term==TRUE){
        slopemeans <- lsmeans(lm.model.object,as.formula(lsformula), at=atval)
        slopemeans <- summary(slopemeans)
        slopemeans[,2] <- as.factor(slopemeans[,2])
        if(length(simple.slopes)==2){
          means.lineplot <- ggplot(data=slopemeans, aes_string(x=as.character(simple.slopes[1]), y=as.character("lsmean"), group=as.character(simple.slopes[2]))) +
          #means.lineplot <- ggplot(data=slopemeans, aes(x=eval(parse(text=simple.slopes[1])), y=lsmean, group=eval(parse(text=simple.slopes[2])))) +
            geom_line(aes_string(linetype=simple.slopes[2])) +
            geom_point() +
            theme_bw() +
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))
          print(means.lineplot)
        }else if(length(simple.slopes)>2){
          wrap <- paste("~ ", simple.slopes[3:length(simple.slopes)], collapse="+")
          means.lineplot <- ggplot(data=slopemeans, aes_string(x=as.character(simple.slopes[1]), y=as.character("lsmean"), group=as.character(simple.slopes[2]))) +
          facet_wrap(as.formula(wrap),strip.position = "bottom") +
          geom_line(aes_string(linetype=simple.slopes[2])) +
          geom_point() +
          theme_bw() +
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))
          print(means.lineplot)
        }
      }else if(length(vars)>2 && int.term==FALSE){


            # atval <- c(min(data[,varnum]),max(data[,varnum]))
            # slopemeans <- lsmeans(lm.model.object, at=atval)
            # slopemeans <- summary(slopemeans)
            # cat(APAhtmlTable(slopemeans, tidy=FALSE))

           for(varnum in vars){
             if(!varnum==vars[1]){
               newvars <- vector()
               new.data <- as.data.frame(matrix(0, nrow = 2))
               new.data <- cbind(c(min(data[,varnum]),max(data[,varnum])))
               newvars <- c(newvars,varnum)
               for (varnum2 in vars){
                 if(!varnum2==vars[1] && !varnum2==varnum){
                   new.data <- cbind(new.data,c(rep(mean(data[,varnum2]),2)))
                   newvars <- c(newvars,varnum2)
                 }
               }

               colnames(new.data) <- newvars
               new.data <- as.data.frame(new.data)
               prediction <- predict(lm.model.object, newdata=new.data)
               new.data <- cbind(new.data,prediction)
               means.lineplot <- ggplot(data=new.data,aes_string(x=as.character(varnum), y=as.character("prediction"))) +
                  geom_line() +
                  geom_point() +
                  theme_bw() +
                  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))
                print(means.lineplot)
             }
           }
      }else if(length(vars)==2){
        slopemeans <- lsmeans(lm.model.object,as.formula(paste("~",vars[2])))
        slopemeans <- summary(slopemeans)
        #slopemeans[,2] <- as.factor(slopemeans[,2])
        #cat(APAhtmlTable(slopemeans,tidy=FALSE))

        means.lineplot <- ggplot(data=data,aes_string(x=as.character(vars[2]), y=as.character(vars[1]))) +
          geom_smooth(method='lm',formula=y~x, se=FALSE, colour="black") +
          theme_bw() +
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))

        # means.lineplot <- ggplot(data=slopemeans, aes_string(x=as.character(vars[2]), y=as.character("lsmean"))) +
        # geom_line()+
        # geom_point()
        print(means.lineplot)
      }
    }
  }
}

VectorIntersect <- function(v,z) {
  unlist(lapply(unique(v[v%in%z]), function(x) rep(x,min(sum(v==x),sum(z==x)))))
}
is.contained <- function(v,z) {length(VectorIntersect(v,z))==length(v)}

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
simpleTest <- function(formula=NULL,simple.slopesV=NULL,data=NULL, digits=3, output="r"){
  for(slopes in simple.slopesV){
    simple.slopes <- unlist(strsplit(slopes,split="*", fixed=TRUE))

    model <- lm(formula=as.formula(formula),data=data)
    for(i in 1:length(simple.slopes)){
      mods <- simple.slopes[-i]
      if(length(simple.slopes)==2){
        #mods <- c(mods,NULL)
        unique_x <- length(unique(data[[mods[1]]]))
        ## Set levels of dichotomous mod_x
        if (unique_x==2){
          mod.level <- unique(data[[mods[1]]])
        }else{
          mod.level <- NULL
        }
        slopes.output <- do.call(jtools::sim_slopes, list(model, pred=as.name(simple.slopes[i]),modx=as.name(mods[1]),modx.values=mod.level, johnson_neyman = FALSE, digits=digits))
        if(output=="r"){
          print(slopes.output)
        }else if(output=="html"){
          title <- paste("Slope of", simple.slopes[i], "at levels of", mods[1])
          slopes.output <- broom::tidy(slopes.output)
          #slopes.output$p.value <- format.pval(slopes.output$p.value, eps = .001, digits = digits)
          slopes.output <- plyr::rename(slopes.output, c("estimate"="Estimate","std.error"="SE","statistic"="t","p.value"="p","conf.low"="LLCI","conf.high"="ULCI","term"="Moderator"))
          slopes.output <- slopes.output[c("Estimate","SE","t","p","LLCI","ULCI","Moderator")]
          cat(simplePsych::APAhtmlTable(slopes.output, tableTitle=title, round=TRUE, digits=digits, tidy=FALSE,leftColAlign = "left"))
        }
      }else if(length(simple.slopes)==3){
        unique_x <- length(unique(data[[mods[1]]]))
        unique_x2 <- length(unique(data[[mods[2]]]))
        ## Set levels of dichotomous mod_x
        if (unique_x==2){
          mod.level <- unique(data[[mods[1]]])
        }else{
          mod.level <- NULL
        }
        ##Mod2 levels
        if (unique_x2==2){
          mod2.level <- unique(data[[mods[2]]])
        }else{
          mod2.level <- NULL
        }
        slopes.output <- do.call(jtools::sim_slopes, list(model, pred=as.name(simple.slopes[i]),modx=as.name(mods[1]),modx.values=mod.level, mod2=as.name(mods[2]), mod2.values=mod2.level, johnson_neyman = FALSE, digits=digits))

        ##Print Output appropriately
        if(output=="r"){
          print(slopes.output)
        }else if(output=="html"){
          title <- paste("Slope of", simple.slopes[i], "at levels of", mods[1], "and", mods[2])
          slopes.output <- broom::tidy(slopes.output)
          #slopes.output$p.value <- format.pval(slopes.output$p.value, eps = .001, digits = digits)
          slopes.output <- plyr::rename(slopes.output,
                                        c("estimate"="Estimate","std.error"="SE","statistic"="t","p.value"="p"
                                          ,"conf.low"="LLCI","conf.high"="ULCI","term"="Moderator1","mod2.term"="Moderator2"))
          slopes.output <- slopes.output[c("Estimate","SE","t","p","LLCI","ULCI","Moderator1","Moderator2")]
          cat(simplePsych::APAhtmlTable(slopes.output, tableTitle=title, round=TRUE, digits=digits, tidy=FALSE,leftColAlign = "left"))
        }

      }
    }
  }
}

PlotSimples <- function(formula=NULL,simple.slopes=NULL,data=NULL, digits=3, output="r"){
  for(slopes in simple.slopesV){
    simple.slopes <- unlist(strsplit(slopes,split="*", fixed=TRUE))

    model <- lm(formula=as.formula(formula),data=data)
    for(i in 1:length(simple.slopes)){
      mods <- simple.slopes[-i]
      if(length(simple.slopes)==2){
        #mods <- c(mods,NULL)
        unique_x <- length(unique(data[[mods[1]]]))
        ## Set levels of dichotomous mod_x
        if (unique_x==2){
          mod.level <- unique(data[[mods[1]]])
        }else{
          mod.level <- NULL
        }
        slopes.output <- do.call(jtools::sim_slopes, list(model, pred=as.name(simple.slopes[i]),modx=as.name(mods[1]),modx.values=mod.level, johnson_neyman = FALSE, digits=digits))
        if(output=="r"){
          print(slopes.output)
        }else if(output=="html"){
          title <- paste("Slope of", simple.slopes[i], "at levels of", mods[1])
          slopes.output <- broom::tidy(slopes.output)
          #slopes.output$p.value <- format.pval(slopes.output$p.value, eps = .001, digits = digits)
          slopes.output <- plyr::rename(slopes.output, c("estimate"="Estimate","std.error"="SE","statistic"="t","p.value"="p","conf.low"="LLCI","conf.high"="ULCI","term"="Moderator"))
          slopes.output <- slopes.output[c("Estimate","SE","t","p","LLCI","ULCI","Moderator")]
          cat(simplePsych::APAhtmlTable(slopes.output, tableTitle=title, round=TRUE, digits=digits, tidy=FALSE,leftColAlign = "left"))
        }
      }else if(length(simple.slopes)==3){
        unique_x <- length(unique(data[[mods[1]]]))
        unique_x2 <- length(unique(data[[mods[2]]]))
        ## Set levels of dichotomous mod_x
        if (unique_x==2){
          mod.level <- unique(data[[mods[1]]])
        }else{
          mod.level <- NULL
        }
        ##Mod2 levels
        if (unique_x2==2){
          mod2.level <- unique(data[[mods[2]]])
        }else{
          mod2.level <- NULL
        }
        slopes.output <- do.call(jtools::sim_slopes, list(model, pred=as.name(simple.slopes[i]),modx=as.name(mods[1]),modx.values=mod.level, mod2=as.name(mods[2]), mod2.values=mod2.level, johnson_neyman = FALSE, digits=digits))

        ##Print Output appropriately
        if(output=="r"){
          print(slopes.output)
        }else if(output=="html"){
          title <- paste("Slope of", simple.slopes[i], "at levels of", mods[1], "and", mods[2])
          slopes.output <- broom::tidy(slopes.output)
          #slopes.output$p.value <- format.pval(slopes.output$p.value, eps = .001, digits = digits)
          slopes.output <- plyr::rename(slopes.output,
                                        c("estimate"="Estimate","std.error"="SE","statistic"="t","p.value"="p"
                                          ,"conf.low"="LLCI","conf.high"="ULCI","term"="Moderator1","mod2.term"="Moderator2"))
          slopes.output <- slopes.output[c("Estimate","SE","t","p","LLCI","ULCI","Moderator1","Moderator2")]
          cat(simplePsych::APAhtmlTable(slopes.output, tableTitle=title, round=TRUE, digits=digits, tidy=FALSE,leftColAlign = "left"))
        }

      }
    }
  }
}
