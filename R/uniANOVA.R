#' A Function to Conduct Univariate Analysis of (Co)Variance.
#'
#' This function conducts a univariate analysis of (co)variance (ANOVA/ANCOVA). It includes options to print descriptives, tests and plots of model assumptions, estimated marginal means, main and simple effects tests, and APA formatted bar plots.
#' @param dv Dependent variable. Defaults to NULL.
#' @param factors Independent variables (converted to factors). All factors will be multiplied in a full-factorial design. Any factors that you do not wish to be crossed with the others should be entere using the "cov=" command. Defaults to NULL.
#' @param cov Independent variables and/or covariates not to be crossed with other factors. These variables will also not be used for grouping in descriptive statistics. You may still utilize these factors in simple and main effects tests using the "comparisons=" command. Defaults to NULL.
#' @param data A data frame containing the observed variables to correlate. Defaults to NULL.
#' @param descriptives TRUE/FALSE indicating whether to display descriptive statistics for the dependent variable grouped by factors. Defaults to TRUE.
#' @param assumptions TRUE/FALSE indicating whether you would like to print tests of model assumptions and plots of residuals.
#' @param ss A number indicating which type sum of squares to use (i.e., 1,2,3). Defaults to 3. Uses the CAR package to conduct ANOVAs with type 3 sums of squares.
#' @param digits Number of digits you want to round numbers too. Leave blank if you don't want to round numbers or if the Data Frame contains strings. Defaults to NULL.
#' @param emmeans TRUE/FALSE indicating whether to display estimated marginal means grouped by factors.
#' @param conf A number indicating the desired confidence interval. Defaults to .95.
#' @param comparisons A list of comparisons for which you would like to conduct main or simple effects tests. Prints estimated marginal means and tests of effects (including effect size) for each comparison specified in the list. Defaults to NULL.
#' @param postcor A string indicating type of correction to apply to main or simple effects tests requested using the "comparisons=" commond. Options include "none", tukey", "scheffe", "sidak", "bonferroni", "dunnettx", "mvt". Defaults to "none".
#' @param plotBar TRUE/FALSE indicating whether a QQPlot bar plot should be returned for estimated marginal means across all factors, with error bars indicating confidence intervals at levels set by "conf". Defaults to FALSE.
#' @param xLabel A string indicating the label for the x-axis in the bar plot if one is requested. Defaults to NULL, which will display variable name.
#' @param yLabel A string indicating the label for the y-axis in the bar plot if one is requested. Defaults to NULL, which will display variable name.
#' @param fillLabel A string indicating the label for the fill variable in the bar plot if one is requested. Defaults to NULL, which will display variable name.
#' @param yLower A number indicating the lower limit of the y-axis in the bar plot if requested. Defaults to NULL, in which case limits are set by QQPlot.
#' @param yUpper A number indicating the upper limit of the y-axis in the bar plot if requested. Defaults to NULL, in which case limits are set by QQPlot.
#' @param yBreaks A number indicating the break value (i.e., interval of numbers and ticks) for the y-axis in the bar plot if requested. Defaults to NULL, in which case limits are set by QQPlot.
#' @param output A string indicating whether results should be displayed as R output (output="r") or HTML output "output="html". Defaults to r output.
#' @keywords ANOVA
#' @export
#' @import broom
#' @import lsmeans
#' @import psych
#' @import car
#' @import MBESS
#' @import ggplot2
#' @import jmvcore
#' @import plyr
#' @import apaTables
#' @examples
#' uniANOVA(dv = "mpg", factors = c("vs","am"), cov=c("hp"), data = mtcars, descriptives=TRUE, assumptions=TRUE, ss=3, digits=3, emmeans=TRUE, conf=.95, comparisons=c("am","vs","vs*am"), plotBar=TRUE, xLabel="Engine", yLabel="Miles per Gallon (US)", fillLabel="Transmission" yLower=10,yUpper=45,yBreaks=5,output = "html").

uniANOVA <- function(dv=NULL,
                     formula=NULL,
                     factors=NULL,
                     cov=NULL,
                     data=NULL,
                     descriptives=TRUE,
                     assumptions=FALSE,
                     ss=3,
                     digits=3,
                     emmeans=FALSE,
                     conf=.95,
                     comparisons=NULL,
                     postCorr="none",
                     plotBar=FALSE,
                     xLabel=NULL,
                     yLabel=NULL,
                     fillLabel=NULL,
                     yLower=NULL,
                     yUpper=NULL,
                     yBreaks=NULL,
                     output="r"){

  #Check for required values and return error if missing (NULL).
  if (is.null(data)) stop('No dataframe provided')
  if (is.null(dv)) stop('No dependent variable provided')
  if (is.null(factors)) stop('No factors provided: Atleast 1 required')

  data <- data
  #set depdendent variable to as.numeric.
  data[[dv]] <- jmvcore::toNumeric(data[[dv]])

  #set variables listed in factors to as.factor.
  for (varName in factors){
    data[[varName]] <- as.factor(data[[varName]])
  }

  #Create formula from factors and covariates.
  by.formula <- paste(factors, collapse="*")

  if(length(cov)>0){
    cov.formula <- paste(cov, collapse="+")
    formula <- paste(dv,"~",by.formula,"+",cov.formula)
  }else{
    formula <- paste(dv,"~",by.formula)
  }
  pformula <- formula
  formula <- as.formula(formula)

  #Compute and print anova
  if(ss==3){
    options(contrasts=c("contr.sum", "contr.poly"))
    model.aov2 <- aov(formula=formula, data=data)
    drop1(model.aov2,~.,test="F")
    results <- car::Anova(model.aov2, type=3, singular.ok=FALSE, silent=TRUE)
  }else{
    model.aov2 <- aov(formula=formula, data=data)
    results <- anova(aov(formula=formula, data=data))
  }

  #Create and Print ANOVA Table
  table.matrix3 <- broom::tidy(results)
  MeanSquare <- table.matrix3$sumsq/table.matrix3$df
  table.matrix3 <- cbind(table.matrix3,MeanSquare)

  #Computer partial eta squared and bind to table
  eta <- table.matrix3$sumsq / (table.matrix3$sumsq + table.matrix3$sumsq[nrow(table.matrix3)])
  table.matrix3 <- cbind(table.matrix3,eta)

  #Compute 95% CI of eta and bind to table
  eta.limit <- NULL
  n <- sum(table.matrix3$df)
  for(i in 1:(nrow(table.matrix3)-1)){
    lower <- round(apaTables::get.ci.partial.eta.squared(F.value=table.matrix3$statistic[i], df1=table.matrix3$df[i], df2=table.matrix3$df[nrow(table.matrix3)], conf.level=.90)$LL,digits)
    upper <- round(apaTables::get.ci.partial.eta.squared(F.value=table.matrix3$statistic[i], df1=table.matrix3$df[i], df2=table.matrix3$df[nrow(table.matrix3)], conf.level=.90)$UL,digits)
    eta.limit <- c(eta.limit,paste("[",lower,", ",upper,"]",sep=""))
  }
  eta.limit <- c(eta.limit,"")
  table.matrix3 <- cbind(table.matrix3,eta.limit)
  SS_error <- table.matrix3[nrow(table.matrix3),2]

  #Reorder and Rename columns in matrix for output
  table.matrix3 <- table.matrix3[,c("term","sumsq","df","MeanSquare","statistic","p.value","eta","eta.limit")]
  table.matrix3 <- plyr::rename(table.matrix3, c("term"="", "sumsq"="Sum of Squares","df"="df","MeanSquare"="Mean Square","statistic"="F","p.value"="p","eta"="ηp2","eta.limit"="ηp2 90% CI"))

  row.names(table.matrix3) <- NULL
  table.matrix3 <- round_df(table.matrix3,digits)
  table.matrix3[nrow(table.matrix3),7] <- ""
  table.matrix3[1:(nrow(table.matrix3)-1), 6] <- format.pval(table.matrix3[1:(nrow(table.matrix3)-1), 6], eps = .001, digits = 3)



  #Compute model Statistics and add to variable for output
  model.lm <- lm(as.formula(formula), data=data)
  sum.model.lm <- summary(model.lm)
  eta <- ef_eta((sqrt(sum.model.lm$fstatistic[1])),df=sum.model.lm$fstatistic[2],dferror=sum.model.lm$fstatistic[3])
  f.sig = pf(sum.model.lm$fstatistic[1], sum.model.lm$fstatistic[2],sum.model.lm$fstatistic[3], lower.tail = FALSE)
  model.params <- paste("*F*(",sum.model.lm$fstatistic[2],",",sum.model.lm$fstatistic[3],") = ",round(sum.model.lm$fstatistic[1],3), ", *p*= ",round(f.sig,3) , ", *R*^2^: ",round(sum.model.lm$r.squared,3), sep="")


  #####---------- Descriptives ----------------######
  if (descriptives==TRUE){
    desc2 <- psych::describeBy(data[[dv]],
                               group = as.list(data[factors]),
                               digits= digits, mat=TRUE)

    desc2 <- desc2[,-which(names(desc2) %in% c("vars","item"))]

    desc2 <- plyr::rename(desc2, c("trimmed"="trim"))
  }


  #####---------- CHECK ASSUMPTIONS OF ANOVA --------------######
  if(assumptions==TRUE){

    #####---------- Levene's test ------------######
    leveneFormula <- paste(dv,"~",by.formula)
    leveneOutput <- car::leveneTest(as.formula(leveneFormula), data=data)
    if (output=="html"){
      leveneOutput <- rename(leveneOutput, c("Df"="df", "F value"="*F*","Pr(>F)"="*p*"))
    }else if (output=="r"){
      leveneOutput <- rename(leveneOutput, c("Df"="df", "F value"="F","Pr(>F)"="p"))
    }
    #Print levene's output
    leveneOutput <- round_df(leveneOutput,digits)


    #####---------- Shapiro-Wilk test ------------######
    shapiro.model <-aov(as.formula(leveneFormula), data=data)
    res <- shapiro.model$residuals
    shapiro.output <- broom::tidy(shapiro.test(res))
    if (output=="html"){
      shapiro.output <- plyr::rename(shapiro.output, c("statistic"="*W*", "p.value"="*p*"))
    }else if (output=="r"){
      shapiro.output <- plyr::rename(shapiro.output, c("statistic"="W", "p.value"="p"))
    }
    shapiro.output <- round_df(shapiro.output[1:2],digits)

  }

  #####--------Save Factors to variable for use in emmeans---------#####
  by4 <- factors

  ####Create GGPlot2 Bar Plot###
  if (plotBar==TRUE){

    if(length(by4)>1){
      factors1 <- as.character(by4[1])
      factors2 <- paste(by4[2:length(by4)], collapse=":")
      factors.formula <- paste("~",factors1,"|",factors2,sep="")
    }else{
      factors1 <- as.character(by4[1])
      factors.formula <- paste("~",factors1,sep="")
    }
    org.lsm4 <- lsmeans::lsmeans(model.aov2, as.formula(factors.formula))

    means<- summary(org.lsm4)

    if(is.null(xLabel)){
      xLabel <- as.character(factors[1])
    }
    if(is.null(yLabel)){
      yLabel <- as.character(dv)
    }
    if(is.null(fillLabel)){
      fillLabel <- as.character(factors[2])
    }

    xlab <- c(attributes(eval(substitute(by4[1]),data))$levels)
    x2lab <- c(attributes(eval(substitute(by4[2]),data))$levels)
    colour <- "Grey30"
    for(i in 2:length(x2lab)){
      grey <- paste("Grey",i*25+15,sep="")
      colour <- c(colour,grey)
    }

    scaleBreaks <- yLower
    for(i in seq(yLower+yBreaks,yUpper,yBreaks)){
      scaleBreaks <- c(scaleBreaks,i)
    }


    if(length(by4)==2){
      f1 <- as.character(by4[2])
      f2 <- as.character(by4[1])
      wrap <- paste("~ ", by4[3:length(by4)], collapse="+")
      means.barplot <- ggplot2::ggplot(means, aes_string(x=as.character(by4[1]), y="lsmean", fill=as.character(by4[2]))) +
        geom_bar(stat="identity", colour="black", position="dodge", width=.75) +
        ylab(yLabel) +
        xlab(xLabel) +
        scale_x_discrete(limits=xlab) +
        scale_y_continuous(breaks=scaleBreaks) +
        scale_fill_manual(name=fillLabel, values = colour) +
        coord_cartesian(ylim = c(yLower, yUpper)) +
        geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.3, position=position_dodge(.75)) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))
    }else if (length(by4)>2){
      wrap <- paste("~ ", by4[3:length(by4)], collapse="+")
      means.barplot <- ggplot2::ggplot(means, aes_string(x=as.character(by4[1]), y="lsmean", fill=as.character(by4[2]))) +
        facet_wrap(as.formula(wrap),strip.position = "bottom") +
        geom_bar(stat="identity", colour="black", position="dodge", width=.75) +
        ylab(yLabel) +
        xlab(xLabel) +
        scale_x_discrete(limits=xlab) +
        scale_y_continuous(breaks=scaleBreaks) +
        scale_fill_manual(name=fillLabel, values = colour) +
        coord_cartesian(ylim = c(yLower, yUpper)) +
        geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.3, position=position_dodge(.75)) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))
    }else if (length(by4)==1){
      f1 <- as.character(by4[1])
      means.barplot <- ggplot2::ggplot(means, aes_string(x=as.character(by4[1]), y="lsmean")) +
        geom_bar(stat="identity", colour="black", position="dodge", width=.75) +
        ylab(yLabel) +
        xlab(xLabel) +
        scale_x_discrete(limits=xlab) +
        scale_y_continuous(breaks=scaleBreaks) +
        coord_cartesian(ylim = c(yLower, yUpper)) +
        geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.3, position=position_dodge(.75)) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(),axis.text = element_text(family="sans", size = 12, colour="black"),axis.title = element_text(family="sans", size = 12, colour="black", face="bold"),legend.text = element_text(family="sans", size = 12, colour="black"), legend.title = element_text(family="sans", size = 12, colour="black", face="bold"),strip.background=element_blank(),strip.placement = "outside",strip.text=element_text(family="sans", size = 12, colour="black", face="bold"))
    }
  }

  if(output=="r"){
    options(width = 200)
    ## Print Title of Analysis
    cat(paste("Univariate Analysis of Variance \n\n Model: ",as.character(pformula),"\n",sep=""))

    ## Print Descriptives and header
    if (descriptives==TRUE){
      title <- paste("\n Descriptives for",as.character(dv),"\n")
      cat(title)
      print(desc2, row.names = FALSE)
    }

    if(assumptions==TRUE){
      cat("\nCheck Model Assumptions\n")
      #Print levene's output
      print(leveneOutput, row.names = FALSE)
      #Print Shapiro-Wilk output
      cat("\nShapiro-Wilk test of normality of residuals \n")
      print(shapiro.output, row.names = FALSE)
      #####---------- Plots of the residuals ------------######
      # cat("\nHistogram of Residuals \n")
      # print(qplot(res, geom = "histogram",
      #             colour = I("black"), fill = I("white"),
      #             xlab = "Residuals", ylab = "Count", binwidth=.5) + theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank()))
      #cat("\nPlots of Residuals \n")
      #par(mfrow=c(2,2))
      #print(plot(model.aov2))
      par(mfrow=c(1,2),mar=c(2, 4, 2, 1) + 0.1)
      plot(model.aov2)
    }

    ## Print Tests of Between-subjects Effects Anova Table
    cat("\n Tests of Between-subjects Effects \n")
    table.matrix3[is.na(table.matrix3)] <- ""
    print(table.matrix3, row.names = FALSE)
    #cat(APAhtmlTable(table.matrix3, tidy=FALSE, digits=digits,  genNote = model.params, tableTitle=caption1, round=FALSE, headerAlign = "center", contentAlign = "right"))

    if (emmeans==TRUE){
      cat("\n Estimated Marginal Means \n")
      emmeans <- suppressWarnings(emMeansOutput(by4=by4,model.aov2=model.aov2,output=output, data=data, digits=digits, conf=conf)) %>% suppressWarnings()
      #print(emmeans, row.names = FALSE)
    }

    ## Print Estimated Margin Means and Simple Effects tests if requested
    if(!is.null(comparisons)){
      for(by6 in comparisons){
        cat(paste("\nEstimated Marginal Means and Pairwise Comparisons:", by6, "\n"))
        model.formula3 <- as.character(formula(model.aov2))
        if(!grepl("*", by6, fixed=TRUE) && any(grepl("*", model.formula3, fixed=TRUE))){
          cat("** NOTE: Results may be misleading due to involvement in interactions \n")
        }
       posthoc.comparisons <- suppressWarnings(postHocComparisons(compfactors=by6,model.aov2=model.aov2,output=output, postCorr=postCorr, data=data, SS_error=SS_error, n=n, digits=digits, conf=conf)) %>% suppressWarnings()
        #print(posthoc.comparisons, row.names = FALSE)
      }
    }

    ## Print bar plot if requested
    if(plotBar==TRUE){
      #cat("\n Bar Graph \n")
      print(means.barplot)
    }

  }else if(output=="html"){
    ## Print Title of Analysis
    cat(paste("<h1>Univariate Analysis of Variance</h1> <br /><h5>Model: ",gsub("\\*", "&#42;", as.character(pformula)),"</h5>",sep=""))
    caption1 <- paste("ANOVA Table: ",  gsub("\\*", "&#42;", pformula))

    ## Print Descriptives and header
    if (descriptives==TRUE){
      title <- paste("Descriptives for",as.character(dv))
      cat(APAhtmlTable(desc2, digits=digits, tableTitle=title, round=TRUE, tidy=FALSE))
    }

    if(assumptions==TRUE){
      cat("<h3>Check Model Assumptions</h3>", sep="")
      #Print levene's output
      cat(APAhtmlTable(leveneOutput, tableTitle="Levene's Test for homogeneity of variance", round=FALSE, tidy=FALSE,leftColAlign = "center"))
      #Print Shapiro-Wilk output
      cat(APAhtmlTable(shapiro.output, tableTitle="Shapiro-Wilk test of normality of residuals", round=FALSE, tidy=FALSE,leftColAlign = "center"))
      #####---------- Plots of the residuals ------------######
      cat("Histogram of Residuals <br />", sep="")
      # print(qplot(res, geom = "histogram",
      #             colour = I("black"), fill = I("white"),
      #             xlab = "Residuals", ylab = "Count", binwidth=.5) + theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank()))
      #
      #par(mfrow=c(2,2))
      #plot(model.aov2)
      par(mfrow=c(1,2),mar=c(2, 4, 2, 1) + 0.1)
      plot(model.aov2)
    }

    ## Print Tests of Between-subjects Effects Anova Table
    cat("<h3>Tests of Between-subjects Effects</h3>")
    cat(APAhtmlTable(table.matrix3, tidy=FALSE, digits=digits,  genNote = model.params, tableTitle=caption1, round=TRUE, headerAlign = "center", contentAlign = "right"))

    if (emmeans==TRUE){
      cat("<h3>Estimated Marginal Means</h3>")
      emmeans <- emMeansOutput(by4=by4,model.aov2=model.aov2,output=output, data=data, digits=digits, conf=conf)
      cat(emmeans)
    }

    ## Print Estimated Margin Means and Simple Effects tests if requested
    if(!is.null(comparisons)){
      for(by6 in comparisons){
        cat(paste("<h3>Estimated Marginal Means and Pairwise Comparisons:", by6, "</h3>"))
        model.formula3 <- as.character(formula(model.aov2))
        if(!grepl("*", by6, fixed=TRUE) && any(grepl("*", model.formula3, fixed=TRUE))){
          cat("** NOTE: Results may be misleading due to involvement in interactions <br/>")
        }
        posthoc.comparisons <- postHocComparisons(compfactors=by6,model.aov2=model.aov2,output=output, postCorr=postCorr, data=data, SS_error=SS_error, n=n, digits=digits, conf=conf)
        #cat(posthoc.comparisons)
      }
    }


    ## Print bar plot if requested
    if(plotBar==TRUE){
      cat("<h3>Bar Graph</h3><br />")
      print(means.barplot)
    }
  }
}


### Calculate Eta-squared ###
ef_eta <- function(t=NULL,f=NULL,sumsq_error=NULL,dfeffect=1,dferror=NULL){
  if (is.null(f) && is.null(t)){
    return("err")
  }else if (is.null(f)){
    f <- t^2
  }
  sumsq_effect <- (f*(sumsq_error/dferror))/dfeffect
  denom <- sumsq_error + sumsq_effect
  return(round((sumsq_effect/denom),digits=3))
}


#Ouput Post Hoc Comparisons
postHocComparisons <- function(compfactors=NULL,model.aov2=NULL,output="html", postCorr="none", data=data, SS_error=NULL, n=NULL, digits=3, conf=.95){
  by4 <- as.list(strsplit(compfactors, "[*]")[[1]])
  by5 <- by4

  termlabels <- labels(terms(model.aov2))
  termlabels <- gsub("\\:","*",termlabels)

  for (by6 in compfactors){
    warn <- paste("comparison term", by6, "not in formula")
    if(!by6 %in% termlabels) stop(warn)
  }

    for(i in 1:length(by4)){
      if(length(by4)>1){
        factors1 <- as.character(by5[1])
        factors2 <- paste(by5[2:length(by5)], collapse="+")
        factors.formula <- paste("~",factors1,"|",factors2,sep="")
        w.t2 <- lsmeans::lsmeans(model.aov2, as.formula(factors.formula))
        levelsof <- paste(as.character(by5[2:length(by5)]), collapse = " and ")
        caption_text <- paste("Pairwise Comparisons of", as.character(by5[1]), "at levels of", levelsof, sep=" ")
      }else{
        factors1 <- as.character(by5[1])
        factors.formula <- paste("~",factors1,sep="")
        w.t2 <- suppressMessages(lsmeans::lsmeans(model.aov2, as.formula(factors.formula)))
        caption_text <- paste("Pairwise Comparisons of", as.character(by5[1]), sep=" ")
      }

      pairwise3 <- suppressMessages(summary(pairs(w.t2), adjust=postCorr)) %>% suppressWarnings()

      #Calculate partial Eta squared for simple effects
      eta <- ef_eta(t=pairwise3$t.ratio,sumsq_error=SS_error,dfeffect=1,dferror=pairwise3$df)
      pairwise3 <- cbind(pairwise3,eta)

      ## Calculate 95% CI of partial eta squared ##
      eta.limit2 <- NULL
      for(y in 1:(nrow(pairwise3))){
        lower<-round(MBESS::ci.pvaf(F.value=(pairwise3$t.ratio[y]^2), df.1=1, df.2=pairwise3$df[y], N=n, conf.level=.90)$Lower.Limit.Proportion.of.Variance.Accounted.for,3)
        upper<-round(MBESS::ci.pvaf(F.value=(pairwise3$t.ratio[y]^2), df.1=1, df.2=pairwise3$df[y], N=n, conf.level=.90)$Upper.Limit.Proportion.of.Variance.Accounted.for,3)
        eta.limit2 <- c(eta.limit2,paste("[",lower,", ",upper,"]",sep=""))
      }
      eta.limit2 <- c(eta.limit2,"")
      pairwise3 <- cbind(pairwise3,eta.limit2[1:nrow(pairwise3)])


      pairwise3[1:(nrow(pairwise3)), "p.value"] <- format.pval(pairwise3[1:(nrow(pairwise3)), "p.value"], eps = .001, digits = 3)

      by6 <- by5
      by5[2:length(by4)] <- by6[1:(length(by4)-1)]
      by5[1] <- by6[length(by4)]

      if (output=="r" && i == 1){
        emmeans <- suppressMessages(lsmeans::lsmeans(model.aov2, as.formula(factors.formula))) #%>% suppressWarnings()
        emmeans <- suppressMessages(summary(emmeans, level=conf)) #%>% suppressWarnings()

        caption1 <- paste("\nEstimated Margin Means:", compfactors,"\n")
        emmeans <-  plyr::rename(emmeans, c("lsmean"="Estimated Mean","lower.CL"="CI Lower","upper.CL"="CI Upper"),warn_missing=FALSE, warn_duplicated=FALSE)
        emmeans <- round_df(emmeans,digits)

        cat(caption1)
        print(emmeans,row.names = FALSE)
      }else if(output=="html" && i == 1){
        emmeans <- suppressMessages(lsmeans::lsmeans(model.aov2, as.formula(factors.formula))) #%>% suppressWarnings()
        avg <-  emmeans@misc$avgd.over

        emmeans <- suppressMessages(summary(emmeans, level=conf)) #%>% suppressWarnings()
        lsmeans.confint <- conf
        caption1 <- paste("Estimated Margin Means:", gsub("\\*", "&#42;", compfactors))
        emmeans <- plyr::rename(emmeans, c("lsmean"="Estimated Mean","lower.CL"="CI Lower","upper.CL"="CI Upper"))
        emmeans <- round_df(emmeans,digits)

        note <-""
        if (!identical(avg, character(0))){
          note <- paste("Results averaged over levels of:", paste(as.character(avg),collapse = ", "),"<br/>")
        }
        note <- paste(note,"Confidence level:", as.character(lsmeans.confint))

        #print(emmeans)
        cat(APAhtmlTable(emmeans, tableTitle=caption1, round=FALSE, digits=digits, tidy=FALSE,leftColAlign = "left", genNote=note)[1])
      }

      #Print Simple Effects Test Table
      if (output=="r"){
        pairwise3 <-  suppressWarnings(plyr::rename(pairwise3, c("contrast"="Contrast","estimate"="Estimate","SE"="SE","df"="df","t.ratio"="t","p.value"="p","eta"="ηp2","eta.limit2[1:nrow(pairwise3)]"="ηp2 90% CI")))
        pairwise3 <- round_df(pairwise3,digits)

        cat("\n",caption_text,"\n")
        print(pairwise3, row.names = FALSE)
      }else if(output=="html"){
        pairwise3 <- plyr::rename(pairwise3, c("contrast"="Contrast","estimate"="Estimate","SE"="SE","df"="df","t.ratio"="*t*","p.value"="*p*","eta"="\u03B7~p~^2^","eta.limit2[1:nrow(pairwise3)]"="\u03B7~p~^2^ 90% CI"))
        cat(APAhtmlTable(pairwise3, tableTitle=caption_text, round=TRUE, digits=digits, tidy=FALSE,leftColAlign = "left"))
      }
    }
}

#Ouput Post Hoc Comparisons
emMeansOutput <- function(by4=NULL,model.aov2=NULL,output="html", data=data, digits=3, conf=.95){
  #by4 <- by5
  by5 <- by4

    if(length(by4)>1){
      factors1 <- as.character(by5[1])
      factors2 <- paste(by5[2:length(by5)], collapse="*")
      factors.formula <- paste("~",factors1,"|",factors2,sep="")
    }else{
      factors1 <- as.character(by5[1])
      factors.formula <- paste("~",factors1,sep="")
    }

    if (output=="r"){
      emmeans <- suppressWarnings(lsmeans::lsmeans(model.aov2, as.formula(factors.formula))) #%>% suppressMessages()
      emmeans <- suppressMessages(summary(emmeans, level=conf))  #%>% suppressWarnings()
      caption1 <- paste("Estimated Margin Means: ", paste(by4, collapse="*"))
      emmeans <- suppressWarnings(plyr::rename(emmeans, c("lsmean"="Estimated Mean","lower.CL"="CI Lower","upper.CL"="CI Upper")))
      emmeans <- round_df(emmeans,digits)

      cat("\n",caption1,"\n")
      print(emmeans, row.names=FALSE)
    }else if(output=="html"){
      emmeans <- lsmeans::lsmeans(model.aov2, as.formula(factors.formula))
      avg <-  emmeans@misc$avgd.over
      emmeans <- suppressWarnings(summary(emmeans, level=conf))  #%>% suppressWarnings()
      caption1 <- paste("Estimated Margin Means: ", gsub("\\*", "&#42;", paste(by4, collapse="*")))
      emmeans <- plyr::rename(emmeans, c("lsmean"="Estimated Mean","lower.CL"="CI Lower","upper.CL"="CI Upper"))
      emmeans <- round_df(emmeans,digits)

      note <-""
      if (!identical(avg, character(0))){
        note <- paste("Results averaged over levels of:", as.character(avg),"<br/>")
      }
      note <- paste(note,"Confidence level:", as.character(conf))
      cat(simplePsych::APAhtmlTable(emmeans, tableTitle=caption1, round=TRUE, digits=digits, tidy=FALSE,leftColAlign = "left"))
    }
}
