#' A Function to Create APA Formatted HTML Tables
#'
#' This function allows you to print the output of data frames and most R functions in an APA formatted table.
#' @param df Data Frame or Function output to print in an APA formatted HTML table. Defaults to NULL.
#' @param digits Number of digits you want to round numbers too. Leave blank if you don't want to round numbers or if the Data Frame contains strings. Defaults to NULL.
#' @param tidy Do you want to convert the object into a data frame (TRUE/FALSE). This is accomplished using the Tidy() function in the package Broom. Defaults to TRUE.
#' @param genNote General note to be displayed below the table (Automatically begins with NOTE. "Your note here"). Defaults to NULL.
#' @param font Sets the font that text is displayed in. Defaults to Times New Roman.
#' @param fontsize Sets the font size of all text in the table. Defaults to 12pt.
#' @param tableNum Adds a table number above the table. Format is "Table XX". If NULL, no table number is displayed. Defualts to NULL.
#' @param tableTitle Adds a title above the table. If NULL, no title will be displayed. Defaults to NULL.
#' @param headerBold Do you want the column headers to be bolded (TRUE/FALSE)? Defaults to FALSE.
#' @param headerAlign Sets the horizontal alignment of the column headers ("left","center","right"). Defaults to "center".
#' @param leftColAlign Sets the horizontal alignment of the first column in the table ("left","center","right"). Defaults to "left".
#' @param contentAlign Sets the horizontal alignment of all other cells in the table besides the header (see headerAlign) and left column (see leftColAlign).
#' @param header Manually define the column headers as a list. The length of the list must match the number of columns in the table. Defaults to NULL. 
#' @keywords APA Tables
#' @export
#' @import broom
#' @import knitr
#' @examples
#' lm.obj <- lm(mpg~am+disp, data=mtcars)
#' lm.summary <- summary(lm.obj)
#' APAhtmlTable(df=lm.summary, digits=3, tableNum=1, tableTitle="Miles per Gallon Regressed on Transmission Type and Displacement", header=c("Term","Estimate","SE","t","p"))


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