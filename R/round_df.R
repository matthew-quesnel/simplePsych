#' A Function that rounds all numeric values in a dataframe to the specifed number of decimal places.
#'
#' This function rounds all numeric values in a dataframe to the specifed number of decimal places.
#' @param df Data Frame you want to round. Defaults to NULL.
#' @param digits Number of digits you want to round numbers too. Leave blank if you don't want to round numbers or if the Data Frame contains strings. Defaults to 3.
#' @keywords round
#' @examples
#' APAhtmlTable()
round_df <- function(df=NULL, digits=3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}