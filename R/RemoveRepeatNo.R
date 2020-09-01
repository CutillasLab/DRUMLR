# roxygen description
#' @title Update Internal Databases for KSEA
#' @name RemoveRepeatNo
#' @usage  RemoveRepeatNo(colnames(df_input))
#' @export RemoveRepeatNo
#' @description Removes the repeat numbers for triplicate data so that datasets can be filtered by cell name
#' @param repeatseperator The string which separates cell name from repeat number. For DRUML data "__" is used as default.


RemoveRepeatNo <- function(x, repeatseperator = "__") {
  if ("dplyr" %in% (.packages()) == FALSE) {
    library(dplyr)
  }
  out <-
    lapply(
      x,
      FUN = function(x) {
        strsplit(x, repeatseperator)[[1]][[1]]
      }
    ) %>% unlist()
  return(out)
}
