
#roxygen description
#' @name  inv_ecdf
#' @export inv_ecdf
#' @title  inv_ecdf

inv_ecdf <- function(f){
  x <- environment(f)$x
  y <- environment(f)$y
  approxfun(y, x)
}
