#roxygen description
#' @name  MinMaxNormalise
#' @export MinMaxNormalise
#' @param x dataframe or matrix to be scaled (all data must be numeric)
#' @param margin can be 1 or 2. x will be scaled by row if .margin =1 or  column .margin=2
#' @title  Normalise data by min max scaler
#' @usage MinMaxNormalise(df = DRUMLRdata:::DRUML_phos, .margin=2)


MinMaxNormalise <- function(x, .margin=2){
  require(dplyr)

  .min_max <- function(x){
    mn <-min(x, na.rm = T)
    mx <-max(x, na.rm = T)
    out <- (x-mn)/(mx-mn)

    return(out)
  }

  out <- base::apply(x, MARGIN = .margin, FUN = .min_max)

  if(.margin==1){
    out <- out %>% t() %>% data.frame(stringsAsFactors = F)
  } else if (.margin==2) {
    out <- data.frame(out, stringsAsFactors = F)
  }

  return(out)
}
