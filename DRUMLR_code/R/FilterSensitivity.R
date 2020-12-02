# Emperical Markers of Sensitivity and resistantce
#' @name FilterSensitivity
#' @export FilterSensitivity
#' @title Filter Drugs by sensitivity data
#' @param sensitivity_database sensitvitiy data with cell lines as rownamed and drugs as columns
#' @param iqr_tolerance IQR range cut off to filter drug list
#' @param n_coverage n_coverage cut off required to ensure models are reliable
#' @param cell_names List of cell names to subset sensitivity data

FilterSensitivity <- function(sensitivity_database = NULL,
                              iqr_tolerance = 0.04,
                              n_coverage = 0.8,
                              cell_names = NULL){

  if("dplyr" %in% (.packages())==FALSE){library(dplyr)}

  #filter sensitivity by cell lines
  if(is.null(cell_names)==FALSE){
    sensitivity_database <- sensitivity_database[rownames(sensitivity_database) %in% cell_names,]
  }

  #filter by n_coverage
  na_cut_off <- nrow(sensitivity_database)*n_coverage
  na_cut_off <- apply(sensitivity_database, MARGIN = 2, FUN = function(x){sum(!is.na(x)) >= na_cut_off}) %>% unlist()

  sensitivity_database <- sensitivity_database[,na_cut_off]

  #filter by iqr
  iqr_cut_off <- apply(sensitivity_database, MARGIN = 2, FUN = function(x){
    out <- IQR(x, na.rm = T) >= iqr_tolerance
    if(is.na(out)){out <- FALSE}
    return(out)
  }) %>% unlist()

  sensitivity_database <- sensitivity_database[,iqr_cut_off]

  drugs <- colnames(sensitivity_database[,-1])
  return(drugs)
}

FilterSensitivity <- compiler::cmpfun(FilterSensitivity)
