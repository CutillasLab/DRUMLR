# Emperical Markers of Sensitivity and resistantce
#' @name FilterSensitivity
#' @export FilterSensitivity
#' @title Filter Drugs by sensitivity data
#' @param sensitivity_database sensitvitiy data with cell lines as rownamed and drugs as columns
#' @param cell_annotations dataframe with cell line annotations cell names must be passed through make.names function and placed in "cell.name" column. Tumour type must be stored in "unique.tissue" column.
#' @param cell_lines Cell lines present in omics data
#' @param tissue_type Tissue types to filter sensitvity by the default is those used in DRUML panel ("haematopoietic_and_lymphoid_tissue", "liver", "oesophagus")
#' @param iqr_tolerance IQR range cut off to filter drug list
#' @param n_coverage n_coverage cut off required to ensure models are reliable

FilterSensitivity <- function(sensitivity_database = NULL,
                             cell_annotations = NULL,
                             cell_lines = NULL,
                             tissue_type = c("AML", "Hepatocellular", "Oesophagus"),
                             iqr_tolerance = 0.04,
                             n_coverage = 0.8){
  require(dplyr)

  if(is.null(cell_annotations)){cell_annotations <- DRUMLR:::DRUMLcellinfo %>% data.frame()}
  if(is.null(sensitivity_database)){sensitivity_database <- DRUMLR:::DRUMLaac %>% data.frame()}
  if(is.null(cell_lines)){
    cell_lines <- colnames(DRUMLR:::DRUMLphos)
    cell_lines <- RemoveRepeatNo(cell_lines)
  }else{
    cell_lines <- RemoveRepeatNo(cell_lines)
    }
  if(is.null(tissue_type)){errorCondition(message = "please add a list of cancer tissue types to use for analysis: those available in the PharmacoDB database are autonomic_ganglia, biliary_tract, bone, breast, central_nervous_system, cervix, endometrium, haematopoietic_and_lymphoid_tissue, kidney, large_intestine, liver, lung, oesophagus, ovary, pancreas, pleura, prostate, skin, soft_tissue, stomach, thyroid, upper_aerodigestive_tract, urinary_tract")
                                          stop()
                                          }

  #filter cell lines by tissue
  cell_names1 <- cell_annotations[cell_annotations$unique.tissue %in% tissue_type, "Cell.name"]
  cell_names <- cell_names1[cell_names1 %in% cell_lines]

  #filter sensitivity by cell lines

  sensitivity_database <- sensitivity_database[rownames(sensitivity_database) %in% cell_names,]

  #filter by n_coverage
    na_cut_off <- nrow(sensitivity_database)*n_coverage
    na_cut_off <- apply(sensitivity_database, MARGIN = 2, FUN = function(x){sum(!is.na(x)) >= na_cut_off}) %>% unlist()

    sensitivity_database <- sensitivity_database[,na_cut_off]

    #filter by iqr
    x <- sensitivity_database[,2]
    iqr_cut_off <- apply(sensitivity_database, MARGIN = 2, FUN = function(x){
      out <- IQR(x, na.rm = T) >= iqr_tolerance
      if(is.na(out)){out <- FALSE}
      return(out)
      }) %>% unlist()

    sensitivity_database <- sensitivity_database[,iqr_cut_off]

    drugs <- list(colnames(sensitivity_database[,-1]),cell_names, cell_names1)
    names(drugs)<- c("drugs", "cell_names", "all_tissue_names")
    return(drugs)
}
