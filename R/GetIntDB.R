# roxygen description
#' @title A score of marker ratios
#' @name GetIntMarkerDB
#' @export GetIntMarkerDB
#' @usage  GetIntMarkerDB("aml_prot_markers")
#' @description draw in internal databases
#' @param marker_database markers database to be drawn into environment. Can be path to csv marker file or the following:
#' \itemize{
#'   \item phospho_aml_markers
#'   \item phospho_solid_markers
#'   \item prot_aml_markers
#'   \item prot_solid_markers
#'   \item rna_aml_markers
#'   \item rna_solid_markers
#'   }

GetIntMarkerDB <- function(marker_database) {
  if ("dplyr" %in% (.packages()) == FALSE) {
    library(dplyr)
  }
  if (marker_database == "phospho_aml_markers") {
    out <-
      DRUMLR:::phospho_aml_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  } else if (marker_database == "phospho_solid_markers") {
    out <-
      DRUMLR:::phospho_solid_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  } else if (marker_database == "prot_aml_markers") {
    out <-
      DRUMLR:::prot_aml_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  } else if (marker_database == "prot_solid_markers") {
    out <-
      DRUMLR:::prot_solid_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  } else if (marker_database == "rna_aml_markers") {
    out <-
      DRUMLR:::rna_aml_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  } else if (marker_database == "rna_solid_markers") {
    out <-
      DRUMLR:::rna_solid_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else{
    out <-
      read.csv(marker_database, stringsAsFactors = F) %>% data.frame(row.names = 1, stringsAsFactors = F)

  }
  out <-
    list("sensitive_markers" = out[, c("m_sens", "sensitive_markers")],
         "resistant_markers" = out[, c("m_res", "resistant_markers")])

  return(out)
}
