# roxygen description
#' @title A score of marker ratios
#' @name GetIntMarkerDB
#' @export GetIntMarkerDB
#' @usage  GetIntMarkerDB("aml_prot_markers")
#' @description draw in internal databases
#' @param database database to draw into environment

GetIntMarkerDB <- function(database){
  if(database == "aml_phos_markers"){out <- DRUMLR:::aml_phos_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "all_phos_markers"){out <- DRUMLR:::all_phos_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "solid_phos_markers"){out <- DRUMLR:::solid_phos_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "aml_prot_markers"){out <- DRUMLR:::aml_prot_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "all_prot_markers"){out <- DRUMLR:::all_prot_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "solid_prot_markers"){out <- DRUMLR:::solid_prot_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "aml_rna_markers"){out <- DRUMLR:::aml_rna_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "all_rna_markers"){out <- DRUMLR:::all_rna_markers %>% data.frame(row.names = 1, stringsAsFactors = F)
  }else if(database == "solid_rna_markers"){out <- DRUMLR:::solid_rna_markers %>% data.frame(row.names = 1, stringsAsFactors = F)}

  out <- list("sensitive_markers" = out[,c("m_sens", "sensitive_markers")],
              "resistant_markers" = out[,c("m_res", "resistant_markers")])

  return(out)
}
