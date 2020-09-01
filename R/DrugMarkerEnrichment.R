#Building models with different markers
# roxygen description
#' @title Generate D values from input data and marker database
#' @name DrugMarkerEnrichment
#' @export DrugMarkerEnrichment
#' @usage  DrugMarkerEnrichment(df = df.ppindex, marker_database = "aml_prot_markers", scale = T, output = c("Distance", "zscore"))
#' @description Get Distance values or zscores between sensitive and resistant markers
#' @param df input data with cell lines as colnames and varaibles as rownames
#' @param marker_database markers database to be used can be the following:
#' \itemize{
#'   \item phospho aml
#'   \item phospho solid_markers
#'   \item prot aml
#'   \item prot solid
#'   \item rna aml
#'   \item rna solid
#'   }
#'  @param marker_database_path path of the directory containing the DRUMLR marker file
#' @param scale Set as TRUE to scale input data
#' @param output The output metric you want to use can be "Distance" or "zscore"

#marker_database_path = "C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data"

DrugMarkerEnrichment <-
  function(df,
           marker_database,
           marker_database_path) {

    if ("dplyr" %in% (.packages()) == FALSE) {
      library(dplyr)
    }

    if ("sites" %in% colnames(df) == F) {
      dat_edges <- paste(rownames(df))
    } else{
      dat_edges <- paste(df$sites)
      df <- df[, colnames(df) != "sites"]
    }

    #dat_edges <- strsplit(dat_edges, ";") %>% unlist() %>% unique()

    #get markers if using DRUML markers

    if(marker_database %in% c("phospho aml", "phospho solid", "prot aml", "prot solid", "rna aml", "rna solid")){

      marker_file = paste(marker_database_path, "Suppl_Dataset_7_sensitivity_markers_all.xlsx", sep = "/")

      dbs <- readxl::read_excel(marker_file, sheet = marker_database)

    }else{dbs  <- marker_database}

    #split marker edges into list
    resistant_markers <- strsplit(as.character(dbs$resistant_markers),"-")
    names(resistant_markers) <- dbs$drugs

    colnames(dbs)
    #split marker edges into list
    sensitive_markers <- strsplit(as.character(dbs$sensitive_markers),"-")
    names(sensitive_markers) <- dbs$drugs

    #filter edges by rownames retain only nodes with more than 3 edges

    sensitive_markers <- lapply(sensitive_markers, FUN= function(x){x <- x[x %in%dat_edges]})
    sensitive_markers <- sensitive_markers[lengths(sensitive_markers)>=3]

    resistant_markers <- lapply(resistant_markers, FUN= function(x){x <- x[x %in%dat_edges]})
    resistant_markers <- resistant_markers[lengths(resistant_markers)>=3]

    #ensure that values for resistant and sensitve markers are present
    com_markers <- dbs$drugs[dbs$drugs %in% names(sensitive_markers)&dbs$drugs %in% names(resistant_markers)]
    sensitive_markers <- sensitive_markers[com_markers]
    resistant_markers <- resistant_markers[com_markers]

    #function for getting inputs
    .GetDistInputs <- function(x){
      df <-df[x,]
      median <- apply(df, MARGIN=2, function(xx){median(xx, na.rm = T)})
      Q3 <- apply(df, MARGIN=2, function(xx){quantile(xx,probs = 0.75, na.rm = T)})
      out <- data.frame("D"= Q3 - median)
      return(out)
    }

    print("calculating resistance values")
    resistant_inputs <- lapply(resistant_markers, .GetDistInputs)

    print("calculating sensitivity values")
    sensitive_inputs <- lapply(sensitive_markers, .GetDistInputs)

    out <- lapply(com_markers, FUN = function(x){
      df <- sensitive_inputs[[x]]-resistant_inputs[x]
      colnames(df) <- x
      return(df)
    })

    out <- data.frame(out) %>%t()

    print("Drug Enrichment Complete")
    return(out)

  }

