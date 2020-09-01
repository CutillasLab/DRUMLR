#roxygen description
#' @name  PredictDRUML
#' @export PredictDRUML
#' @param df_input input dataframe
#' @param .scale if TRUE data will be scaled
#' @param input_type Input data type/model label
#' @param drug Drug for which predictions will be made
#' @param models list of models for which you would like to use for predictions
#' @param model_path_file path to model path file which is created in model dir by build models
#' @param .marker_database markers database which will be used to
#' @title  Make a heatmap of marker enrichment
#' @usage PredictDRUML(markers = GetMarkers("barasetib")$sensitive$.,inputdata = LoadData("PDBphos"),tissuefilter = "haematopoietic_and_lymphoid_tissue",metric = "aac", drug = "ABT-199")


#test_datasets
#model_path_file <- "C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Models/model_paths.csv"
#input_type <- "phospho"
#drugs <- "barasertib"
#cancer_type <- "aml"
#models <- c("svm", "pls")
#path_file <- "DRUML"
#models <- NULL
 # colnames(df_input) <- paste(colnames(df_input), cancer_type, sep = ".")

#df_input <- DrugMarkerEnrichment(df = df_backup, marker_database = "phospho aml", marker_database_path = "C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data")

#make DRUMLR_data path csv internal

PredictDRUML <- function(df_input,
                         input_type,
                         drugs,
                         cancer_type,
                         models,
                         .marker_database,
                         models_paths_path){

  if("tidyverse" %in% (.packages())==FALSE){library(dplyr)}
  if("caret" %in% (.packages())==FALSE){library(caret)}
  if("h2o" %in% (.packages())==FALSE){library(h2o)}
  if("doParallel" %in% (.packages())==FALSE){library(doParallel)}
  if("foreach" %in% (.packages())==FALSE){library(foreach)}

  if(models_paths_path == "DRUML"){
    path_file <- read.csv("C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Models/model_paths.csv", stringsAsFactors = F)
  }else{
    path_file <- read.csv(models_dir, stringsAsFactors = F)
  }

  if(is.null(models)){
    models <- c("bayes", "cubist", "nnet", "pcr", "pls", "rf", "svm")
  }

  #get model path for prediction

  models_path <- path_file[path_file$data_input == input_type&
                             path_file$model %in% models &
                             path_file$type %in% cancer_type &
                             path_file$drug %in% drugs,]


  #generate distance values which will be used for prediction
  df_distance <- df_input
  rownames(df_distance) <- paste(rownames(df_input),cancer_type, sep = ".")

  #median imputation of missing data
  df_distance[is.na(df_distance2)]<- 0

  rownames(df_distance2) <- filter_variables

  #predict by drawing models into environment sequentially
  predictions <- foreach(i = models_path$local_path, .combine = "cbind")%do%{
    if(models_dir == "DRUMLR"){
      model <- readRDS(url(models_path$github_link_actual))
    }else{
      model <- readRDS(paste(models_dir, i, sep = "/"))}

    predicted_vals <- predict(model, t(df_distance))%>% data.frame()
    rownames(predicted_vals) <- colnames(df_distance)
    colnames(predicted_vals) <- paste(models_path[models_path$local_path == i, c("drug", "type", "model")], collapse = "_")
    return(predicted_vals)
  }

  #return predicions
  return(predictions)
}
