#roxygen description
#' @name  PredictDRUML_example
#' @export PredictDRUML_example
#' @param df_input input dataframe
#' @param .scale if TRUE data will be scaled
#' @param input_type Input data type/model label
#' @param drugs Drug for which predictions will be made
#' @param models list of models for which you would like to use for predictions
#' @param model_path_file path to model path file which is created in model dir by build models
#' @param .marker_database markers database which will be used
#' @param shut_h2o if TRUE the h2o will be shutdown after predictions are made
#' @title  Make a heatmap of marker enrichment
#' @usage PredictDRUML(markers = GetMarkers("barasetib")$sensitive$.,inputdata = LoadData("PDBphos"),tissuefilter = "haematopoietic_and_lymphoid_tissue",metric = "aac", drug = "ABT-199")

PredictDRUML_example <- function(df_distance,
                         input_type,
                         drugs="all",
                         cancer_type,
                         models= "all",
                         models_dir,
                         shut_h2o =F,
                         computational_load = NULL){w

  if("tidyverse" %in% (.packages())==FALSE){library(dplyr)}
  if("caret" %in% (.packages())==FALSE){library(caret)}
  if("h2o" %in% (.packages())==FALSE){library(h2o)}
  if("doParallel" %in% (.packages())==FALSE){library(doParallel)}
  if("foreach" %in% (.packages())==FALSE){library(foreach)}

  if(!is.null(computational_load)) {
    cores <-
      as.integer((computational_load * parallel::detectCores()), length = 1)

  }else{
    cores <- 1
  }
  registerDoParallel(cores = cores)
  print(paste("running on", cores, "cores"))

  path_file <- DRUMLR:::DRUMLpaths

  if("all" %in% models){
    models <- c("bayes", "cubist", "nnet", "pcr", "pls", "rf", "svm", "dl")
  }
  if( "all" %in% drugs){
    drugs <- as.factor(path_file$drug) %>% levels() %>% paste()
  }

  #get model path for prediction
  models_path <- path_file[path_file$data_input == input_type&
                             path_file$model %in% models &
                             path_file$type %in% cancer_type &
                             path_file$drug %in% drugs,]

  models_path <- models_path[models_path$local_path != "not found", ]
  models_path$input_id <- paste(models_path$type,models_path$drug, sep = "_")
  models_dl <- models_path[models_path$model == "dl", ]
  models_all <- models_path[models_path$model != "dl", ]

  #generate distance values which will be used for prediction

  if(cancer_type == "solid"){
    rownames(df_distance) <- paste(rownames(df_distance),"st", sep = ".")
  }else{
    rownames(df_distance) <- paste(rownames(df_distance),cancer_type, sep = ".")}

  #reformat distance and get input names for prediction table labels
  df_distance <- df_distance %>% t() %>% data.frame()
  mod_inputs <- DRUMLR:::DRUMLMLinputs
  mod_inputs <- mod_inputs[mod_inputs$drug %in% drugs&
                             mod_inputs$data_input == input_type&
                             mod_inputs$type %in% cancer_type,]

  .trim_df <- function(markers){
    markers <- base::strsplit(x = markers, split = "-", fixed = T)
    markers <- unlist(markers)
    return(markers)
  }

  dist_list <- lapply(paste(mod_inputs$input_markers), FUN = .trim_df)
  names(dist_list) <- mod_inputs$input_id

  #predict by drawing models into environment sequentially
  predictions_all <- foreach(i = models_all$local_path, .combine = "cbind")%dopar%{
    #get model info and filter and scale df
    mod_imp <- models_all[models_all$local_path == i, "input_id"]
    df_mod <- DRUMLR::MinMaxNormalise(df_distance[,dist_list[[mod_imp]]],.margin=1)

    model <- readRDS(paste(models_dir, i, sep = ""))
    predicted_vals <- predict(model, df_mod)%>% data.frame()
    model_variabiles<- paste(models_all[models_all$local_path == i, c("drug", "type", "model")], collapse = "_")
    colnames(predicted_vals) <- model_variabiles
    return(predicted_vals)
  } %>% data.frame(stringsAsFactors = F)

   doParallel::stopImplicitCluster()
  #predict for dl models drawing models into environment sequentially
  #initiate h2o and make df_distance a h2o frame


  if("dl" %in% models){
    h2o.init(nthreads = cores)
    predictions_dl <- foreach(i = models_dl$local_path, .combine = "cbind")%do%{
      #get model info and filter and scale df
      mod_imp <- models_dl[models_dl$local_path == i, "input_id"]
      df_mod <- DRUMLR::MinMaxNormalise(df_distance[,dist_list[[mod_imp]]],.margin=1) %>% as.h2o()

      model <- h2o::h2o.import_mojo(paste(models_dir, i, sep = ""))
      predicted_vals <- h2o::h2o.predict(model, df_mod) %>% as.data.frame()
      model_variabiles <- paste(models_dl[models_dl$local_path == i, c("drug", "type", "model")], collapse = "_")
      colnames(predicted_vals) <- model_variabiles
      return(predicted_vals)
    } %>% data.frame(stringsAsFactors = F)

    #shutdown h2o instance
    if(shut_h2o==T){
    h2o.shutdown()
    "Y"}
    #combine predictions
    if(length(models)>1){
      predictions <- cbind(predictions_all, predictions_dl)

    }else{
      predictions <- predictions_dl
    }
  }else{
    predictions <- predictions_all
  }

  #return predictions
  return(predictions)
}
