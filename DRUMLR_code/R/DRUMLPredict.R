#roxygen description
#' @name  PredictDRUML
#' @export PredictDRUML
#' @param df_distance dataframe of caluclated D values of input data made using DRUMLR::DrugMarkerEnrichment()
#' @param drugs Drug for which predictions will be made
#' @param models list of models for which you would like to use for predictions
#' @param path_file path file containing model information. If none is specified the default path file save location will be used.
#' @title  Make a heatmap of marker enrichment
#' @usage PredictDRUML(df_distance = validation_data, drugs = c("BYL.719", "ABT.119"), models= c("rf", "dl"), models_dir = "~/Downloads/Models", path_file = df_path)

PredictDRUML <- function(df_distance,
                         drugs="all",
                         models= "all",
                         models_dir,
                         path_file = NULL,
                         computational_load = NULL){

  if("tidyverse" %in% (.packages())==FALSE){library(dplyr)}
  if("caret" %in% (.packages())==FALSE){library(caret)}
  if("h2o" %in% (.packages())==FALSE){library(h2o)}
  if("doParallel" %in% (.packages())==FALSE){library(doParallel)}
  if("foreach" %in% (.packages())==FALSE){library(foreach)}

  #get number of cores for parallel processing
  if(!is.null(computational_load)) {
    cores <-
      as.integer((computational_load * parallel::detectCores()), length = 1)

  }else{
    cores <- 1
  }
  registerDoParallel(cores = cores)
  print(paste("running on", cores, "cores"))

  #find path file if none is specified

  if(is.null(path_file)){
    path_file <- list.files(models_dir) %>% grep(pattern = "model_paths.csv", fixed = T, value = T)
    path_file <- read.csv(paste(models_dir, path_file, sep = "/"),row.names = 1, stringsAsFactors = F)

  }

  if(models == "all"){
    models <- c("bayes", "cubist", "nnet", "pcr", "pls", "rf", "svm", "dl")
  }
  if(drugs == "all"){
    drugs <- as.factor(path_file$drug) %>% levels() %>% paste()
  }

  #get model path for prediction
  models_path <- path_file[path_file$model_type %in% models &
                           path_file$drug %in% drugs,]

  models_dl <- models_path[models_path$model_type == "dl", ]
  models_all <- models_path[models_path$model_type != "dl", ]

  #generate distance values which will be used for prediction
    .trim_df <- function(markers){
    markers <- base::strsplit(x = markers, split = ";", fixed = T)
    markers <- unlist(markers)
    return(markers)
  }

  dist_list <- lapply(paste(models_path$res_markers, models_path$sens_markers, sep = ";"), FUN = .trim_df)
  names(dist_list) <- models_path$model_path
  colnames(models_all)

  #reformat distance and get input names for prediction table labels
  df_distance <- df_distance %>% t() %>% data.frame()

  #predict by drawing models into environment sequentially
  predictions_all <- foreach(i = models_all$model_path, .combine = "cbind")%dopar%{
    model <- readRDS(paste(models_dir, i, sep = "/"))
    mod_input <- dist_list[[i]]
    df_mod <- DRUMLR::MinMaxNormalise(df_distance[,mod_input],.margin=2)
    predicted_vals <- predict(model, df_distance)%>% data.frame()
    rownames(predicted_vals) <- rownames(df_distance)
    colnames(predicted_vals) <- paste(models_all[models_all$model_path == i, c("drug", "input_type", "model_type")], collapse = "_")
    return(predicted_vals)
  }

  #terminate doparallel cluster
  doParallel::stopImplicitCluster()


  #predict for dl models drawing models into environment sequentially

  #initiate h2o and make df_distance a h2o frame
  if("dl" %in% models){
    #run h2o on the same numer of cores as doparallel
    h2o.init(nthreads = cores)
    predictions_dl <- foreach(i = models_dl$model_path, .combine = "cbind")%do%{
      #get model info and filter and scale df
      mod_imp <- dist_list[[i]]
      df_mod <- DRUMLR::MinMaxNormalise(df_distance[,dist_list[[mod_imp]]],.margin=2) %>% as.h2o()

      model <- h2o::h2o.import_mojo(paste(models_dir, i, sep = ""))
      predicted_vals <- h2o::h2o.predict(model, df_mod) %>% as.data.frame()
      model_variabiles <- paste(models_dl[models_dl$model_path == i, c("drug", "input_type", "model_type")], collapse = "_")
      colnames(predicted_vals) <- model_variabiles
      return(predicted_vals)
    } %>% data.frame(stringsAsFactors = F)

    #shutdown h2o instance
    if(shut_h2o==T){
      h2o::h2o.shutdown()
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

}

