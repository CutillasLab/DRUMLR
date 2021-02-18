#' Title
#' @title Build DRUML models for a list of drugs
#'
#' @param df_input dataframe of input data.
#' @param .marker_database Markers database which will be used to for Drug Marker Enrichment
#' @param drugs list of drugs from which models will be generated. If left NULL models will be built for all drugs in the marker database.
#' @param df_response Dataframe of drug responses(AAC) with cell lines as row names and column names as drugs.
#' @param models list of models to use for machine learning:
#' \itemize{
#'   \item glm = Baysian GLM model
#'   \item cubist = Cubist model
#'   \item pls = Partial list squares model
#'   \item pcr = Principal components regression
#'   \item nnet = Neural network model
#'   \item scm = Super vector machines
#'   \item rf = Random forests regressor model
#'   \item dl = Deep learning model using h2o package (needs Java to be up to date)
#'   }
#' @param save_path save location for models and results csvs
#' @param save_csv if TRUE csv files containing DRUML results will be made
#' @param computational_load Set as a decimal fraction of the number of cores you which to be recruited for model building e.g a value of 0.8 means 80\% of all available CPU cores will be used to build models. If left as NULL only 1 core will be used.
#' @param method Statistical method for filtering distance markers for DRUML models default is spearman
#' @param max_n_D maximum number of markers to use
#' @param min_n_D Minimum number of  markers to use
#' @param order_var variable used to filter Distance models for DRUML:
#'  \itemize{
#'   \item p = p value
#'   \item r = r value
#'   \item slope = slope of fitted linear model
#'   }
#' @name BuildDRUMLs
#' @export BuildDRUMLs
#'
#' @examples

BuildDRUMLs <- function(df_input,
                        .marker_database,
                        df_response,
                        drugs = NULL,
                        input_type = "DRUML_model",
                        save_path = "DRUML_models",
                        save_csv = T,
                        scale = T,
                        computational_load = 0.8,
                        partition.ratio = 0.7,
                        models = c("glm", "cubist", "pls", "pcr", "nnet", "svm", "rf", "dl"),
                        corr_method = "spearman",
                        max_n_D = 30,
                        min_n_D = 7 ,
                        order_var = "r") {
  if ("dplyr" %in% (.packages()) == FALSE) {
    library(dplyr)
  }
  if ("caret" %in% (.packages()) == FALSE) {
    library(caret)
  }
  if ("h2o" %in% (.packages()) == FALSE) {
    library(h2o)
  }

  if ("foreach" %in% (.packages()) == FALSE) {
    library(foreach)
  }

  if ("doParallel" %in% (.packages()) == FALSE) {
    library(doParallel)
  }

  if (is.null(drugs)) {
    drugs <- colnames(df_response)
  }

  #Initiate multicore processing if specified
  if (!is.null(computational_load)) {
    cores <-
      as.integer((computational_load * parallel::detectCores()), length = 1)
    registerDoParallel(cores = cores)
    print(paste("running on", cores, "cores"))
  }

  #make save_path directory if required
  if (!dir.exists(save_path)) {
    dir.create(save_path)
    print(paste(save_path, "created"))
  }

  #Generate distance values for input data
  df_distance <- DRUMLR::DrugMarkerEnrichment(df = df_input,
                                                  marker_database = .marker_database)

  #transpose distances for ML models
  df_distance <- t(df_distance) %>% data.frame(stringsAsFactors = F)

  #Get list of cell lines by removing repeatnumber
  cell_lines <- DRUMLR::RemoveRepeatNo(rownames(df_distance))

  #.BuildallDRUML is an internal function for building ML models for each drug
  .BuildallDRUML <- function(.drug) {
    #add response data and remove na results
    df_input <- df_distance %>% data.frame()
    response <- df_response[cell_lines, .drug, drop = F]
    df_input[, "response"] <- response
    df_input <- na.omit(df_input)
    response <- df_input[, "response"]
    #identify top 10 sensitive and resistant markers

    #optimise the distance values which will be used as a ML input
    op_markers <- DRUMLR::OptimiseMLInput(
      df_train = df_input,
      method = corr_method,
      max_n_D = max_n_D,
      min_n_D = min_n_D,
      partition.ratio = partition.ratio,
      return_corr_analysis = F,
      order_var = order_var
    )

    print(paste(.drug, "markers optimised"))

    #combine markers to build training input data frame
    df_input <- df_input[, c(op_markers[["sensitive drug markers"]],
                             op_markers[["resistant drug markers"]],
                             "response")]

    #Scale the inputs by cell line
    df_input <- df_input[,colnames(df_input)!="response"] %>% DRUMLR::MinMaxNormalise(.margin = 1)
    df_input[, "response"] <- response

    #split data into training and testing subsets
    train.set <- caret::createDataPartition(df_input[, "response"],
                                            p = partition.ratio,
                                            list = FALSE,
                                            times = 1)


    #make training dataset for model
    df_train <- df_input[train.set, ]
    df_test <-df_input[-train.set[, 1],]

    #build models using the .BuildModel function
    DRUML_mods <-
      lapply(
        models,
        FUN = function(x) {
          #DRUMLR::
          DRUMLR::.BuildModel(
            df_train = df_train,
            df_test = df_test,
            model_type = x,
            drug = .drug,
            input_type = input_type,
            save_path = save_path,
          )
        }
      )


    #Reorganize ML building output to build performance results csvs DRUML models
    n_list <- c(1:length(DRUML_mods))
    mod_infos <-
      lapply(
        n_list,
        FUN = function(x) {
          DRUML_mods[[x]]["model_info"][[1]]
        }
      ) %>% dplyr::bind_rows() %>% data.frame(stringsAsFactors = F)
    mod_results <-
      lapply(
        n_list,
        FUN = function(x) {
          DRUML_mods[[x]]["model_results"][[1]]
        }
      ) %>% dplyr::bind_rows() %>% data.frame(stringsAsFactors = F)

    colnames(mod_infos) <- gsub(x = colnames(mod_infos), pattern = "model_info.",replacement =  "", fixed = T)
    colnames(mod_results) <- gsub(x = colnames(mod_results), pattern = "model_results.",replacement =  "", fixed = T)

    mod_infos$res_markers <-
      paste(op_markers$`resistant drug markers`, collapse = ";") %>% unlist()
    mod_infos$sens_markers <-
      paste(op_markers$`sensitive drug markers`, collapse = ";") %>% unlist()


    output <- list(mod_infos, mod_results)
    names(output) <- c("mod_info", "mod_results")

    return(output)
  }

  #if dl models are being h2o will be initiated
  if ("dl" %in% models){
    h2o.init()
  }

  print("Building models")
  #drugs <- drugs[1:5]
  if (is.null(computational_load)) {
    DRUML_output <-
      foreach::foreach(
        i = drugs,
        .inorder = T,
        .errorhandling = "pass",
        .combine = "c",
        .packages = c(
          "dplyr",
          "caret",
          "Cubist",
          "pls",
          "h2o",
          "glmnet",
          "kernlab",
          "randomForest"
        )
      ) %do% {
        out <- .BuildallDRUML(.drug = i)
        return(out)
      }
  } else{
    DRUML_output <-
      foreach::foreach(
        i = drugs,
        .inorder = T,
        .errorhandling = "pass",
        .combine = "c",
        .packages = c(
          "dplyr",
          "arm",
          "h2o",
          "caret",
          "Cubist",
          "pls",
          "glmnet",
          "kernlab",
          "randomForest"
        )
      ) %dopar% {
        out <- .BuildallDRUML(.drug = i)
        return(out)
      }
  }

  print("models built")

  #Organise DRUML output into results and paths dfs
  DRUML_infos <-
    DRUML_output[names(DRUML_output) %in% "mod_info"] %>% dplyr::bind_rows(.id = "") %>% data.frame(stringsAsFactors = F)
  DRUML_results <-
    DRUML_output[names(DRUML_output) %in% "mod_results"] %>% dplyr::bind_rows(.id = "") %>% data.frame(stringsAsFactors = F)


  #make names for path and result outputs
  paths <- paste(input_type, "model_paths", sep = "_")
  results <- paste(input_type, "model_results", sep = "_")


  #write csv files and save in model directory
  if (save_csv) {
    write.csv(file = paste(save_path, paste(paths, "csv", sep = "."), sep = "/"), x = DRUML_infos)
    write.csv(file = paste(save_path, paste(results, "csv", sep = "."), sep = "/"), x = DRUML_results)
    print(paste("csv files saved in", save_path))
  }

  output <- list(DRUML_infos, DRUML_results)
  names(output) <- c(paths, results)

  #close mulitcore processing
  if (is.null(computational_load) == F) {
    doParallel::stopImplicitCluster()
  }

  #if DL models were built user will be given the option to close h2o connection
  if ("dl" %in% models){
    h2o.shutdown()

  }

  return(output)
}

