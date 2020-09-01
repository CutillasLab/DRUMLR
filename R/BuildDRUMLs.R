#' Title
#' @title Build DRUML models for a list of drugs
#'
#' @param df_input Distance values to use for DRUML training.
#' @param .marker_database Marker database to make distance values with.
#' @param drugs list of drugs from df_response dataframe to generate models from if left NULL all column names in df_input will be used as drugs.
#' @param df_response Dataframe of response values with cell line in rownames and columnames as drugs.
#' @param models list of models to use for machine learning:
#' \itemize{
#'   \item glm = Baysian GLM model
#'   \item cubist = Cubist model
#'   \item pls = Partial list squares model
#'   \item pcr = Principal components regression
#'   \item nnet = Neural network model
#'   \item scm = Super vector machines
#'   \item rf = Random forests regressor model
#'   }
#' @param save_path save location for models and results csvs
#' @param save_csv if TRUE csv files containinf DRUML results will be made
#' @param computational_load Set as a decimal fraction of cores you want to use. If left as NULL only 1 core will be used.
#' @param method Statistical method for filtering distance markers for DRUML models default is spearman
#' @param top_pos_n maximum ammount of sensitive distance markers for DRUML models
#' @param top_neg_n  maximum ammount of resistant distance markers for DRUML models
#' @param order_var variable used to filter Distance models for DRUML:
#'  \itemize{
#'   \item p = pvalue
#'   \item r = rvalue
#'   \item slope = slope of fitted linear model
#'   }
#' @name BuildDRUMLs
#' @export BuildDRUMLs
#'
#' @examples

#df_input <- readxl::read_excel("C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data/Suppl_Dataset_2_Phosphoproteomics_of_48_cell_lines_in_triplicate.xlsx", sheet = "ppIndex") %>% data.frame(row.names = 1, stringsAsFactors = F)
#df_response <-read.csv("C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data/Suppl_Dataset_4_Drugs_with_normalised_AAC_and_class_info.csv", row.names = 1)
#.marker_database <- readxl::read_excel("C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data/Suppl_Dataset_7_sensitivity_markers_all.xlsx", sheet = "phospho aml") %>% data.frame(row.names = 1, stringsAsFactors = F)
#drugs = NULL
#input_type = "DRUML_model"
#save_path = "DRUML_models"
#save_csv = T
#scale = T
#computational_load = 0.8
#partition.ratio = 0.7
#models = c("glm", "cubist", "pls", "pcr", "nnet", "svm", "rf")
#corr_method = "spearman"
#top_pos_n = 10
#top_neg_n = 10
#order_var = "r"


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
                        models = c("glm", "cubist", "pls", "pcr", "nnet", "svm", "rf"),
                        corr_method = "spearman",
                        top_pos_n = 10,
                        top_neg_n = 10 ,
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

  df_response <- t(df_response) %>% data.frame(stringsAsFactors = F)

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

  #Build df_distance

  df_distance <- DRUMLR::DrugMarkerEnrichment(df = df_input,
                                              marker_database = .marker_database)


  #transpose distances and get cell lines
  if (scale) {
    df_distance <- base::apply(
      df_distance,
      MARGIN = 1,
      FUN = function(x) {
        x <- (x - min(x)) / (max(x) - min(x))
      }
    )
    print("Distance dataframe scaled")
  } else{
    df_distance <- t(df_distance) %>% data.frame(stringsAsFactors = F)
  }

  cell_lines <- DRUMLR::RemoveRepeatNo(rownames(df_distance))

  .BuildallDRUML <- function(.drug) {
    #add response data and remove na results
    df_train <- df_distance %>% data.frame()
    df_train[, "response"] <- df_response[cell_lines, .drug]
    df_train <- na.omit(df_train)

    #identify top 10 sensitive and

    markers <- DRUMLR::OptimiseMLInput(
      df_train = df_train,
      method = corr_method,
      top_pos_n = top_pos_n,
      top_neg_n = top_neg_n,
      partition.ratio = partition.ratio,
      return_corr_analysis = F,
      order_var = order_var
    )


    df_train <- df_train[, c(markers[["sensitive drug markers"]],
                             markers[["resistant drug markers"]],
                             "response")]
    #build models

    DRUML_mods <-
      lapply(
        models,
        FUN = function(x) {
          #DRUMLR::
            .BuildModel(
            df_train = df_train,
            partition.ratio = partition.ratio,
            make_train_set = T,
            model_type = x,
            drug = .drug,
            input_type = input_type,
            save_path = save_path,
          )
        }
      )

    #combine results for output
    n_list <- c(1:length(DRUML_mods))
    mod_infos <-
      lapply(
        n_list,
        FUN = function(x) {
          DRUML_mods[[x]]["model_info"]
        }
      ) %>% dplyr::bind_rows()
    mod_results <-
      lapply(
        n_list,
        FUN = function(x) {
          DRUML_mods[[x]]["model_results"]
        }
      ) %>% dplyr::bind_rows()

    mod_infos$res_markers <-
      paste(markers$`resistant drug markers`, collapse = ";")
    mod_infos$sens_markers <-
      paste(markers$`sensitive drug markers`, collapse = ";")

    output <- list(mod_infos, mod_results)
    names(output) <- c("mod_info", "mod_results")

    return(output)
  }


  print("Building models")
  if (is.null(computational_load)) {
    DRUML_output <-
      foreach::foreach(
        i = drugs,
        .combine = "c",
        .packages = c(
          "dplyr",
          "caret",
          "Cubist",
          "pls",
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
        .combine = "c",
        .packages = c(
          "dplyr",
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

  #organise foreach loop outputs into 2 dataframse
  DRUML_infos <-
    DRUML_output[names(DRUML_output) %in% "mod_info"] %>% dplyr::bind_rows()
  DRUML_results <-
    DRUML_output[names(DRUML_output) %in% "mod_results"] %>% dplyr::bind_rows()

  #make names for path and result outputs
  paths <- paste(input_type, "model_paths", sep = "_")
  results <- paste(input_type, "model_results", sep = "_")

  #write csv files and place in save directory
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

  return(output)

}
