#' Title Internal function for building machine learning models
#'
#' @title Internal Build model function
#' @name .BuildModel
#' @export .BuildModel
#'
#' @param df_train dataframe to use for model training. Response value must be in column named response.
#' @param partition.ratio partition ratio for making testing dataset
#' @param model_type string of model type to build e.g cubist, glm, nnet, pcr, pls, rf and svm
#' @param drug drug to use for naming model
#' @param input_type input type to use for naming model
#' @param save_path save file location
#' @param metric error metric to use to build model
#'
#'
#' @examples .BuildModel(df_train = df_distances, model_type = "rf", drug = "ABT-199", save_path = "my_model_dir")

.BuildModel <- function(df_train,
                        partition.ratio = 0.7,
                        make_train_set = T,
                        model_type,
                        drug,
                        input_type,
                        save_path,
                        metric = "RMSE") {
  #set seed
  set.seed(123)

  #set build start time
  build_start <- Sys.time()

  #if called make the training set
  if (make_train_set) {
    train.set <- caret::createDataPartition(df_train[, "response"],
                                            p = partition.ratio,
                                            list = FALSE,
                                            times = 1)
  }


  #make training dataset for model
  df.t <- df_train[train.set, ]

  if (model_type == "glm") {
    #build training control
    control <- trainControl(method = "repeatedcv",
                            number = 5,
                            repeats = 1)
    #build tuning grid
    grid <-
      expand.grid(C = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 5))

    #Build Bayes glm model
    model <- caret::train(
      response ~ .,
      data = df.t,
      method = 'bayesglm',
      #tunegrid=grid,
      trControl = control
    )

  } else if (model_type == "cubist") {
    #generate seeds
    seeds <- vector(mode = "list", length = nrow(df.t) + 1)
    seeds <- lapply(seeds, function(x)
      1:20)

    #build training control
    train_control <-
      caret::trainControl(method = "repeatedcv",
                          number = 3,
                          seeds = seeds)

    # cubist tuning grid
    grid <- expand.grid(committees = c(1, 10, 50, 100),
                        neighbors = c(0, 1, 5, 9))
    #build cubist model
    model <- caret::train(
      response ~ .,
      data = df.t,
      method = "cubist",
      metric = metric,
      trControl = train_control
    )
  } else if (model_type == "nnet") {
    #build training control
    control <- trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 3)

    #make tuning grid
    tunegrid <- expand.grid(expand.grid(
      size = c(1, 5, 10, 100),
      decay = c(0.1, 0.2, 0.5, 1)
    ))

    #make neural net model
    model <- train(
      response ~ .,
      data = df.t,
      method = "nnet",
      metric = metric,
      #tuneGrid=tunegrid,
      trControl = control,
      trace = F,
      maxit = 2000,
      linout = FALSE
    )

  } else if (model_type == "pcr") {
    #build training control
    train_control <- trainControl(method = "LOOCV")

    #build Principal Components regression model
    model <- train(
      response ~ .,
      data = df.t,
      method = "pcr",
      metric = metric,
      tuneLength = 7,
      maximize = TRUE,
      trControl = train_control
    )

  } else if (model_type == "pls") {
    #build training control
    train_control <- trainControl(method = "LOOCV")

    #build Partial List Squares model
    model <- train(
      response ~ .,
      data = df.t,
      method = "pls",
      metric = "RMSE",
      tuneLength = 7,
      maximize = TRUE,
      trControl = train_control
    )


  } else if (model_type == "rf") {
    #make training control
    control <-
      trainControl(method = "repeatedcv",
                   number = 10,
                   repeats = 3)

    #make tuning grid for random forests
    mtry <- c(1:5)
    tunegrid <- expand.grid(.mtry = mtry)
    #tunegrid <- expand.grid(.mtry=c(1:15), .ntree=c(1000, 1500, 2000, 2500))

    #build random forest model
    model <- train(
      response ~ .,
      data = df.t,
      method = "rf",
      metric = metric,
      tuneGrid = tunegrid,
      trControl = control,
      ntree = 1000,
      importance = TRUE
    )



  } else if (model_type == "svm") {
    #build training control
    control <-
      trainControl(method = "repeatedcv",
                   number = 10,
                   repeats = 3)

    #build tuning grid
    mtry <- c(1:5)
    tunegrid <- expand.grid(C = c(0, 0.000001, 0.00001,  0.0001))

    #build svm model
    model <- train(
      response ~ .,
      data = df.t,
      #method="svmRadial",
      method = "svmLinear",
      metric = metric,
      #tuneGrid=tunegrid,
      Cost = 0.0001,
      trControl = control,
      tuneLengh = 100
    )


  }

  #Get true responses
  df_test_resp <- df_train[-train.set[, 1], "response", drop = F]
  df_train_resp <- df_train[train.set[, 1], "response", drop = F]


  #predict responses using model
  df_test <- df_train[-train.set[, 1], -(ncol(df_train))]
  predicted_response_test <-
    predict(model, newdata = df_test[,!colnames(df.t) %in% "response"])
  predicted_response_train <-
    predict(model, newdata = df.t[,!colnames(df.t) %in% "response"])

  #output results into dataframe
  df_err_tr <-
    data.frame(
      sample = DRUMLR::RemoveRepeatNo(rownames(df.t)),
      measured_response = df_train_resp,
      predicted_response = predicted_response_train,
      data.set = "training"
    )

  df_err_ts <-
    data.frame(
      sample = DRUMLR::RemoveRepeatNo(rownames(df_test)),
      measured_response = df_test_resp,
      predicted_response = predicted_response_test,
      data.set = "testing"
    )

  df_results <- rbind(df_err_tr, df_err_ts)

  #rename columns
  colnames(df_results) <-
    c("sample"  ,  "measured_response", "predicted_response")

  #add additional error calculations
  df_results$accuracy <-
    (df_results$measured_response - df_results$predicted_response) / df_results$measured_response *
    100
  df_results$SE <-
    (df_results$measured_response - df_results$predicted_response) ^ 2
  df_results$AbsError <-
    abs(df_results$measured_response - df_results$predicted_response)


  #get errors data frame
  mod_rmse = data.frame(val.rmse = RMSE(df_err_ts[, "response"], df_err_ts[, "predicted_response"], na.rm = T))

  #save model in correct location
  model_name = paste(input_type, drug, model_type, sep = "_")
  model_path = paste(save_path, "/", model_name, ".rds", sep = "")
  saveRDS(object = model,
          file = model_path,
          compress = "xz")

  #build dataframe with model info
  model_info <- data.frame(
    model_path = model_path,
    model_type = model_type,
    input_type = input_type,
    drug = drug,
    rmse = mod_rmse
  )

  #return dataframes of model info, results table and results error
  out <- list(model_info, df_results)
  names(out) <- c("model_info", "model_results")

  comp_time <- paste(round(Sys.time() - build_start, 1), "seconds")
  print(paste(model_type, "completed for", drug, "in", comp_time))
  print(paste("model stored at:", model_path))
  print(paste("model error is", round(mod_rmse, 3)))


  return(out)

}

.BuildModel <- compiler::cmpfun(.BuildModel)
