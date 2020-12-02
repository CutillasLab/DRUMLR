#' Title Internal function for building machine learning models
#' @title Internal Build model function
#' @name .BuildModel
#' @export .BuildModel
#' @param df_input dataframe to use for model training. Response value must be in column named response.
#' @param partition.ratio partition ratio for making testing dataset
#' @param model_type string of model type to build e.g cubist, glm, nnet, pcr, pls, rf and svm
#' @param drug drug to use for naming model
#' @param input_type input type to use for naming model
#' @param save_path save file location
#' @param metric error metric to use to build model
#' @examples .BuildModel(df_input = df_distances, model_type = "rf", drug = "ABT-199", save_path = "my_model_dir")


.BuildModel <- function(df_train,
                        df_test,
                        model_type,
                        drug,
                        input_type,
                        save_path,
                        metric = "RMSE",
                        initiate_h2o = T) {
  #set seed
  set.seed(123)

  #set build start time
  build_start <- Sys.time()

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
      data = df_train,
      method = 'bayesglm',
      #tunegrid=grid,
      trControl = control
    )

  } else if (model_type == "cubist") {
    #generate seeds
    seeds <- vector(mode = "list", length = nrow(df_train) + 1)
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
      data = df_train,
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
      data = df_train,
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
      data = df_train,
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
      data = df_train,
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
      data = df_train,
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
      data = df_train,
      #method="svmRadial",
      method = "svmLinear",
      metric = metric,
      #tuneGrid=tunegrid,
      Cost = 0.0001,
      trControl = control,
      tuneLengh = 100
    )


  }else if(model_type == "dl"){
    if(initiate_h2o==T){
      h2o::h2o.init()
    }
    #convert training and testing df to h2o frames
    df_train.h2o <- df_train %>% as.h2o()
    df_test.h2o <- df_test  %>% as.h2o()
    #build model
    model = h2o.deeplearning(x = setdiff(colnames(df_train),
                                         c("response")),
                             y = "response",
                             training_frame=df_train.h2o,
                             validation_frame=df_test.h2o,
                             #hidden=c(32,32,32),                  ## small network, runs faster
                             epochs=10000,                      ## hopefully converges earlier...
                             score_validation_samples=10000,      ## sample the validation dataset (faster)
                             stopping_rounds=2,
                             stopping_metric="RMSE",
                             stopping_tolerance=0.01,
                             variable_importances=T)
  }

  #Get true responses
  df_test_resp <- df_test[, "response", drop = F]
  df_train_resp <- df_train[, "response", drop = F]
  df_test <- df_test[,!colnames(df_test) %in% "response"]
  df_train <- df_train[,!colnames(df_train) %in% "response"]

  #make model save names
  model_name = paste(input_type, drug, model_type, sep = "_")
  model_path = paste(save_path, "/", model_name, ".rds", sep = "")

  #predict and save models- h2o models require different code to all others
  if(model_type=="dl"){

    #make tested dataframe
    df_test.h2o <- as.h2o(df_test)
    df_train.h2o <-as.h2o(df_train)

    #predict responses using model
    predicted_response_test <-
      predict(model, newdata = df_test.h2o) %>% as.data.frame()
    predicted_response_train <-
      predict(model, newdata = df_train.h2o) %>% as.data.frame()

    #save DL model
    dest <- h2o::h2o.save_mojo(
      object = model,
      path = model_path)

    file.rename(dest, paste(save_path, "/", model_name, ".zip", sep = ""))
    file.remove(model_path)

  }else{

  #predict responses using model
  predicted_response_test <-
    predict(model, newdata = df_test[,!colnames(df_train) %in% "response"])
  predicted_response_train <-
    predict(model, newdata = df_train[,!colnames(df_train) %in% "response"])

  #save model in correct location
  saveRDS(object = model,
          file = model_path,
          compress = "xz")
  }

  #output results into dataframe
  df_err_tr <-
    data.frame(
      sample = DRUMLR::RemoveRepeatNo(rownames(df_train)),
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
  colnames(df_err_tr) <- c("sample", "measured_response", "predicted_response", "data.set")
  colnames(df_err_ts) <- c("sample", "measured_response", "predicted_response", "data.set")
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
  df_results$model <- model_type


  #get errors data frame
  mod_rmse = data.frame(val.rmse = RMSE(df_err_ts[, "measured_response"], df_err_ts[, "predicted_response"], na.rm = T))

  #build dataframe with model info
  model_info <- data.frame(
    model_path = model_name,
    model_type = model_type,
    input_type = input_type,
    drug = drug,
    rmse = mod_rmse
  )

  #return dataframes of model info, results table and results error
  out <- list(model_info, df_results)
  names(out) <- c("model_info", "model_results")

  #print model information
  comp_time <- paste(round(Sys.time() - build_start, 1), "seconds")
  print(paste(model_type, "completed for", drug, "in", comp_time))
  print(paste("model stored at:", model_path))
  print(paste("model error is", round(mod_rmse, 3)))


  return(out)

}
