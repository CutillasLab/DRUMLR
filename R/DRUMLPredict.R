#roxygen description
#' @name  predict.drug.sensitivity.based.on.distance.markers
#' @title  Make a heatmap of marker enrichment
#' @usage predict.drug.sensitivity.based.on.distance.markers(markers = GetMarkers("barasetib")$sensitive$.,inputdata = LoadData("PDBphos"),tissuefilter = "haematopoietic_and_lymphoid_tissue",metric = "aac", drug = "ABT-199")
#' @param predict.drug.sensitivity.based.on.distance.markers List of markers to filter the input data by
#' @param dataframe to use for the heatmap valuesdf.ppindex is phosphoproteomics, proteomics or transcriptomics data with rownames being identifiers
#' @param rmse.cutoff error cut off to define models as reliable
#' @param solid.or.aml which models to predict sensitvity
#' @param type "phospho","protein","RNA","ksea edges","reactome prot"

predict.drug.sensitivity.based.on.distance.markers <- function(df_distance,
                                                               df.ppindex,
                                                               rmse.cutoff=0.15,
                                                               solid.or.aml="aml",
                                                               model_type=c("phospho",
                                                                      "protein",
                                                                      "RNA"),
                                                               drugs){



  if("tidyverse" %in% (.packages())==FALSE){library(dplyr)}
  if("caret" %in% (.packages())==FALSE){library(caret)}
  if("h2o" %in% (.packages())==FALSE){library(h2o)}
  if("doParallel" %in% (.packages())==FALSE){library(doParallel)}
  if("foreach" %in% (.packages())==FALSE){library(foreach)}

  set.seed(123)
  h2o.init()

  if (input_type=="phospho"){

      if (model_type=="aml"){
      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/1wgs4djqpvhx6sr/Paths%20to%20models%20from%20phospho%20dist%20aml.csv?dl=1")
      # File where markers per drug are listed
      df_m <- DRUMLR:::phospho_aml_markers %>% data.frame(stringsAsFactors = F, row.names = 1)

      }else if(model_type=="solid"){

      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/zitaykih78g9osi/Paths%20to%20models%20from%20phospho%20dist%20solid.csv?dl=1")
      # File where markers per drug are listed
      df_m <- DRUMLR:::phospho_solid_markers %>% data.frame(stringsAsFactors = F, row.names = 1)    }


  }else if(input_type=="protein"){

    if(model_type=="aml"){
      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/sircuo6tveuox0g/Paths%20to%20models%20from%20proteomics%20dist%20aml.csv?dl=1")
      # File where markers per drug are listed
      df_m <- DRUMLR:::prot_aml_markers %>% data.frame(stringsAsFactors = F, row.names = 1)    }else if(solid.or.aml=="solid"){

      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/hh96xwxyunsfx8h/Paths%20to%20models%20from%20proteomics%20dist%20solid.csv?dl=1")
      # File where markers per drug are listed
      df_m <- DRUMLR:::prot_solid_markers %>% data.frame(stringsAsFactors = F, row.names = 1)    }
  }else if(input_type=="RNA"){
    if (model_type=="aml"){
      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/qao4dvr2ny5d6q0/Paths%20to%20models%20from%20rna%20dist%20aml.csv?dl=1")
      # File where markers per drug are listed
      df_m <- DRUMLR:::rna_aml_markers %>% data.frame(stringsAsFactors = F, row.names = 1)    }else if(solid.or.aml=="solid"){
      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/0vb5xco6vwt90v1/Paths%20to%20models%20from%20rna%20dist%20solid.csv?dl=1")
      # File where markers per drug are listed
      df_m <- DRUMLR:::rna_solid_markers %>% data.frame(stringsAsFactors = F, row.names = 1)
    }
  }

  file.models$paths <- gsub("/","\\",file.models$paths, fixed=T)
  file.models <- unique(file.models)

  # transpose the marker data so that it is in the right format for predictions
  df.tt <- data.frame(t(df_distance))

  drugs <- unique(file.models$drug)

  .PredictModelCaret <- function(file.models,model.name,drug, df_input){
    rmse.of.model.per.drug <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"val.rmse"][1]

    if(length(rmse.of.model.per.drug)>0|is.na(rmse.of.model.per.drug)==F|rmse.of.model.per.drug<rmse.cutoff){
      myfile <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"paths"][1]

      if(file.exists(myfile)){
        mymodel <- readRDS(myfile)
        tryCatch({out <- predict(mymodel,df_input)},error=function(e){})}
    }else{out<- data.frame(rep(NA, each = nrow(df_input)))}

    colnames(out)<- model.name
    rownames(out) <- rownames(df_input)
    return(out)
  }

  .PredictModelh2o <- function(file.models,model.name,drug, df_input){
    rmse.of.model.per.drug <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"val.rmse"][1]

    if(length(rmse.of.model.per.drug)>0|!is.na(rmse.of.model.per.drug)|rmse.of.model.per.drug<rmse.cutoff){
      myfile <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"paths"][1]

      if(file.exists(myfile)){
        mymodel <-h2o::h2o.loadModel(myfile)
        tryCatch({df.predicted.aac.dl[drug,] <- t(as.data.frame(predict(mymodel,as.h2o(df.t))))},error=function(e){})}
    }else{out<- data.frame(rep(NA, each = nrow(df_input)))}

    colnames(out)<- model.name
    rownames(out) <- rownames(df_input)
    return(out)
  }

  # Predict for each ML algorithm if the RMSE is less than the cut-off value
  df_pred <- foreach::foreach(i = drugs, combine = .rbind)%dopar%{

    dl <- .PredictModelh2o(file.models = file.models,
                           df_input = df_input,
                           drug = i,
                           model.name = "dl")

    nnet <- .PredictModelCaret(file.models = file.models,
                               df_input = df_input,
                               drug = i,
                               model.name = "nnet")

    pcr <- .PredictModelCaret(file.models = file.models,
                              df_input = df_input,
                              drug = i,
                              model.name = "pcr")

    pls <- .PredictModelCaret(file.models = file.models,
                             df_input = df_input,
                             drug = i,
                             model.name = "pls")

    rf <- .PredictModelCaret(file.models = file.models,
                             df_input = df_input,
                             drug = i,
                             model.name = "rf")

    svm <- .PredictModelCaret(file.models = file.models,
                             df_input = df_input,
                             drug = i,
                             model.name = "svm")

    out<- data.frame("drug" = drug,
                     "dl" = dl,
                     "nnet" = nnet,
                     "pcr" = pcr,
                     "pls" = pls,
                     "rf" = rf,
                     "svm" = svm)

    }



  return(list(distance.marker.data=df.distance.marker.data,
              predicted.based.on.nnet=df.predicted.aac.nnet,
              predicted.based.on.pls=df.predicted.aac.pls,
              predicted.based.on.rf=df.predicted.aac.rf,
              predicted.based.on.pcr=df.predicted.aac.pcr,
              predicted.based.on.svm=df.predicted.aac.svm,
              predicted.based.on.dl=df.predicted.aac.dl))
}
