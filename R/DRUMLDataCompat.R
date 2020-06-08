#roxygen description
#' @name  TestExternalDataComp
#' @title  Recaluclate Error paramaters for external data to account for missing variables
#' @usage TestExternalDataComp(markers = GetMarkers("barasetib")$sensitive$.,inputdata = LoadData("PDBphos"),tissuefilter = "haematopoietic_and_lymphoid_tissue",metric = "aac", drug = "ABT-199")
#' @param TestExternalDataComp List of markers to filter the input data by
#' @param dataframe to use for the heatmap valuesdf.ppindex is phosphoproteomics, proteomics or transcriptomics data with rownames being identifiers
#' @param rmse.cutoff error cut off to define models as reliable
#' @param solid.or.aml which models to predict sensitvity
#' @param type "phospho","protein","RNA","ksea edges","reactome prot"

#missing_pep <- runif(n = 10000, min = 1, max = 22804)%>% round(digits = 0)
#df_input <- DRUMLR:::DRUMLphos[missing_pep,info$cell_names]
#model_type<- "solid"
#input_type <- "phospho"

TestExternalDataComp <- function(df_input,
                                 model_type,
                                 input_type,
                                 drugs = NULL
                            ){
  require(tidyverse)
  require(caret)
  require(doParallel)
  require(caret)
  require(h2o)

  #Get Marker, Validation dataset and models list

  if(model_type == "solid"){
    info <- DRUMLR::FilterSensitivity(sensitivity_database = DRUMLR:::DRUMLaac,
                              cell_annotations = DRUMLR:::DRUMLcellinfo ,
                              cell_lines = colnames(df_input),
                              tissue_type = c("liver", "oesophagus"),)

  }else if(model_type == "aml"){
    info <- DRUMLR::FilterSensitivity(sensitivity_database = DRUMLR:::DRUMLaac,
                              cell_annotations = DRUMLR:::DRUMLcellinfo ,
                              cell_lines = colnames(df_input),
                              tissue_type = "haematopoietic_and_lymphoid_tissue",)
    }

  if (input_type=="phospho"){
    df_val <- DRUMLR:::DRUMLphos[,info$cell_names]
    if (model_type=="aml"){
      # File were model names and paths are located:
      file.models <- read.csv(
        "https://www.dropbox.com/s/1wgs4djqpvhx6sr/Paths%20to%20models%20from%20phospho%20dist%20aml.csv?dl=1")
      # import markers
      df_m <- "phospho_aml_markers"

      }else if(model_type=="solid"){
        # File were model names and paths are located:
        file.models <- read.csv(
          "https://www.dropbox.com/s/zitaykih78g9osi/Paths%20to%20models%20from%20phospho%20dist%20solid.csv?dl=1")

        # import markers
        df_m <- "phospho_solid_markers"
        }

    }else if(input_type=="protein"){

      df_val <- DRUMLR:::DRUMLprot[,info$cell_names]

      if(model_type=="aml"){

        # File were model names and paths are located:
        file.models <- read.csv(
        "https://www.dropbox.com/s/sircuo6tveuox0g/Paths%20to%20models%20from%20proteomics%20dist%20aml.csv?dl=1")
        # import markers
        df_m <- "prot_aml_markers"

        }else if(solid.or.aml=="solid"){

          # File were model names and paths are located:
          file.models <- read.csv("https://www.dropbox.com/s/hh96xwxyunsfx8h/Paths%20to%20models%20from%20proteomics%20dist%20solid.csv?dl=1")
          # import markers
          df_m <- "prot_solid_markers"
        }

      }else if(input_type=="RNA"){
        df_val <- DRUMLR:::DRUMLrna[,info$cell_names]
        if (model_type=="aml"){
          # File were model names and paths are located:
          file.models <- read.csv("https://www.dropbox.com/s/qao4dvr2ny5d6q0/Paths%20to%20models%20from%20rna%20dist%20aml.csv?dl=1")
          # import markers
          df_m <- "rna_aml_markers"
          }else if(solid.or.aml=="solid"){
            # File were model names and paths are located:
            file.models <- read.csv("https://www.dropbox.com/s/0vb5xco6vwt90v1/Paths%20to%20models%20from%20rna%20dist%20solid.csv?dl=1")
            # import markers
            df_m <- "rna_solid_markers"
          }
        }

  ############################################################################################
  #test compatibility


  #remove missing Variables from dataset
  df_val <- df_val[rownames(df_val) %in%rownames(df_input),]

  #get distances with missing variables
  df_val <- DrugMarkerEnrichment(df=df_val, scale = T, marker_database = df_m)[[1]] %>% t()
  df_val <- df_val[[1]]
  #internal prediction function for caret models
  .PredictModelCaret <- function(file.models,model.name,drug, df_input){
    rmse.of.model.per.drug <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"val.rmse"][1]

    if(length(rmse.of.model.per.drug)>0|is.na(rmse.of.model.per.drug)==F){
      myfile <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"paths"][1]

      if(file.exists(myfile)){
        mymodel <- readRDS(myfile)
        tryCatch({out <- predict(mymodel,df_input)},error=function(e){})}
    }else{out<- data.frame(rep(NA, each = nrow(df_input)))}

    colnames(out)<- model.name
    rownames(out) <- rownames(df_input)
    return(out)
  }

  #internal prediction function for h2o models
  .PredictModelh2o <- function(file.models,model.name,drug, df_input){
    rmse.of.model.per.drug <- file.models[file.models$drug== drug & tolower(file.models$model)==model.name ,"val.rmse"][1]

    if(length(rmse.of.model.per.drug)>0|!is.na(rmse.of.model.per.drug)){
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


    out<- data.frame("dl" = dl,
                     "nnet" = nnet,
                     "pcr" = pcr,
                     "pls" = pls,
                     "rf" = rf,
                     "svm" = svm)

    #caluclate new error values for model dur to missing variables
    out <- apply(out, MARGIN = 2, FUN = function(x){
      rmse <- caret::RMSE(pred = x, obs = df_response, na.rm = T)
      caret::R2(pred = x, obs = df_response, na.rm = T)
      return(rmse)
      })

  }

  return(df_pred)

}
