# Emperical Markers of Sensitivity and resistantce
#' @name GenerateEMSR
#' @export GenerateEMSR
#' @title Emperical Markers of Sensitivity and resistance
#' @param df_input input data with cell lines as colnames and varaibles as rownames defaults to our dataset
#' @param df_response response values dataframe with colnames as cell lines and drugs as column. defaults to our dataset
#' @param drug drug generate markers for
#' @param computational_load decimal fraction of cores you want to use. If left as NULL only 1 core will be used.

GenerateEMSR <- function(df_input= NULL,
                           df_response = NULL,
                           drug,
                           p.cutoff= 0.05,
                           fold.cut.off=0.8,
                           resampling.times= 5,
                           return_limma = F,
                           scale_input = F,
                           computational_load = NULL){

    #Ensure relevent libraries are in the environment
    if("dplyr" %in% (.packages())==FALSE){library(dplyr)}
    if("limma" %in% (.packages())==FALSE){library(dplyr)}
    if("foreach" %in% (.packages())==FALSE){library(foreach)}
    if("doParallel" %in% (.packages())==FALSE){library(doParallel)}

    #scale input if required
    if(scale_input == T){
      cols <- colnames(df_input)
      df_input <- apply(df_input, 1, FUN = function(x){scale(x, center = T, scale = T)}) %>% t() %>% data.frame()
      colnames(df_input) <- cols
    }

    #Verify ronames are variables
    if(is.character(rownames(df_input))== F&is.factor(rownames(df_input))== F){
      print("rownames are not variable names")
      errorCondition(message = "rownames are not variable names")
      stop()}

    #Initiate multicore processing if specified
    if(is.null(computational_load)==F){
      cores <- round((computational_load * parallel::detectCores()), 0)
      registerDoParallel(cores = cores)
      print(paste("running on", cores, "cores"))
    }

    #NULL input variables
    if(is.null(df_input)){df_input <- data.frame(DRUMLR:::DRUMLphos, row.names = 1)}
    if(is.null(df_response)){df_response <- DRUMLR:::DRUMLaac}

    #convert column names into cell lines
    col_cell_lines <- RemoveRepeatNo(colnames(df_input))

    #Split Cell lines into resistant and sensitive groups
    df_response <- df_response[rownames(df_response)%in% col_cell_lines, ]
    df_response <- df_response[,drug, drop = F] %>% na.omit()
    median_filter <- df_response[,drug] >= median(df_response[,drug], na.rm = T)

    #make sensitive data.frame
    sensitive_cells <- rownames(df_response)[median_filter] %>% toupper()
    df_sens <- df_input[,col_cell_lines %in% sensitive_cells]

    #make resistant data.frame
    resistant_cells <- rownames(df_response)[!median_filter] %>% toupper()
    df_res <- df_input[,col_cell_lines %in% resistant_cells]

    print("Cells split into sensitive and resistant groups")

    # split sensitive and resistant cells into k*rs groups
    fm.up <- caret::createMultiFolds(colnames(df_sens),k=2,times = resampling.times)
    fm.do <- caret::createMultiFolds(colnames(df_res),k=2,times = resampling.times)
    grid <- expand.grid(1:length(fm.up), 1:length(fm.up))
    sumfolds <- nrow(grid)

    print(paste(Sys.time(), ":Initiating", drug, "EMSR analysis. Analysing", sumfolds, "folds"))
    startt <- Sys.time()

    if(is.null(computational_load)){
      EMSR_analysis <-  foreach::foreach(i = rownames(grid),
                                         .inorder = T,
                                         .packages = c("limma", "dplyr"),
                                         .combine = "cbind")%do%{
                                           #get sensitivity and resistance inputs from folds
                                           sensitive_cells <- colnames(df_sens)[fm.up[grid[i, "Var1"]][[1]]]
                                           resistant_cells <- colnames(df_res)[fm.do[grid[i, "Var2"]][[1]]]

                                           df1 <- df_sens[,sensitive_cells, drop= F]
                                           df2 <- df_res[,resistant_cells, drop= F]
                                           sens_cell <- RemoveRepeatNo(sensitive_cells)
                                           res_cell <- RemoveRepeatNo(resistant_cells)

                                           if(ncol(df1) <=3|ncol(df2)<=3){
                                             print(paste(i, "not enough n values"))

                                             out <- data.frame("folds" = rep(0,each = nrow(df1)),
                                                               "pvalue" = rep(1,each = nrow(df2)))

                                             rownames(out)<- rownames(df1)

                                             colnames(out) <- paste(i,c("folds","pvalue"), sep = ".")

                                           }else if(t.test(df_response[sens_cell,],df_response[res_cell,])$p.value >0.05){
                                             out <- data.frame("folds" = rep(0,each = nrow(df1)),
                                                               "pvalue" = rep(1,each = nrow(df2)))
                                             print(paste(i, "no significant difference between groups"))

                                             rownames(out)<- rownames(df1)

                                             colnames(out) <- paste(i,c("folds","pvalue"), sep = ".")

                                           }else{

                                             sens_cellname_fold<-colnames(df1)
                                             res_cellname_fold<-colnames(df2)

                                             dfall <- cbind(df1, df2)

                                             sensitivity <- c(2,1)[as.factor(unclass(colnames(dfall) %in% sens_cellname_fold))]
                                             #df1 <- data.frame(sensitivity = as.factor(1), df1)
                                             #df2 <- data.frame(sensitivity = as.factor(2), df2)

                                             sens_cellname_fold<-sens_cellname_fold %>% paste(collapse = ";")
                                             res_cellname_fold<-res_cellname_fold %>% paste(collapse = ";")

                                             design <- stats::model.matrix(~ 0+factor(c(sensitivity)))
                                             colnames(design) <- c("sensitive","resistant")
                                             contrast.matrix <- limma::makeContrasts(sensitive-resistant,
                                                                                     levels=design)

                                             fit <- limma::lmFit(dfall,design)
                                             fit2 <- limma::contrasts.fit(fit, contrast.matrix)
                                             fit2 <- limma::eBayes(fit2)
                                             pvals <- data.frame(fit2$p.value)
                                             fvals <- data.frame(fit2$coefficients)

                                             out <- cbind(pvals, fvals)
                                             colnames(out)<- paste(i,c("pvalue","folds"), sep = ".")
                                           }
                                           print(paste(i, "completed"))
                                           return(out)
                                         }

    }else{
      EMSR_analysis <-  foreach::foreach(i = rownames(grid),
                                         .inorder = T,
                                         .packages = c("limma", "dplyr"),
                                         .combine = "cbind")%dopar%{

                                           RemoveRepeatNo <- function(x, repeatseperator ="__"){
                                             out <- lapply(x, FUN = function(x){strsplit(x, repeatseperator)[[1]][[1]]})%>% unlist()
                                             return(out)
                                           }


                                           #get sensitivity and resistance inputs from folds
                                           sensitive_cells <- colnames(df_sens)[fm.up[grid[i, "Var1"]][[1]]]
                                           resistant_cells <- colnames(df_res)[fm.do[grid[i, "Var2"]][[1]]]

                                           df1 <- df_sens[,sensitive_cells, drop = F]
                                           df2 <- df_res[,resistant_cells, drop = F]

                                           sens_cell <- RemoveRepeatNo(sensitive_cells)
                                           res_cell <- RemoveRepeatNo(resistant_cells)

                                           if(ncol(df1) <=3|ncol(df2)<=3){
                                             print(paste(i, "not enough n values"))

                                             out <- data.frame("folds" = rep(0,each = nrow(df1)),
                                                               "pvalue" = rep(1,each = nrow(df2)))

                                             rownames(out)<- rownames(df1)

                                             colnames(out) <- paste(i,c("folds","pvalue"), sep = ".")

                                           }else if(t.test(df_response[sens_cell,],df_response[res_cell,])$p.value >0.05){

                                             out <- data.frame("folds" = rep(0,each = nrow(df1)),
                                                               "pvalue" = rep(1,each = nrow(df2)))

                                             rownames(out)<- rownames(df1)
                                             colnames(out) <- paste(i,c("folds","pvalue"), sep = ".")

                                           }else{

                                             sens_cellname_fold<-colnames(df1)
                                             res_cellname_fold<-colnames(df2)

                                             dfall <- cbind(df1, df2)

                                             sensitivity <- c(2,1)[as.factor(unclass(colnames(dfall) %in% sens_cellname_fold))]
                                             #df1 <- data.frame(sensitivity = as.factor(1), df1)
                                             #df2 <- data.frame(sensitivity = as.factor(2), df2)

                                             sens_cellname_fold<-sens_cellname_fold %>% paste(collapse = ";")
                                             res_cellname_fold<-res_cellname_fold %>% paste(collapse = ";")

                                             design <- stats::model.matrix(~ 0+factor(c(sensitivity)))
                                             colnames(design) <- c("sensitive","resistant")
                                             contrast.matrix <- limma::makeContrasts(sensitive-resistant,
                                                                                     levels=design)

                                             fit <- limma::lmFit(dfall,design)
                                             fit2 <- limma::contrasts.fit(fit, contrast.matrix)
                                             fit2 <- limma::eBayes(fit2)
                                             pvals <- data.frame(fit2$p.value)
                                             fvals <- data.frame(fit2$coefficients)

                                             out <- cbind(pvals, fvals)
                                             colnames(out)<- paste(i,c("pvalue","folds"), sep = ".")
                                           }

                                           return(out)
                                         }

    }

    print(paste(Sys.time(), ":Processing EMSR results"))
    rownames(EMSR_analysis)<- rownames(df_sens)
    #split data into fold and pvalue data
    EMSR_analysis <- list("df_folds" = data.frame(EMSR_analysis[,grepl(".folds", colnames(EMSR_analysis), fixed = T), drop = F]),
                          "df_pvalues" = data.frame(EMSR_analysis[,grepl(".pvalue", colnames(EMSR_analysis), fixed = T), drop = F]))


    #get information for cell lines comprising each iteration and iteration sequence
    fold_groups <- list("sensitive" = lapply(fm.up, FUN = function(x){colnames(df_sens)[x]%>% paste(collapse = ";")}) %>% unlist() %>% data.frame(),
                        "resistant" = lapply(fm.do, FUN = function(x){colnames(df_res)[x]%>% paste(collapse = ";")}) %>% unlist() %>%data.frame())
    fold_groups <- merge.data.frame(fold_groups$sensitive, fold_groups$resistant, by = "row.names", all = T)
    fold_groups <- data.frame(rownames(fold_groups), fold_groups)
    colnames(fold_groups) <- c("fold_group_number", "repeat", "sensitive", "resistant")

    grid <- data.frame(rownames(grid), data.frame(grid))
    colnames(grid)<- c("repeat", "fold_group_number_sensitive", "fold_group_number_resistant")


    #count number of iterations which satisfied pvalue and fold thresholds
    ###############################################################################################################
    #binary analysis of folds and pvalue
    df_folds <- EMSR_analysis$df_folds
    #get average fold values
    mean_fold <- apply(df_folds, 1, FUN = function(x){mean(x,na.rm = T)})

    #get fold up and down regulations
    df_folds[df_folds>=fold.cut.off] <- "u"
    df_folds[df_folds<=-(fold.cut.off)] <- "d"
    n_increase <- apply(df_folds,1, function(x) sum(x=="u"))
    n_decrease <- apply(df_folds,1, function(x) sum(x=="d"))

    #get average p value
    mean_p <- apply(EMSR_analysis$df_pvalues, 1, FUN = function(x){mean(x,na.rm = T)})

    #Filter by p values
    df_folds <- EMSR_analysis$df_folds
    df_folds[EMSR_analysis$df_pvalues > p.cutoff]<- 0

    #get average fold values
    mean_fold_filt <- apply(df_folds, 1, FUN = function(x){mean(x,na.rm = T)})

    #get significant fold up and down regulations
    df_folds[df_folds>=fold.cut.off] <- "u"
    df_folds[df_folds<=-(fold.cut.off)] <- "d"
    n_increase_filt <- apply(df_folds,1, function(x) sum(x=="u"))
    n_decrease_filt <- apply(df_folds,1, function(x) sum(x=="d"))



    #get the number of significant values
    nsig <- apply(df_folds, 1, FUN = function(x){length(na.omit(x))})

    #get number of significant up and down regulations
    df_folds <- EMSR_analysis$df_folds

    #label significant pvalues
    out <- data.frame( "Variable" = rownames(df_input),
                       "n_pvalue_cutoff" =  nsig,
                       "n_increase" =  n_increase,
                       "n_decrease" = n_decrease,
                       "ratio" = (n_increase - n_decrease),
                       "n_increase_pfilt" =  n_increase_filt,
                       "n_decrease_pfilt" = n_decrease_filt,
                       "ratio_pfilt" = (n_increase_filt - n_decrease_filt),
                       "mean_fold" = mean_fold,
                       "mean_fold_pfilt" = mean_fold_filt,
                       "mean_p" = mean_p
    )

    if(return_limma ==T){out <- list("summary" = out,
                                     "folds" = EMSR_analysis$df_folds,
                                     "pvalues" = EMSR_analysis$df_pvalues,
                                     "Fold info"= list("fold groups" = fold_groups,
                                                       "folds" = grid))}


    if(is.null(computational_load)==F){
      doParallel::stopImplicitCluster()
      return(out)
    }
  }


