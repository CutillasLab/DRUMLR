#GenerateEMSR
# Emperical Markers of Sensitivity and resistantce
#' @name BuildEMSR
#' @export BuildEMSR
#' @title Emperical Markers of Sensitivity and resistance
#' @param df_input input data with cell lines as colnames and varaibles as rownames
#' @param df_response response values dataframe with colnames as cell lines and drugs as column
#' @param drug drug generate markers for
#' @param p.cutoff P value threshold to define a repeat comparison as significant
#' @param fold.cut.off fold value threshold to define a repeat comparison direction
#' @param resampling.times number of resampling times
#' @param return_limma Set as TRUE if for limma analysis results to be returned
#' @param scale_input Set as TRUE to scale input data
#' @param computational_load Set as a decimal fraction of cores you want to use. If left as NULL only 1 core will be used.


#make the cell lines thing required



#test data
#df_input <- readxl::read_excel("C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data/Suppl_Dataset_2_Phosphoproteomics_of_48_cell_lines_in_triplicate.xlsx", sheet = "ppIndex") %>% data.frame(row.names = 1, stringsAsFactors = F)
#df_response <-read.csv("C:/Users/Henry_gerdes/Documents/GitHub/DRUML_Data/Input_data/Suppl_Dataset_4_Drugs_with_normalised_AAC_and_class_info.csv", row.names = 1)
#p.cutoff = 0.05
#fold.cut.off = 0.8
#resampling.times = 5
#computational_load = NULL
#nfolds <-  10
#maxmarker_res <- 100
#maxmarker_sens <- 100

BuildEMSR <- function(df_input = NULL,
                         df_response = NULL,
                         drugs,
                         p.cutoff = 0.05,
                         fold.cut.off = 1,
                         resampling.times = 5,
                         return_limma = F,
                         computational_load = NULL,
                         nfolds = 10,
                         maxmarker_res = 100,
                         maxmarker_sens = 100) {

  #Ensure relevent libraries are in the environment
  if ("dplyr" %in% (.packages()) == FALSE) {
    library(dplyr)
  }
  if ("limma" %in% (.packages()) == FALSE) {
    library(dplyr)
  }
  if ("foreach" %in% (.packages()) == FALSE) {
    library(foreach)
  }
  if ("doParallel" %in% (.packages()) == FALSE) {
    library(doParallel)
  }

  #Verify rownames are variables
  if (is.character(rownames(df_input)) == F &
      is.factor(rownames(df_input)) == F) {
    print("rownames are not variable names")
    errorCondition(message = "rownames are not variable names")
    stop()
  }

  #Initiate multicore processing if specified
  if (is.null(computational_load) == F) {
    cores <- round((computational_load * parallel::detectCores()), 0)
    registerDoParallel(cores = cores)
    print(paste("running on", cores, "cores"))
  }

  #convert column names into cell lines
  col_cell_lines <- DRUMLR::RemoveRepeatNo(colnames(df_input))

  #remove cell lines not present in database
  df_response_ori <-
    df_response[,colnames(df_response) %in% col_cell_lines]
  if (is.null(drugs)){
    drugs <- rownames(df_response)
    }
  ######################################################################################################################################
  #Drug loop will need to begin here

  ######################################################################################################################################
  marker_db <- foreach(i = drugs,
          .inorder = T,
          .combine = "rbind",
          .packages = c("dplyr", "DRUMLR", "caret", "foreach"))%dopar%{
  drug <- i
  #Split Cell lines into resistant and sensitive groups
  df_response <- df_response_ori[drug, ,drop = F] %>% t() %>% na.omit()


  median_filter <-
    df_response[, drug] >= median(df_response[, drug], na.rm = T)

  #make sensitive data.frame
  sensitive_cells <-
    rownames(df_response)[median_filter] %>% toupper()
  df_sens <- df_input[, col_cell_lines %in% sensitive_cells]

  #make resistant data.frame
  resistant_cells <-
    rownames(df_response)[!median_filter] %>% toupper()
  df_res <- df_input[, col_cell_lines %in% resistant_cells]

  print("Cells split into sensitive and resistant groups")

  # split sensitive and resistant cells into k*rs groups
  fm.up <-
    caret::createMultiFolds(colnames(df_sens), k = 2, times = resampling.times)
  fm.do <-
    caret::createMultiFolds(colnames(df_res), k = 2, times = resampling.times)
  grid <- expand.grid(1:length(fm.up), 1:length(fm.do))
  sumfolds <- nrow(grid)

  print(
    paste(
      Sys.time(),
      ":Initiating",
      drug,
      "EMSR analysis. Analysing",
      sumfolds,
      "folds"
    )
  )
  # filter >3

  bfold <- lengths(fm.up) >=3 & lengths(fm.do) >=3
  fm.up <- fm.up[bfold]
  fm.do <- fm.do[bfold]
  #Internal Function for limma analysis

  .LimmaEMSR <- function(i){
    #get sensitivity and resistance inputs from folds
    sensitive_cells <-
      fm.up[grid[i, "Var1"]][[1]]
    resistant_cells <-
      fm.do[grid[i, "Var2"]][[1]]
    df1 <-
      df_sens[, sensitive_cells, drop = F]
    df2 <-
      df_res[, resistant_cells, drop = F]

      dfall <-
        cbind(df1, df2)

      sensitivity <-
        c(2, 1)[as.factor(unclass(colnames(dfall) %in% colnames(df_sens)))]

      design <-
        stats::model.matrix( ~ 0 + factor(c(sensitivity)))
      colnames(design) <-
        c("sensitive", "resistant")
      contrast.matrix <-
        limma::makeContrasts(sensitive - resistant,
                             levels =
                               design)

      fit <-
        limma::lmFit(dfall, design)
      fit2 <-
        limma::contrasts.fit(fit, contrast.matrix)
      fit2 <-
        limma::eBayes(fit2)

      fvals <-
        data.frame(fit2$coefficients)

      print(paste("fold ", i, " analysed"))
      out <- fvals
      colnames(out) <- i

    return(out)
  }


  startt <- Sys.time()


    df_folds <-  foreach::foreach(
      i = rownames(grid),
      .inorder = T,
      .packages = c("limma", "dplyr"),
      .combine = "cbind"
    ) %do% {
      .LimmaEMSR(i)
    }

  print(paste(Sys.time(), ":Processing EMSR results"))
  #split data into fold and pvalue data

  #count number of iterations which satisfied pvalue and fold thresholds
  ###############################################################################################################

  #get average fold values
  mean_fold <-
    apply(
      df_folds,
      1,
      FUN = function(x) {
        mean(x, na.rm = T)
      }
    )

  #get fold up and down regulations
  df_folds[df_folds < fold.cut.off & df_folds > -(fold.cut.off)] <- NA
  df_folds[df_folds >= fold.cut.off] <- 1
  df_folds[df_folds <= -(fold.cut.off)] <- -1

  df_limma <- data.frame(
    "ratio" = apply(df_folds, 1, function(x){sum(x, na.rm = T)}),
    "mean_fold" = mean_fold
  )


  #build dataframe of sensitive markers
  #alter resistant fold values so they increase with significance
  res_m <- df_limma[df_limma$ratio <0,] * -(1)
  sens_m <- df_limma[df_limma$ratio >0,]

  #reorder markers by mean fold
  colnames(res_m)
  res_m <- res_m[order(res_m[, "mean_fold"], res_m[, "ratio"], decreasing = T), ]
  res_m_names <- rownames(res_m)
  sens_m <- sens_m[order(sens_m[, "mean_fold"], sens_m[, "ratio"], decreasing = T), ]
  sens_m_names <- rownames(sens_m)
  #identify testing paramaters for marker lengths
  if (is.null(maxmarker_res)) {
    maxmarker_res <- nrow(res_m)
  }
  if (is.null(maxmarker_sens)) {
    maxmarker_sens <- nrow(sens_m)
  }

  sens_test <-
    seq(
      from = 5,
      to = maxmarker_sens,
      by = round(maxmarker_sens / nfolds, 0)
    ) #%>% list()
  res_test <-
    seq(
      from = 5,
      to = maxmarker_res,
      by = round(maxmarker_res / nfolds, 0)
    ) #%>% list()

  print("generating markers lists")
  #generate list of markers using lengths of sens_test and res_test
  sens_markers <-
    lapply(
      sens_test,
      FUN = function(x) {
        sens_m_names[1:x]
      }
    )
  res_markers <-
    lapply(
      res_test,
      FUN = function(x){
        res_m_names[1:x]

      }
    )

  #build tuning grid for foreachfunction
  grid <- expand.grid(1: length(sens_markers), 1:length(res_markers))

  #Internal function for generating marker distances


  .GetDistances <- function(i) {
    x <- grid[i, "Var1"]
    y <- grid[i, "Var2"]

    sens_m <- sens_markers[[x]]
    res_m <- res_markers[[y]]

    sens <- df_input[sens_m,]
    s_med <- apply(
      sens,
      2,
      FUN = function(x)
        (median(x, na.rm = T))
    )
    s_q3 <-
      apply(
        sens,
        2,
        FUN = function(x)
          (quantile(
            x, probs = 0.75, na.rm = T
          ))
      )

    res <- df_input[res_m,]
    r_med <- apply(
      res,
      2,
      FUN = function(x)
        (median(x, na.rm = T))
    )
    r_q3 <-
      apply(
        res,
        2,
        FUN = function(x)
          (quantile(
            x, probs = 0.75, na.rm = T
          ))
      )

    D <- data.frame((s_med - r_med) + (s_q3 - r_q3))
    colnames(D) <- i
    return(D)
  }

  #Internal function of analysing marker distance correlations with drug sensitivity
  .SpearmanMarkerAnalysis <- function(i) {
    cc <-
      cor.test(df_dist[,i],
               df_dist[,"response"],
               method = "spearman",
               complete.cases = TRUE)
    out <-
      data.frame("pvalue" = cc$p.value,
                 "estimate" = cc$estimate)
    out$rep <- c(1:nrow(out))
    return(out)
  }

  #get distances using marker lists on one core
    df_dist <- foreach(
      i = rownames(grid),
      .inorder = T,
      .combine = "cbind",
      .packages = "dplyr"
    ) %do% {

      out <- .GetDistances(i) %>% data.frame(stringsAsFactors = F)
      print(paste(
        "Distances complete for comparison: ",
        i,
        "/",
        nrow(grid),
        sep = ""
      ))
      return(out)
    }

    print("begining spearman correlation analysis of marker groups")

    #get response values for optimisation
    df_dist[, "response"] <- data.frame(df_response)[col_cell_lines,]

    #analyse correlation between distances and drug sensitivity using spearman ranking on multiple cores
    marker_analysis <- foreach(
      i = 1:(ncol(df_dist)-1),
      .inorder = T,
      .combine = "rbind",
      .packages = c("dplyr")
    ) %do% {
      out <- .SpearmanMarkerAnalysis(i)
      print(paste(
        "Spearman analysis complete for comparison: ",
        i,
        "/",
        ncol(df_dist),
        sep = ""
      ))
      return(out)

    }

  #get numbers of markers used for each comparison
  marker_analysis$rep <- 1:nrow(marker_analysis)
  marker_analysis$n_sens <-
    lapply(
      marker_analysis$rep,
      FUN = function(x) {
        sens_test[grid[x, "Var1"]]
      }
    ) %>% unlist() %>% as.double()
  marker_analysis$n_res <-
    lapply(
      marker_analysis$rep,
      FUN = function(x) {
        res_test[grid[x, "Var2"]]
      }
    ) %>% unlist() %>% as.double()

  #identify optimal marker length combination
  optimal <-
    marker_analysis[marker_analysis$pvalue == min(marker_analysis$pvalue),]

  #if multiple combinations generate optimal marker combinations use the one with the fewest total markers
  if (nrow(optimal) > 1) {
    tot <- optimal$n_sens + optimal$n_res
    optimal <- optimal[tot %in% min(tot),]
  }

  #generate marker output
  print("analysis complete")
  out <- data.frame(
    "drugs" = drug,
    "m_sens" = optimal$n_sens,
    "m_res" = optimal$n_res,
    "senstive_markers" = paste(sens_markers[[grid[optimal$rep, "Var1"]]], collapse = "-"),
    "resistant_markers" = paste(res_markers[[grid[optimal$rep, "Var2"]]], collapse = "-"),
    stringsAsFactors = F
  )
  colnames(out) <- c("drugs", "m_sens", "m_res", "sensitive_markers", "resistant_markers")

  return(out)
}

  if (is.null(computational_load) == F) {
    doParallel::stopImplicitCluster()
  }

  return(marker_db)
}
