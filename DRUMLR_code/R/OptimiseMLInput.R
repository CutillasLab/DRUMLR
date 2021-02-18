#' Title
#' @title Optimise ML input distance markers
#' @name OptimiseMLInput
#' @param df_train training dataset for machine learning. df_input with response variables in a column named response
#' @param method statistical method used to determine correlation
#' @param max_n_D maximum number of markers to use
#' @param min_n_D Minimum number of  markers to use
#' @param return_corr_analysis TRUE results of correlation analysis will be returned in addition to markers lists
#' @param order_var variable used to order which markers to use
#' @export OptimiseMLInput

OptimiseMLInput <- function(df_train,
                            method = "spearman",
                            max_n_D = 30,
                            min_n_D = 7,
                            partition.ratio = 0.7,
                            return_corr_analysis = FALSE,
                            order_var = "r") {
  train.set <-
    createDataPartition(
      na.omit(df_train$response),
      p = partition.ratio,
      list = FALSE,
      times = 1
    )


  #get training input and response datasets

  df_trr <- df_train[train.set, "response"]
  df_tri <- df_train[train.set, !colnames(df_train) %in% "response"]

  #correlation analysis of distances against response values

  df_corr <- lapply(
    colnames(df_tri),
    FUN = function(x) {
      #get x and y values for analysis
      marker <- x
      x <- df_tri[, x]
      y <- df_trr
      #analyse correlation
      d_corr <- cor.test(x = x, y = y, method = method)
      #draw r and p values from analysis
      p <- d_corr$p.value
      r <- d_corr$estimate
      #calculate slope
      if (!is.na(r)) {
        slope <- coef(lm(formula = y ~ x))[2]
      } else{
        slope <- 0
      }

      #generate output data
      output <- data.frame(
        "marker" = marker,
        "p" = p,
        "r" = r,
        "slope" = slope
      )
      return(output)
    }
  ) %>% dplyr::bind_rows() %>% as.data.frame()


  r_positive <- df_corr[order(-df_corr[, order_var]), ]
  r_negative <- df_corr[order(df_corr[, order_var]), ]

  top_r_positive <- r_positive[r_negative$p <= 0.05 & r_negative[, order_var]>0, ]
  top_r_negative <- r_negative[r_negative$p <= 0.05 & r_negative[, order_var]<0, ]

  #function to determine how many markers to use
  .n_check <- function(x){
    x <- nrow(x)
    if(x>max_n_D){
      out <- max_n_D
    }else if(x<min_n_D){
      out <- min_n_D
    }else{
        out <- x
      }
  }

  #get n values for each marker
  top_pos_n <- .n_check(top_r_positive)
  top_neg_n <- .n_check(top_r_negative)

  #Get D values for ML learning
  top_r_positive <-  r_positive[1:top_pos_n, "marker"] %>% paste()
  top_r_negative <- r_negative[1:top_neg_n, "marker"] %>% paste()

  output <- list(top_r_positive, top_r_negative)
  names(output) <-
    c("sensitive drug markers", "resistant drug markers")

  if (return_corr_analysis) {
    output <- c(output, df_corr)
  }
  return(output)
}
