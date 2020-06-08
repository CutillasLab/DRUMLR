#Building models with different markers
# roxygen description
#' @title Generate D values from input data and marker database
#' @name DrugMarkerEnrichment
#' @export DrugMarkerEnrichment
#' @usage  DrugMarkerEnrichment(df = df.ppindex, marker_database = "aml_prot_markers", scale = T, output = c("Distance", "zscore"))
#' @description Get zscores between sensitive and resistant markers
#' @param df dataframe of peptide intensities cant be used for DMSO controls yet
#' @param marker_database markers database to be used can be aml_phos_markers, all_phos_markers, solid_phos_markers, aml_prot_markers, all_prot_markers, solid_prot_markers, aml_rna_markers, all_rna_markers or solid_rna_markers
#' @param scale if TRUE data will be scaled before putting Drug Enrichment
#' @param output The output metric you want to use can be "Distance" or "zscore"

DrugMarkerEnrichment <- function(df, scale= T, marker_database, output = "Distance"){
  if("dplyr" %in% (.packages())==FALSE){library(dplyr)}
  print(paste("running Drug Marker Enrichment of ", marker_database))
  dbs  <- GetIntMarkerDB(marker_database)

  #Build database dictionary using rownames or sites column
  if("sites" %in% colnames(df)==F){dat_edges <- paste(rownames(df))
      }else{dat_edges <- paste(df$sites)
      df <- df[,colnames(df)!="sites"]}

  #build df dictionary using data input
  dat_edges_dict<- lapply(dat_edges, function(x){out <- strsplit(x, ";")%>% unlist()
    names(out) <- rep(x,length(out))
    return(out)})
  dat_edges_dict <- data.frame(database_edge = unlist(dat_edges_dict), df_edge = names(unlist(dat_edges_dict)))
  dat_edges_list <- strsplit(dat_edges,";") %>% unlist()
  rownames(df) <- dat_edges

  #scale data if needed
  if(scale == T){df <-  t(df) %>% scale(center = T, scale = T) %>% t() %>% data.frame()}

 DMEnrich <-function(df, database, scale=T){

    db <- dbs[[database]]

    #make nodes list
    nodes <- paste(rownames(db))

    #get all edges that are present in the db
    nodes_edges <- lapply(nodes, FUN = function(x){
      edges <- strsplit(paste(db[x,2]), ";") %>% unlist()
      edges <- edges[edges %in% dat_edges_list]
      return(edges)
    })
    names(nodes_edges) <- nodes
    print("nodes indentified")

    #rename nodes_edges name

    nodes_edges <-nodes_edges[lengths(nodes_edges)>2]
    nodes <- names(nodes_edges)
    print("edges compiled")
    #change node edges to their corresponding rownames in the input dataframe

    nodes_edges <- lapply(nodes, function(x){
      matches <-  dat_edges_dict[dat_edges_dict$database_edge %in% nodes_edges[[x]], "df_edge"] %>% paste() %>% unique()
      if(length(matches)<=2|length(matches)==0){out <- NA}
      return(matches)
    })
    names(nodes_edges) <- nodes
    print("edges reformatted")

    # Filter nodes based off of m value and make m list
    nodes_edges <- nodes_edges[lengths(nodes_edges)>2]
    nodes_m <- lengths(nodes_edges)
    nodes <- names(nodes_edges)

    #Apply functions
    .Getmed <-function(xx){
      d <- df[nodes_edges[[xx]],]
      d <- apply(d, MARGIN = 2, function(x){median(x,na.rm = T)})
      return(d)
    }

    .Getq3 <-function(xx){
      d <- df[nodes_edges[[xx]],] %>% na.omit()
      d <- apply(d, MARGIN = 2, function(x){quantile(x = x, probs = 0.75)})
      return(d)
    }

    .Getsd <- function(xx){
      d <- df[nodes_edges[[xx]],]
      d <- apply(d, MARGIN = 2, function(x){sd(x,na.rm = T)})
      return(d)
    }

    .Getmean <-function(xx){
      d <- df[nodes_edges[[xx]],]
      d <- apply(d, MARGIN = 2, function(x){mean(x,na.rm = T)})
      return(d)
    }

    reformat <- function(df, scale = F){
      names(df) <- nodes
      df <- data.frame(df)
      if(scale == F){
        df <-t(df) %>% data.frame()
      }else{
        df <- scale(df, center = T, scale = T) %>% data.frame()
      }
      return(df)
    }


    #build dataframes

    dfmed <- lapply(nodes, .Getmed) %>% reformat()
    print("median values calculated")
    dfmean <- lapply(nodes, .Getmean) %>% reformat()
    print("mean values calculated")
    dfq3 <- lapply(nodes, .Getq3) %>% reformat()
    print("Q3 values calculated")
    dfsd <- lapply(nodes, .Getsd) %>% reformat()
    print("Standard Deviations calculated")

    out <- list(
      "median" = dfmed,
      "mean" = dfmean,
      "sd" = dfsd,
      "q3" = dfq3,
      "m" = data.frame(nodes = nodes, m = nodes_m)
    )
    names(out) <- c(paste(database,"median", sep = "_"),
                    paste(database,"mean", sep = "_"),
                    paste(database,"sd", sep = "_"),
                    paste(database,"q3", sep = "_"),
                    paste(database,"m", sep = "_"))

    print(paste("finished", database, "analysis"))
    return(out)
  }

  ############################################################################################


  sens_enr <- DMEnrich(df = df, database = "sensitive_markers", scale = scale)
  res_enr <- DMEnrich(df = df, database = "resistant_markers", scale = scale)
  print("KSEA enrichment finished")
  node_list <-  merge.data.frame(res_enr$resistant_markers_m,  sens_enr$sensitive_markers_m, by="nodes")$nodes
  resistant_med <- res_enr$resistant_markers_median[node_list,]
  sensitivie_med <- sens_enr$sensitive_markers_median[node_list,]
  resistant_q3 <- res_enr$resistant_markers_q3[node_list,]
  sensitive_q3 <- sens_enr$sensitive_markers_q3[node_list,]
  resistant_sd <-  res_enr$resistant_markers_sd[node_list,]
  sensitive_sd <-  sens_enr$sensitive_markers_sd[node_list,]
  m_res <- res_enr$resistant_markers_m[node_list,]
  m_sens <- sens_enr$sensitive_markers_m[node_list,]

  z <- (sensitivie_med-resistant_med)/(sqrt((((resistant_sd^2)/m_res$m)+((sensitive_sd^2)/m_sens$m))))
  D <- (sensitivie_med -sensitive_q3)-(resistant_med -resistant_q3)

  out <- list("Distance" = D, "zscore" = z)
  out <- out[output]

  return(out)
}
