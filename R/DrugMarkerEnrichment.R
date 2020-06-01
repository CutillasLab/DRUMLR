#Building models with different markers
# roxygen description
#' @title A score of marker ratios
#' @name DrugMarkerEnrichment
#' @export DrugMarkerEnrichment
#' @usage  MarkerZscores(df, marker_database = "aml.markers", scale = T)
#' @description Get zscores between sensitive and resistant markers
#' @param df dataframe of peptide intensities cant be used for DMSO controls yet
#' @param marker_database markers database to be used
#' @param input datatype input can be phospho or prot
#' @param scale if TRUE data will be scaled before putting in ksea

DrugMarkerEnrichment <- function(df, scale= T, marker_database){
  require(dplyr)
  print(paste("running KSEA enrichment of ", marker_database))

  ############################################## adapt ksear to just get mean median and sd

  KsearMMSM <-function(df, database, scale=T){
    #identify database
    datadir <- "D:/OneDrive - Queen Mary, University of London/r_packages/KSEAR-v2_annotated/databases/"
    path <- read.csv(paste(datadir,"database_paths.csv",sep =""), row.names = 1, stringsAsFactors = F)[database, "filename"] %>% as.character()
    path <-paste(datadir,path,sep ="")
    print(paste(gsub(datadir,"",path), "imported"))

    if("sites" %in% colnames(df)==F){dat_edges <- paste(rownames(df))
    }else{dat_edges <- paste(df$sites)
    df <- df[,colnames(df)!="sites"]}

    if(scale == T){df <-  t(df) %>% scale(center = T, scale = T) %>% t() %>% data.frame()}

    #make list of edges
    #if(database %in% phos_dbs){dat_edges<- lapply(dat_edges, FUN = function(x){substr(x, 1, (nchar(x)-1))})%>%unlist()}
    #build dictonary for converting peptides to rownames
    dat_edges_dict<- lapply(dat_edges, function(x){out <- strsplit(x, ";")%>% unlist()
    names(out) <- rep(x,length(out))
    return(out)})
    dat_edges_dict <- data.frame(database_edge = unlist(dat_edges_dict), df_edge = names(unlist(dat_edges_dict)))
    dat_edges_list <- strsplit(dat_edges,";") %>% unlist()
    rownames(df) <- dat_edges

    #format database
    db <- read.csv(path, stringsAsFactors = F)
    rownames(db) <- make.names(db[,1], unique = T)
    db <- db[,-1]

    #make nodes list
    nodes <- paste(rownames(db))

    #get all edges that are present in the db
    nodes_edges <- lapply(nodes, FUN = function(x){
      edges <- strsplit(db[x,2], ";") %>% unlist()
      edges <- edges[edges %in% dat_edges_list]
      return(edges)
    })
    print("nodes indentified")

    #rename nodes_edges name
    names(nodes_edges) <- nodes
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
    #get z

    #Apply functions
    .Getmed <-function(xx){
      d <- df[nodes_edges[[xx]],]
      d <- apply(d, MARGIN = 2, function(x){median(x,na.rm = T)})
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
    dfsd <- lapply(nodes, .Getsd) %>% reformat()

    out <- list(
      "median" = dfmed,
      "mean" = dfmean,
      "sd" = dfsd,
      "m" = data.frame(nodes = nodes, m = nodes_m)
    )
    names(out) <- c(paste(database,"median", sep = "_"),
                    paste(database,"mean", sep = "_"),
                    paste(database,"sd", sep = "_"),
                    paste(database,"m", sep = "_"))

    time<- round(Sys.time()-start, 0)
    print(paste("finished in:", time))
    return(out)
  }

  ############################################################################################


  ksearl <- KsearMMSM(df = df, database = marker_database, scale = scale)
  print("KSEA enrichment finished")
  med <- ksearl[[paste(marker_database,"mean", sep = "_")]]
  sd <-ksearl[[paste(marker_database,"sd", sep = "_")]]
  m <- ksearl[[paste(marker_database,"m", sep = "_")]]
  resistant_med <- med[grepl("_resistant", rownames(med)),]
  sensitivie_med <- med[grepl("_sensitive", rownames(med)),]
  resistant_sd <-  sd[grepl("_resistant", rownames(sd)),]
  sensitive_sd <-  sd[grepl("_sensitive", rownames(sd)),]
  m_res <- m[grepl("_resistant", rownames(m)),]
  m_sens <- m[grepl("_sensitive", rownames(m)),]

  z <- (sensitivie_med-resistant_med)/(sqrt((((resistant_sd^2)/m_res$m)+((sensitive_sd^2)/m_sens$m))))

  rownames(z) <- gsub("_sensitive", "", rownames(z))

  return(z)
}
