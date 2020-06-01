#roxygen description
#' @name  Drug.Marker.Enrichment
#' @export Drug.Marker.Enrichment
#' @title  Drug marker Enrichment
#' @usage Drug.Marker.Enrichment(df.ppindex = dat, marker_db = df.aml)
#' @param df.ppIndex df.ppindex is phosphoproteomics, proteomics or transcriptomics data with rownames being identifiers
#' @param marker_db marker database to use to get enrichment

Drug.Marker.Enrichment <- function(df.ppIndex, marker_db){


  # marker_db, database of sensitivity markers
  # df.ppindex = ppindex or protquant data

  df.db <- marker_db
  nc <- ncol(df.ppIndex)
  nr <- nrow(df.db)

  results.pvalues <-data.frame(matrix(1,nrow = nr, ncol = nc))
  rownames(results.pvalues) <- make.names(df.db$drug, unique = T)
  colnames(results.pvalues) <- colnames(df.ppIndex)#[2:nc]

  results.ks.statistic <-data.frame(matrix(nrow = nr, ncol = nc))
  rownames(results.ks.statistic) <- make.names(df.db$drug,unique = T)
  colnames(results.ks.statistic) <- colnames(df.ppIndex)#[2:nc]

  results.median.up <-data.frame(matrix(nrow = nr, ncol = nc))
  rownames(results.median.up) <- make.names(df.db$drug,unique = T)
  colnames(results.median.up) <- colnames(df.ppIndex)#[2:nc]

  results.median.do <-data.frame(matrix(nrow = nr, ncol = nc))
  rownames(results.median.do) <- make.names(df.db$drug,unique = T)
  colnames(results.median.do) <- colnames(df.ppIndex)#[2:nc]

  results.distance <-data.frame(matrix(nrow = nr, ncol = nc))
  rownames(results.distance) <- make.names(df.db$drug,unique = T)
  colnames(results.distance) <- colnames(df.ppIndex)#[2:nc]

  results.area <-data.frame(matrix(nrow = nr, ncol = nc))
  rownames(results.area) <- make.names(df.db$drug,unique = T)
  colnames(results.area) <- colnames(df.ppIndex)#[2:nc]

  results.q <-data.frame(matrix(nrow = nr, ncol = nc))
  rownames(results.q) <- make.names(df.db$drug, unique=T)
  colnames(results.q) <- colnames(df.ppIndex)#[2:nc]

  results.m <-data.frame(matrix(nrow = nr, ncol = 2))
  rownames(results.m) <- make.names(df.db$drug,unique = T)
  colnames(results.m) <- c("up in sensitive","up in resistant")

  results.sites <-data.frame(matrix(nrow = nr, ncol = 2))
  rownames(results.sites) <- make.names(df.db$drug,unique = T)

  all.sites <- data.frame(matrix(nrow=1,ncol=2))
  colnames(all.sites) <- c("site","drug")

  r=1
  for (r in 1:nr) {
    start.time <- Sys.time()
    m.s <- df.db$n.sites.sensitive[r]
    m.r <- df.db$n.sites.resistant[r]
    drug <- df.db$drug[r]
    print(drug)
    if (is.na(m.s) == F & is.na(m.r) == F){
      if(m.s>4 & m.r>4){
        markers.s <- as.character(df.db$up.in.sensitive[r])
        markers.r <- as.character(df.db$up.in.resistant[r])

        ss <- c(unlist(strsplit(markers.s,"-")))
        rr <- c(unlist(strsplit(markers.r,"-")))

        #df.xx <- subset(df.fold,df.fold[,1] %in% paste(ss,";",sep=""))
        # df.xxp<- subset(df.pvalue,df.pvalue[,1] %in% paste(ss,";",sep=""))
        # df.xxf<- subset(df.fdr,df.fdr[,1] %in% paste(ss,";",sep=""))

        df.xx.s <- df.ppIndex[rownames(df.ppIndex) %in% ss,]
        df.xx.r <- df.ppIndex[rownames(df.ppIndex) %in% rr,]
        m.up <- 0
        m.do <- 0
        m.up <- nrow(df.xx.s)
        m.do <- nrow(df.xx.r)

        #sites <-paste(unlist(df.xx[,1]),collapse="-")


        if (length(m.up)>0 & length(m.do)>0){
          if (m.up>2 & m.do>2){
            results.m[r,1] <- m.up
            results.m[r,2] <- m.do
            end.time <- Sys.time()

            c=1
            for (c in 1:nc){
              values.all <- as.numeric(subset(df.ppIndex[,c], df.ppIndex[,c]!=0))
              myvalues.s <- as.numeric(subset(df.xx.s[,c],df.xx.s[,c]!=0))
              myvalues.r <- as.numeric(subset(df.xx.r[,c],df.xx.r[,c]!=0))
              pval <- 1
              ks.e <- 0
              tryCatch({
                myks <- ks.test(myvalues.s,myvalues.r)
                pval <- myks$p.value
                ks.e <- myks$statistic}
                , error=function(e){}
              )
              mysd <- sd(values.all)
              mymedian.r <- median(myvalues.r, na.rm = T)
              mymedian.s <- median(myvalues.s, na.rm = T)
              mymedian.all <- median(values.all, na.rm = T)
              cc <- mymedian.s
              b <- mymedian.r
              f <-   as.numeric(quantile(myvalues.r, na.rm = T)[2])
              e <-   as.numeric(quantile(myvalues.s, na.rm = T)[2])
              b <-   as.numeric(quantile(myvalues.r, na.rm = T)[4])
              a <-   as.numeric(quantile(myvalues.s, na.rm = T)[4])
              a1 <- ((a-b+cc-b)/2)*.25
              a2 <- ((f-e+cc-b)/2)*.25
              myarea <- a1+a2
              results.pvalues[r,c] <- pval
              results.ks.statistic[r,c] <- ks.e
              results.area[r,c] <- myarea
              results.median.do[r,c] <- mymedian.r
              results.median.up[r,c] <- mymedian.s
              results.distance[r,c] <- (mymedian.s - mymedian.r)+(a-b) # the sum of distances at median and quantiles
            }
          }

          print (end.time-start.time)
        }
      }
    }
  }



  return(list(results.distance=results.distance,
              results.area=results.area,
              results.ks.statistic=results.ks.statistic,
              results.pvalues=results.pvalues,
              results.m=results.m,
              results.median.up=results.median.up,
              results.median.do=results.median.do,
              results.sites=results.sites))

}

