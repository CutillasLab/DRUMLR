# roxygen description
#' @title Update Internal Databases for KSEA
#' @name Update_databases
#' @export Update_databases
#' @usage  Update_databases()
#' @description Update databases from databases files folder
#' @param datadir database directory which leads to database_paths.csv file
#' }

Update_databases <- function(datadir= "D:/OneDrive - Queen Mary, University of London/r_packages/KSEAR-v2_annotated/databases/"){
  require(dplyr)
  require(usethis)
  load_database<- function(database){
    path <- read.csv(paste(datadir,"database_paths.csv",sep =""), row.names = 1, stringsAsFactors = F)[database, "filename"] %>% as.character()
    path <-paste(datadir,path,sep ="")
    db <- read.csv(path, stringsAsFactors = F) %>% data.frame(stringsAsFactors = F)
    return(db)
  }

  ctams <- load_database("ctams")
  edges <- load_database("edges")
  function. <-load_database("function.")
  location <- load_database("location")
  nci <- load_database("nci")
  PDBmodules <- load_database("PDBmodules")
  pdts <- load_database("pdts")
  pdts.function <- load_database("pdts.function")
  pdts.location <- load_database("pdts.location")
  pdts.nci <- load_database("pdts.nci")
  pdts.process <- load_database("pdts.process")
  pdts.reactome <- load_database("pdts.reactome")
  process <- load_database("process")
  psite <- load_database("psite")
  reactome <- load_database("reactome")
  signor <- load_database("signor")
  PDBaac <- henryspackage::LoadData("PDBaac")
  PDBphos <- henryspackage::LoadData("PDBphos")
  PDBprot <- henryspackage::LoadData("PDBprot")
  PDBrna <- henryspackage::LoadData("PDBrna")
  PDBcellinfo <- henryspackage::LoadData("PDBcellinf")[,c("cell.name", "unique.tissue")]
  PDBcellinfo$cell.name <- make.names(PDBcellinfo$cell.name, unique = T)
  PDBaac$cell.name <- make.names(PDBaac$cell.name, unique = T)

  usethis::use_data(ctams,
                    edges,
                    function.,
                    location,
                    nci,
                    PDBmodules,
                    pdts,
                    pdts.function,
                    pdts.location,
                    pdts.nci,
                    pdts.process,
                    pdts.reactome,
                    process,
                    psite,
                    reactome,
                    signor,
                    PDBaac,
                    PDBphos,
                    PDBprot,
                    PDBrna,
                    PDBcellinfo,
                    internal = T, overwrite = T, compress = "xz")

  print("all databases updated")
}

Update_databases()
