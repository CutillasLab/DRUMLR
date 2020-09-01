# roxygen description
#' @title Update Internal Databases for KSEA
#' @name Update_databases
#' @usage  Update_databases()
#' @description Update databases from databases files folder



Update_Databases <- function(Inputdir="Input_data"){

  require("dplyr")
  require("readxl")
  #Input Databases

  .LoadDatabase<- function(database){

    path <- paste(Inputdir,"/",database ,sep ="/")

    file_type <- substr(database, start = nchar(database)-3, stop = nchar(database))

    if(file_type == "xlsx"){
      db <- readxl::read_excel(path=path, sheet="ppIndex")%>% data.frame(row.names = 1, stringsAsFactors = F)
    }else if(file_type==".csv"){
      db <- read.csv(file = path, stringsAsFactors = F) %>% data.frame(row.names = 1, stringsAsFactors = F)
    }
    return(db)
  }

  files <- list.files(Inputdir)
  DRUMLphos <- .LoadDatabase(files[2])
  DRUMLprot <- .LoadDatabase(files[3])
  DRUMLprot <- DRUMLprot[, !colnames(DRUMLprot)%in%c("N.A",	"Name",	"Accessions",	"Mascot.Score",	"No.Unique.Peptides",	"No.Pept.Identifications")]
  DRUMLaac <- .LoadDatabase(files[4]) %>% t() %>% data.frame(stringsAsFactors = F)
  DRUMLdruginfo <- .LoadDatabase(files[5])
  DRUMLcellinfo <- .LoadDatabase(files[6])
  #Marker Databases

  .ReadMarker <-function(sheet){
    path = paste(Inputdir, files[7], sep = "/")
    out <- readxl::read_excel(path = path, sheet = sheet) %>% data.frame(stringsAsFactors = F, row.names = 1)
    return(out)
  }

  phospho_aml_markers <- .ReadMarker("phospho aml")
  phospho_solid_markers <- .ReadMarker("phospho solid")
  prot_aml_markers <- .ReadMarker("prot aml")
  prot_solid_markers <- .ReadMarker("prot solid")
  rna_aml_markers <- .ReadMarker("rna aml")
  rna_solid_markers <- .ReadMarker("rna solid")

  usethis::use_data(phospho_aml_markers,
                    phospho_solid_markers,
                    prot_aml_markers,
                    prot_solid_markers,
                    rna_aml_markers,
                    rna_solid_markers,
                    DRUMLaac,
                    DRUMLphos,
                    DRUMLprot,
                    DRUMLcellinfo,
                    DRUMLdruginfo,
                    internal = T,
                    overwrite = T,
                    compress = "xz")

  print("all databases updated")
}

