# Make Markers example

library(doParallel)
library(foreach)
library(dplyr)
library(limma)

tissue_types<-list(solid = c("liver", "oesophagus"),
                        aml = "haematopoietic_and_lymphoid_tissue",
                        all= c("liver", "oesophagus", "haematopoietic_and_lymphoid_tissue"))

variables <- expand.grid(c("aml", "all", "solid"), c("prot", "phos", "rna")) %>% data.frame()
colnames(variables) <- c("tissues", "input")

for(row_no in rownames(variables)){

  input <- variables[row_no, "input"]
  tissue <- variables[row_no, "tissues"]

  if(input == "prot"){df_input <- DRUMLR:::PDBprot}
  if(input == "phos"){df_input <- DRUMLR:::PDBphos}
  if(input == "rna"){df_input <- DRUMLR:::PDBrna}

  input_info <- DRUMLR::FilterSensitivity(iqr_tolerance = 0.15,
                                          cell_annotations = DRUMLR:::PDBcellinfo,
                                          cell_lines = colnames(df_input),
                                          sensitivity_database = DRUMLR:::PDBaac,
                                          n_coverage = 0.8,
                                          tissue_type = tissue_types[[tissue]]
                                          )

  doParallel::registerDoParallel(cores= 13)

  markers <- foreach::foreach(i = input_info$drugs, .inorder = T, .combine = "rbind",.packages = c("limma", "dplyr"))%dopar%{
    EMSR <-  GenerateEMSR(df_input = df_input[,input_info$cell_names],
                 df_response =DRUMLR:::PDBaac[input_info$cell_names,],
                 drug = i,
                 resampling.times = 10,
                 p.cutoff = 0.05,
                 fold.cut.off = 0.8,
                 marker.n.cutoff = 0.3,
                 direction_percentage = 0.8,
                 computational_load = NULL,
                 pfilt = T)
    return(EMSR$markers_pfilt %>% na.omit())
    }

  doParallel::stopImplicitCluster()
  markers <- data.frame(markers)
  write.csv(x = markers, paste(tissue, input, "markers.csv", sep = "_"))
}
