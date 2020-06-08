# Make Markers example

library(doParallel)
library(foreach)
library(dplyr)
library(limma)

tissue_types<-list(solid = c( "Esophageal", "Hepatocellular"),
                        aml = "AML",)


  input_info <- FilterSensitivity(iqr_tolerance = 0.05,
                                          cell_annotations = DRUMLR:::DRUMLcellinfo,
                                          cell_lines = colnames(df_input),
                                          sensitivity_database = DRUMLR:::DRUMLaac,
                                          n_coverage = 0.8,
                                          tissue_type = "AML"
                                          )

  doParallel::registerDoParallel(cores= 13)

  #filter input data by cancer type cell names

  df_input <- df_input[, RemoveRepeatNo(colnames(df_input))%in% input_info$cell_names]

  markers <- foreach::foreach(i = input_info$drugs,
                              .inorder = T,
                              .combine = "rbind",
                              .packages = c("limma", "dplyr"))%dopar%{


    EMSR <-  GenerateEMSR(df_input = df_input,
                 df_response =DRUMLR:::DRUMLaac,
                 drug = "barasertib",
                 resampling.times = 10,
                 p.cutoff = 0.05,
                 fold.cut.off = 0.8,
                 computational_load = 0.8
                 )

    out <- OptimiseMarkers(EMSRoutput = EMSR,
                    drug = "barasertib",
                    df_response = DRUMLR:::DRUMLaac,
                    df_input = df_input,
                    scale_input = F,
                    maxmarker_res = 200,
                    maxmarker_sens = 200,
                    nfolds = 50,
                    computational_load = 0.8
                      )

    return(out)
    }
