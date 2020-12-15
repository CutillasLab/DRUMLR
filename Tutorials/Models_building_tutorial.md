    knitr::opts_chunk$set(echo = TRUE)
    #load libraries
    library(dplyr)
    library(DRUMLR)

Import and scale data
---------------------

Any dataset can be used for DRUML. Inputs need to be dataframes with
variables as rownames.

    #import datasets
    df_input <- DRUMLRdata::DRUMLphos %>% scale(center = T, scale = T)
    df_response <- DRUMLRdata::DRUMLaac

    #get list of AML cell lines
    aml_cells <-
      DRUMLRdata::DRUMLcellinfo[DRUMLRdata::DRUMLcellinfo$unique.tissue == "AML", "Cell.name"]
    input_cells <- DRUMLR::RemoveRepeatNo(colnames(df_input))

    #filter input by cell line
    df_input_aml <- df_input[, input_cells %in% aml_cells]

    #filter drugs by biological range
    drugs <-
      DRUMLR::FilterSensitivity(
        sensitivity_database = df_response,
        iqr_tolerance = 0.1,
        n_coverage = 0.8,
        cell_names = aml_cells
      )

Generate markers
----------------

use the processed data to build markers for DRUMLR model building

    markers <-
      DRUMLR::BuildEMDR(
        df_input = df_input_aml,
        df_response = t(df_response),
        drugs = drugs,
        computational_load = 0.8
      )

\#\#Build Models function .BuildModels contains the paramaters used for
the models

    Model_info <-
      DRUMLR::BuildDRUMLs(
        df_input = df_input_aml,
        .marker_database = markers,
        input_type = "aml_phospho",
        save_path = "~/DRUML_models",
        save_csv = T,
        models = "all",
        df_response = df_response,
        drugs = drugs,
        computational_load = 0.8
      )

Predict\_validation data
------------------------

    df_validation <- external_data

    df_validatiaon_dist <-
      DRUMLR::DrugMarkerEnrichment(df = df_validation, marker_database = markers)

    df_predictions <-
      DRUMLR::PredictDRUML(df_distance = df_validatiaon_dist,
                           models = "all",
                           models_dir = "~/DRUML_models")
