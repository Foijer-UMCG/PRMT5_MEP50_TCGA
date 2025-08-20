rm(list = ls())
here::i_am("scripts/check_data_integrity.R")

# need to library for %>%
library(dplyr)

source(file.path(here::here(),
                 "scripts/util_funcs.R"))

data_dirs <- list.dirs(file.path(here::here(),
                                 "data"),
                       recursive = FALSE)

# force recalculating aneuploidy, even if file is already present
force = FALSE

for (dir in data_dirs) {
  # this doesn't seem to work for some reason
  if(file.exists(file.path(dir,
                           "filtered_data.Rds")) &
     !force){
    dataset <- basename(dir)
    cat("Dataset", dataset,  "has already been processed, skipping to next\n")
    next
  }


  if (grepl("TARGET-", dir)){
    print("TARGET detected - need some other considerations for data processing.
          Setting flags to ensure proper functioning.")
    grep = TRUE
    matching_string <- "Primary Blood Derived Cancer"
  }else{
    print("Normal TCGA dataset detected, setting normal flags.")
    grep = FALSE
    matching_string <- "Primary Tumor"
  }

  # load in the data from the associated directories
  # we immediately subset to the first entry of the first row
  # as that contains the actual metadata
  RNA_info <- readRDS(file.path(dir,
                                "transcriptome_query.Rds"))[[1]][[1]]
  CNV_info <- readRDS(file.path(dir,
                                "CNV_query.Rds"))[[1]][[1]]
  clinical_info <- readRDS(file.path(dir,
                                     "clinical.Rds"))

  # adds patient ID for easier control
  RNA_info$patient <- sapply(RNA_info$cases,
                             FUN = get_GDC_idents)
  CNV_info$patient <- sapply(CNV_info$cases,
                             FUN = get_GDC_idents)
  clinical_info$patient <- sapply(clinical_info$submitter_id,
                                  FUN = get_GDC_idents)
  # collate all unique IDs
  all_IDs <- append(append(
    clinical_info$patient,
    CNV_info$patient),
    RNA_info$patient)
  unique_patients <- unique(all_IDs)

  # checks the presence and location for each data modality
  transcriptome_presence <-
    check_data_presence(dataframe = RNA_info,
                        matching_string = matching_string,
                        base_path = file.path(dir,
                                              "Transcriptome_Profiling",
                                              "Gene_Expression_Quantification"),
                        col_name = "RNA_path",
                        grep = grep)
  primary_CNV_presence <-
    check_data_presence(dataframe = CNV_info,
                        matching_string = matching_string,
                        base_path = file.path(dir,
                                              "Copy_Number_Variation",
                                              "Copy_Number_Segment"),
                        col_name = "primary_CNV_path",
                        grep = grep)

  # not doing anything with the healthy CNV atm, gonna leave it out
  # healthy_CNV_presence <-
  #   check_data_presence(dataframe = CNV_info,
  #                       matching_string = "Blood Derived Normal",
  #                       base_path = file.path(dir,
  #                                             "Copy_Number_Variation",
  #                                             "Copy_Number_Segment"),
  #                       col_name = "healthy_CNV_path")

  # merges the 3 separate frames into 1, only if ID is present in all frames
  data_presence <- Reduce(function(x, y) merge(x, y, by = "ID", all = FALSE),
                       list(transcriptome_presence,
                         primary_CNV_presence)
                         # healthy_CNV_presence)
                       )

  # calculates the aneuploidy for each sample
  aneuploidy_scores <- purrr::map_dbl(data_presence$primary_CNV_path, function(path) {
    cnv_data <- load_CNV(path)
    calc_aneuploidy(seg_means = cnv_data$Segment_Mean,
                    lengths = cnv_data$End - cnv_data$Start,
                    base_ploidy = 2,
                    rounded = FALSE)
  })
  data_presence$aneuploidy <- aneuploidy_scores

  # save the combined object for later use in plotting
  saveRDS(data_presence,
          file = file.path(dir,
                           "filtered_data.Rds"))
  sprintf("Finished preprocessing of %s", basename(dir))
}
