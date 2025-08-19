rm(list = ls())
here::i_am("scripts/EDA_TARGET_dataset.R")

# gets some utils
source(file.path(here::here(),
                 "scripts/util_funcs.R"))

# regex to warm up my coding
dir <- grep(pattern = "TARGET-",
                  x = list.dirs(file.path(here::here(),
                                          "data"),
                                recursive = FALSE),
                  value = TRUE)

# if multiple matches grab the first - this can be part of the grep most likely
if (length(dir) > 1){
  dir <- dir[1]
}

# actual exploratory data analysis (EDA)
RNA_info <- readRDS(file.path(dir,
                              "transcriptome_query.Rds"))[[1]][[1]]
CNV_info <- readRDS(file.path(dir,
                              "CNV_query.Rds"))[[1]][[1]]
clinical_info <- readRDS(file.path(dir,
                                   "clinical.Rds"))

RNA_sample_types <- unique(RNA_info$sample_type)
CNV_sample_types <- unique(CNV_info$sample_type)

sample_matching_string <- "Primary Blood Derived Cancer"

RNA_info$patient <- sapply(RNA_info$cases,
                           FUN = get_GDC_idents)
CNV_info$patient <- sapply(CNV_info$cases,
                           FUN = get_GDC_idents)
clinical_info$patient <- sapply(clinical_info$submitter_id,
                                FUN = get_GDC_idents)

primary_CNV_paths <- check_data_presence(dataframe = CNV_info,
                                 matching_string = sample_matching_string,
                                 base_path = file.path(dir,
                                                       "Copy_Number_Variation",
                                                       "Copy_Number_Segment"),
                                 col_name = "path",
                                 grep = TRUE)
primary_trascriptome_presence <-

# calculates the aneuploidy for each sample in the CNV dataset
aneuploidy_scores <- purrr::map_dbl(CNV_paths$path, function(path) {
  cnv_data <- load_CNV(path)
  calc_aneuploidy(seg_means = cnv_data$Segment_Mean,
                  lengths = cnv_data$End - cnv_data$Start,
                  base_ploidy = 2,
                  rounded = FALSE)
})

# ran the dataset through the analysis pipeline, inspecting the resulting data here
rm(list = ls())
filtered_data <- readRDS(file.path(here::here(),
                  "data",
                  "TARGET-AML",
                  "filtered_data.Rds"))
