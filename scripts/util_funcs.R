library(dplyr)

get_GDC_idents <- function(x){
  ident <- strsplit(x, split = "-")[[1]][3]
  return(ident)
}


create_load_path <- function(row, central_dir){
  load_path <- file.path(central_dir, row$file_id, row$file_name)
  return(load_path)
}


load_transcriptome <- function(filepath){
  RNAseq_data <- read.csv(file = filepath,
                          sep = "\t",
                          header = TRUE,
                          comment.char = "#")
  return(RNAseq_data)
}


load_CNV <- function(filepath){
  CNV_data <- read.csv(file = filepath,
                       sep = "\t")
  return(CNV_data)
}


calc_aneuploidy <- function(seg_means,
                            lengths,
                            rounded = FALSE,
                            base_ploidy = 2) {
  # converts the seg means to copy numbers
  CN <- (2^seg_means) * 2

  # rounds to discrete if flagged - more true to the biology
  if (rounded) {
    CN <- round(CN)
  }

  # calc aneuploidy per bin and length of a bin as total
  bins_aneu <- CN - base_ploidy
  lengths_rescaled <- lengths / sum(lengths)

  # measures the absolute deviation from base ploidy
  aneuploidy <- abs(bins_aneu) * lengths_rescaled
  aneuploidy <- sum(aneuploidy)
  return(aneuploidy)
}


check_data_presence <- function(dataframe,
                                matching_string,
                                base_path,
                                col_name,
                                grep = FALSE){
  # only selects samples matching the given string, possibly with grep
  if (grep){
    selected_samples <- dataframe[grepl(pattern = matching_string,
                                        dataframe$sample_type), ]
  }else{
    selected_samples <- dataframe[dataframe$sample_type == matching_string, ]
  }

  selected_samples <- selected_samples %>%
    dplyr::group_by(patient) %>%
    dplyr::filter(n() == 1) %>%
    dplyr::ungroup()

  data_locs <- data.frame(
    ID = selected_samples$patient,
    load_path = FALSE
  )

  # now need to apply the function so that we get all the load paths for each ID
  data_paths <- apply(selected_samples, 1, function(row){
    file.path(base_path, row["file_id"], row["file_name"])
  })

  # collate together ID and paths
  data_locs$load_path <- data_paths
  colnames(data_locs) <- c("ID", col_name)

  return(data_locs)
}
