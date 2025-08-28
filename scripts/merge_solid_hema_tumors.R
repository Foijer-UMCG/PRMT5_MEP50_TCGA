rm(list = ls())
here::i_am("scripts/merge_solid_hema_tumors.R")

source(file.path(here::here(),
                 "scripts/util_funcs.R"))
source(file.path(here::here(),
                 "scripts/plotting_funcs.R"))

data_dirs <- list.dirs(file.path(here::here(),
                                 "data"),
                       recursive = FALSE)


load_and_tag <- function(dir) {
  # pick the file inside the directory (adjust pattern if needed)
  file <- list.files(dir,
                     full.names = TRUE,
                     pattern = "^expression_data.Rds$")

  df <- readRDS(file)
  # classify type by directory name
  dataset_type <- ifelse(grepl("TARGET|MP2PRT", dir), "hematological", "solid")

  # add metadata columns
  df$dataset_dir  <- basename(dir)
  df$dataset_type <- dataset_type
  df
}

# 2. Load all datasets and merge into one big dataframe
all_data <- lapply(data_dirs, load_and_tag) %>% bind_rows()

for (data_type in c("solid", "hematological", "all")){
  if (data_type ==  "all"){
    data <- all_data
  }else{
    data <- all_data[all_data$dataset_type == data_type, ]
  }
  cancer_type <- data_type
  plot_dir <- file.path(here::here(),
                        "plots",
                        data_type)

  # want to reassign the aneuploidy groups to general
  # defines the aneu_division
  # maybe do this per
  q <- stats::quantile(data$aneuploidy,
                       probs = c(0.30, 0.70),
                       na.rm = TRUE)

  # and adds as factors in the dataframe
  data$aneu_factor <- as.integer(data$aneuploidy>= q[2]) +
    as.integer(data$aneuploidy>= q[1])
  data$aneu_factor <- factor(data$aneu_factor,
                                      levels = c(0, 1, 2),
                                      labels = c("low", "middle", "high"))


  for (gene in c("MEP50", "PRMT5")){
    # also makes the boxplots
    plot_exp_boxplots(data = data,
                      gene = gene,
                      cancer_type = cancer_type,
                      plot_name = file.path(plot_dir,
                                            sprintf("%s_boxplot_aneuploidy_split_%s.pdf",
                                                    gene,
                                                    cancer_type)))

    plot_gene_distro(data = data,
                     gene = gene,
                     cancer_type = cancer_type,
                     plot_name = file.path(plot_dir,
                                           sprintf("%s_distribution_%s.pdf",
                                                   gene,
                                                   cancer_type)))

    plot_aneu_gene_scatter(data = data,
                           gene = gene,
                           aneuploidy = "aneuploidy",
                           cancer_type = cancer_type,
                           plot_name = file.path(plot_dir,
                                                 sprintf("%s_against_aneuploidy_%s.pdf",
                                                         gene,
                                                         cancer_type)))



  }

  plot_genes(data,
             gene1 = "PRMT5",
             gene2 = "MEP50",
             cancer_type = cancer_type,
             plot_name = file.path(plot_dir,
                                   sprintf("PRMT5_MEP50_corr_%s.pdf", cancer_type)))
}

