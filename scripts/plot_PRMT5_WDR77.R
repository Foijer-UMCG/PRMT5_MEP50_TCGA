rm(list = ls())
here::i_am("scripts/plot_PRMT5_WDR77.R")

source(file.path(here::here(),
                 "scripts/util_funcs.R"))
source(file.path(here::here(),
                 "scripts/plotting_funcs.R"))


data_dirs <- list.dirs(file.path(here::here(),
                                 "data"),
                       recursive = FALSE)

for (dir in data_dirs) {
  # if (grepl("TARGET-AML", dir)){
  #   print("TARGET-AML detected - need some other considerations for data processing")
  #   next
  # }
  if (!(grepl("TARGET-AML", dir))){
    next
  }

  # reads the pre-filtered dataset
  filtered_data <- readRDS(file.path(dir, "filtered_data.Rds"))

  # extracts both PRMT5 & WDR77 expression
  PRMT5_expression <- purrr::map_dbl(filtered_data$RNA_path, function(path) {
    exp_data <- load_transcriptome(path)
    exp_data[exp_data$gene_name == "PRMT5", ]$fpkm_unstranded
  },
  .progress = "PRMT5_exp_extraction")

  WDR77_expression <- purrr::map_dbl(filtered_data$RNA_path, function(path) {
    exp_data <- load_transcriptome(path)
    exp_data[exp_data$gene_name == "WDR77", ]$fpkm_unstranded
  },
  .progress = "WDR77_exp_extraction")

  filtered_data$PRMT5 <- PRMT5_expression
  filtered_data$WDR77 <- WDR77_expression

  # defines the aneu_division
  q <- stats::quantile(filtered_data$aneuploidy,
                       probs = c(0.30, 0.70),
                       na.rm = TRUE)

  # and adds as factors in the dataframe
  filtered_data$aneu_factor <- as.integer(filtered_data$aneuploidy>= q[2]) +
    as.integer(filtered_data$aneuploidy>= q[1])
  filtered_data$aneu_factor <- factor(filtered_data$aneu_factor,
                                      levels = c(0, 1, 2),
                                      labels = c("low", "middle", "high"))

  cancer_type <- basename(dir)
  plot_dir <- file.path(here::here(),
                        "plots",
                        cancer_type)
  dir.create(plot_dir,
             recursive = TRUE)

  for (gene in c("PRMT5", "WDR77")){
    # makes the density plots and saves them - see plotting_funcs.R for details
    plot_gene_distro(data = filtered_data,
                     gene = gene,
                     cancer_type = cancer_type,
                     plot_name = file.path(plot_dir,
                                           sprintf("%s_distribution_%s.pdf",
                                                   gene,
                                                   cancer_type)))

    # also makes the boxplots
    plot_exp_boxplots(data = filtered_data,
                      gene = gene,
                      cancer_type = cancer_type,
                      plot_name = file.path(plot_dir,
                                            sprintf("%s_boxplot_aneuploidy_split_%s.pdf",
                                                    gene,
                                                    cancer_type)))

    # and the aneuploidy vs. gene plots
    plot_aneu_gene_scatter(data = filtered_data,
                           gene = gene,
                           aneuploidy = "aneuploidy",
                           cancer_type = cancer_type,
                           plot_name = file.path(plot_dir,
                                                 sprintf("%s_against_aneuploidy_%s.pdf",
                                                         gene,
                                                         cancer_type)))
  }

  # finally the PRMT5 vs. WDR77 expression - bit more hardcoded than the rest
  plot_genes(filtered_data,
             gene1 = "PRMT5",
             gene2 = "WDR77",
             cancer_type = cancer_type,
             plot_name = file.path(plot_dir,
                                   sprintf("PRMT5_WDR77_corr_%s.pdf", cancer_type)))
  sprintf("Finished plotting of %s", cancer_type)
}
