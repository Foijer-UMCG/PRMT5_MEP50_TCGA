here::i_am("scripts/plot_testing.R")

# loaded in a dataset to try some new ways of presenting plots
plot_genes(filtered_data,
           gene1 = "PRMT5",
           gene2 = "WDR77",
           cancer_type = cancer_type,
           plot_name = file.path(plot_dir,
                                 sprintf("PRMT5_WDR77_corr_%s.pdf", cancer_type)))

p <- ggpubr::ggscatter(filtered_data,
                       x = "PRMT5",
                       y = "WDR77",
                       color = "blue",
                       facet.by = "aneu_factor",
                       size = 0.7,
                       add = "reg.line",
                       conf.int = TRUE,
                       cor.coef = TRUE,
                       cor.coeff.args = list(method = "pearson"),
                       add.params = list(color = "red"),
                       title = sprintf("%s expression against %s expression, %s", "PRMT5", "WDR77", cancer_type),
                       xlab = sprintf("PRTM5 (FPKM)"),
                       ylab = sprintf("WDR77 (FPKM)")) +
  ggpubr::stat_cor(ggplot2::aes(color = aneu_factor)) +
  ggplot2::coord_fixed(ratio = 1)

p

# defines the aneu_division
q <- stats::quantile(filtered_data$aneuploidy,
                     probs = c(0.40, 0.60),
                     na.rm = TRUE)

# and adds as factors in the dataframe
filtered_data$aneu_factor <- as.integer(filtered_data$aneuploidy>= q[2]) +
  as.integer(filtered_data$aneuploidy>= q[1])
filtered_data$aneu_factor <- factor(filtered_data$aneu_factor,
                                    levels = c(0, 2),
                                    labels = c("low", "high"))
filtered_data <- tidyr::drop_na(filtered_data)

# might be better to plot/report statistics this way:
p <- ggpubr::ggscatter(filtered_data,
                       x = "PRMT5",
                       y = "WDR77",
                       color = "aneu_factor",
                       size = 0.7,
                       )

# testing some different plotting looks with already loaded data
cancer_type <- basename(dir)

my_lims <- range(with(filtered_data, c(WDR77, PRMT5)))

p <- ggpubr::ggscatter(filtered_data,
                       x = "PRMT5",
                       y = "WDR77",
                       color = "blue",
                       size = 0.7,
                       add = "reg.line",
                       conf.int = TRUE,
                       cor.coef = TRUE,
                       cor.coeff.args = list(method = "pearson"),
                       add.params = list(color = "red"),
                       title = sprintf("%s expression against %s expression, %s", "PRMT5", "WDR77", cancer_type),
                       xlab = sprintf("PRTM5 (FPKM)"),
                       ylab = sprintf("WDR77 (FPKM)")) +
  ggpubr::stat_cor(ggplot2::aes(color = aneu_factor)) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::coord_cartesian(xlim = my_lims, ylim = my_lims)
p

filtered_data$MEP50 <- filtered_data$WDR77
p <- ggpubr::ggscatter(filtered_data,
                       x = "PRMT5",
                       y = "MEP50",
                       color = "blue",
                       size = 0.7,
                       add = "reg.line",
                       conf.int = TRUE,
                       cor.coef = TRUE,
                       cor.coef.args = list(methods = "pearson"),
                       add.params = list(color = "red")) +
  ggplot2::coord_cartesian(xlim = my_lims, ylim = my_lims)
ggplot2::ggsave(plot = p,
                filename = file.path(here::here(),
                          "test_plot.pdf"),
                unit = "in",
                height = 7,
                width = 7)
p


data_dirs <- list.dirs(file.path(here::here(),
                                 "data"),
                       recursive = FALSE)

# quickly renaming all the WDR77 cols to MEP50
for (dir in data_dirs) {
  filtered_data <- readRDS(file.path(dir,
                    "expression_data.Rds"))
  filtered_data$MEP50 <- filtered_data$WDR77
  filtered_data <- subset(filtered_data, select = -c(WDR77))
  saveRDS(object = filtered_data,
          file = file.path(dir,
                           "expression_data.Rds"))
}

