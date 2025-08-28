rm(list = ls())
here::i_am("scripts/pres_plots.R")

source(file.path(here::here(),
                 "scripts/util_funcs.R"))
source(file.path(here::here(),
                 "scripts/plotting_funcs.R"))

data_dirs <- list.dirs(file.path(here::here(),
                                 "data"),
                       recursive = FALSE)

for (dir in data_dirs){
  file <- list.files(dir,
                     full.names = TRUE,
                     pattern = "^expression_data.Rds$")

  data <- readRDS(file)

  # make all the individual plots required for the cowplot
  p1 <- plot_genes(data,
                   gene1 = "PRMT5",
                   gene2 = "MEP50",
                   plot_name = "",
                   cancer_type = basename(dir),
                   return = TRUE) +
    cowplot::theme_cowplot() +
    ggplot2::ggtitle(NULL)

  p2 <- plot_aneu_gene_scatter(data,
                               gene = "PRMT5",
                               aneuploidy = "aneuploidy",
                               plot_name = "",
                               cancer_type = basename(dir),
                               return = TRUE) +
    cowplot::theme_cowplot() +
    ggplot2::ggtitle(NULL)

  p3 <- plot_aneu_gene_scatter(data,
                               gene = "MEP50",
                               aneuploidy = "aneuploidy",
                               plot_name = "",
                               cancer_type = basename(dir),
                               return = TRUE) +
    cowplot::theme_cowplot() +
    ggplot2::ggtitle(NULL)

  p4 <- cowplot::ggdraw() +
    cowplot::draw_label(basename(dir),
                        fontface = "italic",
                        size = 24)

  compiled_plot <- cowplot::plot_grid(plotlist = list(p1, p2, p3, p4),
                     ncol = 2,
                     nrow = 2,
                     labels = c("A", "B", "C", ""))

  ggplot2::ggsave(filename = file.path(here::here(),
                                       "plots",
                                       basename(dir),
                                       "combined_plot.pdf"),
                  plot = compiled_plot,
                  width = 14,
                  height = 14,
                  unit = "in")
  ggplot2::ggsave(filename = file.path(here::here(),
                                       "plots",
                                       basename(dir),
                                       "combined_plot.png"),
                  plot = compiled_plot,
                  width = 14,
                  height = 14,
                  unit = "in")
}
