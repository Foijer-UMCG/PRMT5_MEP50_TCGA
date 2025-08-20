rm(list = ls())
here::i_am("scripts/merge_solid_hema_tumors.R.R")

source(file.path(here::here(),
                 "scripts/util_funcs.R"))
source(file.path(here::here(),
                 "scripts/plotting_funcs.R"))

data_dirs <- list.dirs(file.path(here::here(),
                                 "data"),
                       recursive = FALSE)

