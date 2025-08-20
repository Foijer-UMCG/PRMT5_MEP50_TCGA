rm(list = ls())
here::i_am("scripts/download_TCGA.R")

# comment out all the projects that are already downloaded,
# otherwise they will be downloaded again!
projects_wanted <- c(
  "TARGET-ALL-P1",
  "TARGET-ALL-P2"
#  "TCGA-PRAD",
#  "TCGA-HNSC",
#  "TCGA-SKCM"
#  "TARGET-AML",
#  "TCGA-LUAD",
#  "TCGA-BRCA",
#  "TCGA-COAD"
)

for (project in projects_wanted){
  project_dir <- file.path(here::here(),
                           "data")

  dir.create(project_dir,
             recursive = TRUE)

  # download the transcriptome information
  query_transcriptome <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification")

  TCGAbiolinks::GDCdownload(
    query = query_transcriptome,
    method = "api",
    directory = project_dir,
    files.per.chunk = 1)

  saveRDS(query_transcriptome,
          file = file.path(project_dir,
                           project,
                           "transcriptome_query.Rds"))

  # clinical information
  clinical_data <- TCGA_BRCA_clinical <- TCGAbiolinks::GDCquery_clinic(
    project = project,
    type = "clinical",
    save.csv = FALSE)

  saveRDS(clinical_data,
          file = file.path(project_dir,
                           project,
                           "clinical.Rds"))

  # and finally the copy number information
  query_CNV <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Copy Number Segment")

  TCGAbiolinks::GDCdownload(
    query = query_CNV,
    method = "api",
    directory = project_dir,
    files.per.chunk = 1)

  saveRDS(query_CNV,
          file = file.path(project_dir,
                           project,
                           "CNV_query.Rds"))
}
