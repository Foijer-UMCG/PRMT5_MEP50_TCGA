rm(list = ls())
here::i_am("scripts/download_TCGA.R")

# comment out the projects you don't want to download
# there's a check in place to prevent double downloading
projects_wanted <- c(
  "MP2PRT-ALL",
  "TARGET-ALL-P1",
  "TARGET-ALL-P2",
  "TCGA-PRAD",
  "TCGA-HNSC",
  "TCGA-SKCM",
  "TARGET-AML",
  "TCGA-LUAD",
  "TCGA-BRCA",
  "TCGA-COAD"
)

force = FALSE
for (project in projects_wanted){
  project_dir <- file.path(here::here(),
                           "data")

  dir.create(project_dir,
             recursive = TRUE,
             showWarnings = FALSE)

  dest_dir <- file.path(project_dir,
                        project,
                        "Transcriptome_Profiling",
                        "Gene_Expression_Quantification")

  # only downloads if the directory doesn't exist or forced
  if ( !(dir.exists(dest_dir)) | force){
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
  }

  # clinical information, small so no extra checks
  clinical_data <- TCGA_BRCA_clinical <- TCGAbiolinks::GDCquery_clinic(
    project = project,
    type = "clinical",
    save.csv = FALSE)

  saveRDS(clinical_data,
          file = file.path(project_dir,
                           project,
                           "clinical.Rds"))

  dest_dir <- file.path(project_dir,
                        project,
                        "Copy_Number_Variation",
                        "Copy_Number_Segment")

  # only downloads if the directory doesn't exist or forced
  if ( !(dir.exists(dest_dir)) | force){
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

  # and also downloads protein expression data (if available - it's rare)
  dest_dir <- file.path(project_dir,
                        project,
                        "Proteome_Profiling",
                        "Protein_Expression_Quantification")

  # only downloads if the directory doesn't exist or forced
  if (!(dir.exists(dest_dir)) || force) {
    tryCatch({
      query_protein <- TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Proteome Profiling",
        data.type = "Protein Expression Quantification"
      )

      TCGAbiolinks::GDCdownload(
        query = query_protein,
        method = "api",
        directory = project_dir,
        files.per.chunk = 1
      )

      saveRDS(query_protein,
              file = file.path(project_dir,
                               project,
                               "protein_query.Rds"))
    },
    error = function(e) {
      message(sprintf("⚠️ No protein array data available for project %s: %s",
                      project, e$message))
    },
    warning = function(w) {
      message(sprintf("⚠️ Warning while querying protein data for %s: %s",
                      project, w$message))
      invokeRestart("muffleWarning")  # suppresses repeated warnings
    })
  }

}
