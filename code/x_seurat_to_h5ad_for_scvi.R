#!/usr/bin/Rscript

# Install required packages (if needed)
if (!require("formr")) remotes::install_github("rubenarslan/formr")
if (!require("sceasy")) devtools::install_github("cellgeni/sceasy")

# if conda is not installed
#reticulate::install_miniconda()
#reticulate::use_condaenv("r-reticulate")
#py_install("numpy", method = "conda")

# Set the conda environment
#SCVI_PATH <- ""
use_condaenv(SCVI_PATH) 

# Load libraries
library("anndata")
library("formr")
library("Seurat")
library("sceasy")

filepath.in <- "/home/rstudio/mnt_out/Rds/"
filepath.out <- "/home/rstudio/mnt_out/h5ad/"

sample.list <- read.table("/home/rstudio/mnt_out/CancReg/cancer_regression_paper/metadata/internal_cohort_file_names.csv", 
                          header = TRUE)
samples <- sample.list$file_name

# Load Rds file - convert Seurat to annData object - save .h5ad file
for (i in samples) {
  print(i)
  filepath <- file.path(filepath.in,i)
  print(filepath)  
  seurat <- readRDS(filepath)
  adata <- convertFormat(seurat,from="seurat", 
                         to="anndata", 
                         main_layer="counts", 
                         drop_single_values=FALSE)
  write_h5ad(adata, file.path(paste0(filepath.out, i,".h5ad")))
}
