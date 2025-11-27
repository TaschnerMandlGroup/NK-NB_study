# Extract Cell IDs and corresponding cell type annotation

fetahu2023 <- readRDS(file = "/home/rstudio/mnt_out/NK_project/Rds/fetahu2023_withNewMetadata.Rds")


cells <- colnames(fetahu2023)
md <- fetahu2023@meta.data  
head(colnames(md))   

table(md$cell_type_assignment)    
table(Idents(fetahu2023))         

Idents(fetahu2023) <- "cell_type_assignment"  

out <- data.frame(
  submitter_id = md$fetahu_et_al_2023,
  cell_barcode = cells,
  cell_type_assignment = md$orig.ident,
  CL_ontology_id = NA,
  row.names = NULL,
  check.names = FALSE
)

write.table(out,"/home/rstudio/mnt_out/NK_project/tables/fetahu2023_cell_ids_and_celltype.csv", 
          row.names = FALSE,quote = FALSE,sep = "\t")


write.table(out,"/home/rstudio/mnt_out/CancReg/fetahu2023/scPCA_pipeline/config_files/cell-metadata.tsv", 
          row.names = FALSE,quote = FALSE,sep = "\t")
