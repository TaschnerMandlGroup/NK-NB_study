#!/usr/bin/env Rscript 

# ___R packages___
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library("SCpubr"))
suppressPackageStartupMessages(library("scCustomize"))
suppressPackageStartupMessages(library("ggVennDiagram"))
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("miloR"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("SingleR"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stats"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("DElegate")) #
suppressPackageStartupMessages(library("clusterProfiler")) #
suppressPackageStartupMessages(library("org.Hs.eg.db")) # 
suppressPackageStartupMessages(library("enrichplot")) #
suppressPackageStartupMessages(library("AnnotationDbi")) #
suppressPackageStartupMessages(library("pathview")) #
suppressPackageStartupMessages(library("knitr")) #
suppressPackageStartupMessages(library("pathview")) #
suppressPackageStartupMessages(library("aggregateBioVar")) #
suppressPackageStartupMessages(library("SummarizedExperiment")) 
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("magrittr")) #
suppressPackageStartupMessages(library("pheatmap")) # 
suppressPackageStartupMessages(library("ggtext")) #
suppressPackageStartupMessages(library("ggbeeswarm")) #
suppressPackageStartupMessages(library("ggthemes")) #
suppressPackageStartupMessages(library("slingshot")) #
suppressPackageStartupMessages(library("reshape2")) #
suppressPackageStartupMessages(library("EnhancedVolcano")) 
suppressPackageStartupMessages(library("msigdbr")) #

# ___NK sub-set markers___

## ___NK functional markers from: Yang et al., 2019___

yang2019_clust <- 
  read_excel("/home/rstudio/mnt_resources/NK_cells_gene_sets/Yang2019_supp_Table_1.xlsx")

yang2019_clust_mark <- yang2019_clust$gene
length(yang2019_clust_mark)


## ____NK functional markers from: Crinier et al., 2021____

adaptive_NK_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/adaptive_NK_marker_genes.txt", header = FALSE)
adaptive_NK_genes <- adaptive_NK_genes$V1
conventional_NK_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/conventional_NK_marker_genes.txt", header = FALSE)
conventional_NK_genes <- conventional_NK_genes$V1
mature_NK_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/mature_NK_marker_genes.txt", header = FALSE)
mature_NK_genes <- mature_NK_genes$V1
NKP_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/NKP_marker_genes.txt", header = FALSE)
NKP_genes <- NKP_genes$V1

crinier_et_al_genes <- list(adaptive_NK = adaptive_NK_genes,
                            conventional_NK = conventional_NK_genes,
                            mature_NK = mature_NK_genes,
                            progenitor_NK = NKP_genes)


## ____NK functional markers from: Hannah et al., 2004____

library(readr)
CD56bright_CD16neg_NK_vs_CD56dim_CD16pos_NK_genes <- 
  read_delim("~/mnt_resources/NK_cells_gene_sets/CD56bright_CD16neg_NK_vs_CD56dim_CD16pos_NK_genes.csv",
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE, skip = 2)


activated_CD16pos_NK_vs_CD56dim_CD16pos_NK_genes <- 
  read_delim("~/mnt_resources/NK_cells_gene_sets/activated_CD16pos_NK_vs_CD56dim_CD16pos_NK_genes.csv",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE, skip = 2)

activated_CD16pos_NK_vs_CD56bright_CD16neg_NK_genes <- 
  read_delim("~/mnt_resources/NK_cells_gene_sets/activated_CD16pos_NK_vs_CD56bright_CD16neg_NK_genes.csv",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE, skip = 2)


# !he following features were not be found: 
# STK12, KIAA0101, E2-EPF, CDW52, TSC, CHC1L, CDC2, HSPC150, DAF, FLJ10517, FLJ10461, HN1, SC1, MGC14833, MGC861, TUBB2, KIAA0042
# these are the correct ones:
# AURKB, PCLAF, UBE2S, CD52, TSC1,TSC2, RCBTB2, CDK1, UBE2T, CD55, ASPM, ECT2, JPT1, SPARCL1, UQCC2, CENPM, TUBB2B, KIF14

activated_CD16pos_NK_vs_CD56dim_CD16pos_NK_genes <- rbind(activated_CD16pos_NK_vs_CD56dim_CD16pos_NK_genes, 
                                                          "AURKB", "PCLAF", "UBE2S", "CD52", "TSC1","TSC2", 
                                                          "RCBTB2", "CDK1", "UBE2T", "CD55", "ASPM", "ECT2", 
                                                          "JPT1", "SPARCL1", "UQCC2", "CENPM", "TUBB2B", "KIF14")

## ____NK functional markers from: Blanquart et al., 2024___

inflammed_NK_genes <- read_delim("~/mnt_resources/NK_cells_gene_sets/inflammed_NK_genes.csv",
                                 delim = "\t", escape_double = FALSE,
                                 col_names = FALSE, trim_ws = TRUE, skip = 2)

# The following features were not be found:
# LINC-PINT, LINC01578
# these are the correct ones:
# PINT, CHASERR <- still not found!

inflammed_NK_genes <- rbind(inflammed_NK_genes, "PINT", "CHASERR")


## ____NK other markers____

# direct_R = c("TNFRSF10B","TNFSF6")
NK_functional_markers <- list(canonical_markers = c("EOMES","ID2","TBX21"),
                              NK_subtype_markers = c("FCGR3A","NCAM1","CD160",
                                                     "CD52","IL2RB", "CD27", 
                                                     "CD69", "IL7R"),
                              indirect = c("PRF1","GZMB"),   # Indirect killing of tumor cells
                              direct = c("TNFSF10","FASLG"), # Direct killing of tumor cells
                              recruitment = c("IFNG","TNF")) # Recruitment of other immune cell types such as DCs, T-cells


NK_activating_cytokines <- c("IL2", "IL15", "IL12A", 
                            "IL12B","IL18")

NK_maturation_genes <- c("IFI44L", "IFI6", "IFIT3", "IFI44", "HLA-DPB1", 
                         "HLA-DPA1", "HLA-DRB5", "HLA-DRB1", "ZEB2", "KLF2")

cytokines <- c("IL7","KITLG","FLT3","LTA","LTB","IL15","IL17A","TGFB1")

TFs <- c("EOMES", "ID2", "TBX21", "NFIL3", "HNF1A", "ETS1", "STAT5A", "STAT5B", 
         "TBX21", "TOX", "TOX2","PRDM1","ZEB2","GATA3","SMAD4","FOXO1")

# ___NK subset genes from___
# Activated CD16+ NK vs. CD56dimCD16+ NK DEGs (top 50)
activated_NK_cells <- c("LTB","GZMK","CTSW","AURKB","ANXA2","UBE2C","TK1","CKS2",
                        "PCLAF","STMN1","FEN1","KPNA2","NCF4","UBE2S","CD52",
                        "BIRC5","MAD2L1","SLC1A5","DHCR24","TESC","MAL","RCBTB2",
                        "HLA-DRA","CDK1","DPP4","UBE2T","TYMS","LGALS3","LGALS1",
                        "CD55" ,"ZWINT","SHMT2","PRDX1","PCNA","MMP25","HLA-DMA",
                        "ASPM" ,"ECT2" ,"CRIP1","NUDT1","JPT1","TCF19","UQCC2",
                        "CENPM","RANBP1","TUBB4B","SDF2L1","PSMD1","KIF14", 
                        "HTATSF1")


# Mei et al., 2024
# CD56dim-2 represented a terminally differentiated cell state 
# CD56dim-2 also show high expression of KLRK1, which is one of the main NK cell 
# activating receptor involved in anti-tumor activity 

NK_Mei_et_al_markers <- list(CD56bright = c("SELL", "TCF7", "SPTSSB", "CCR7","NCAM1"),
                             CD56dim_1 = c("GZMB","FCGR3A", "SPON2","PRF1"),
                             CD56dim_2 = c("PPM1L","AOAH","SYNE2", "ZEB2","PDE4D", "KLRK1","HAVCR2"),
                             Immature_NK = c("CCL3","XCL1","CXCR6","IL2RB","LY9","CD7","KLRB1"),
                             NKT_cells = c("GZMH", "CD3D", "CD8A", "CD8B", "NKG7", "NCAM1", "TBX21"),
                             Active_NK_cells = c("CD7", "CXCR4", "IL2RB", "CD69", "KLRB1"))
                             
NK_Mei_et_al_markers_short <- list(CD56bright = c("SELL", "TCF7", "CCR7"),
                             CD56dim_1 = c("FCGR3A", "SPON2"),
                             CD56dim_2 = c("ZEB2", "KLRK1","HAVCR2"),
                             NKT_cells = c("GZMH", "CD8A", "CD8B", "NKG7", "NCAM1", "TBX21"),
                             Active_NK_cells = c("CD7", "CXCR4", "IL2RB", "CD69", "KLRB1"))

# NK and NB cell secretion genes

secreted_factors_low <- c("C5", "CD40LG", "CSF3","CSF2", "CXCL1", "CCL1",
                          "IL1A", "IL1B", "IL1RN", "IL2","IL4","IL5","IL6",
                           "IL10","IL12B", "IL13", "IL17A","IL25",
                          "IL18","IL21","IL27", "CXCL10","CXCL11","CCL2",
                          "CXCL12", "SERPINE1", "TREM1")


secreted_factors_high <- c("ICAM1","IFNG","IL16","IL32",
                           "MIF", "CCL3", "CCL4", "CCL5","CXCL8", "TNF")


secreted_factors <- c("C5", "CD40LG", "CSF3","CSF2", "CXCL1", "CCL1",
                          "IL1A", "IL1B", "IL1RN", "IL2","IL4","IL5","IL6",
                          "IL10","IL12B", "IL13", "IL17A","IL25",
                          "IL18","IL21","IL27", "CXCL10","CXCL11","CCL2",
                          "CXCL12", "SERPINE1", "TREM1",
                      "ICAM1","IFNG","IL16","IL32",
                      "MIF", "CCL3", "CCL4", "CCL5","CXCL8", "TNF")

secreted_factors_MinusMIF <- c("C5", "CD40LG", "CSF3","CSF2", "CXCL1", "CCL1",
                      "IL1A", "IL1B", "IL1RN", "IL2","IL4","IL5","IL6",
                      "IL10","IL12B", "IL13", "IL17A","IL25",
                      "IL18","IL21","IL27", "CXCL10","CXCL11","CCL2",
                      "CXCL12", "SERPINE1", "TREM1",
                      "ICAM1","IFNG","IL16","IL32",
                      "CCL3", "CCL4", "CCL5","CXCL8", "TNF")

coculture_secreted_factors <- c("CSF2","CCL1","IFNG", "IL1B", "IL2",  "MIF",
                                "CCL3","CCL5", "CXCL12", "TNF")

## ___Color pallets___
COLOR_CODE_ATAC = c("NB" = "palegreen4",
                    "T_cell" = "#cc0000",   
                    "NKT" = "#0B5394",
                    "Monocytes" = "#39BEB1",
                    "Memory_B_cell" = "#992669",  
                    "B_cell" = "#fc58a9",   
                     "Erythroblasts" = "#113556",   
                     "SC (14)" = "#7f5b22",  
                     "ND" = "darkgray")


COLOR_CODE_FETAHU2023 = c("NB (8)" = "palegreen4",
                          "T (5)" = "#cc0000",   
                          "T (6)" = "darkred", 
                          "T (9)" = "#f44336",
                          "T (18)" = "#ac2020",
                          "NK (4)" = "#0B5394",
                          "M (1)" = "#1f78b4", 
                          "M (2)" = "#66a8de",    
                          "M (10)" = "#39BEB1",
                          "M (15)"= "#124e77",
                          "B (19)" = "#992669",  
                          "B (3)" = "#fc58a9",   
                          "B (7)" = "#9e3b75",   
                          "B (11)" = "#7f0d71",   
                          "B (16)" = "#8a1bca",   
                          "pDC (12)" = "#15acbc",   
                          "E (13)" = "#113556",   
                          "SC (14)" = "#7f5b22",  
                          "SC (17)" = "#997439",
                          "SC (20)" = "#573807",  
                          "other (21)" = "#000000")

COLOR_CODE_FETAHU2023_META = c("NB (8)" = "darkgreen",
                          "T-cells" = "darkblue",   
                          "NK-cells" = "#b4d1ee",
                          "Myeloid cells" = "#FF0080",
                          "B-cells" = "purple",   
                          "Dendritic cells" = "#FAA21B",   
                          "Erythroid cells" = "#CDA646",   
                          "Progenitor cells" = "brown",  
                          "Unknwon" = "#000000")

COLOR_CODE_FETAHU2023_META_v2 = c("T-cells" = "darkblue",   
                               "NK-cells" = "#b4d1ee",
                               "Myeloid cells" = "#FF0080",
                               "B-cells" = "purple",   
                               "Dendritic cells" = "#FAA21B",   
                               "Erythroid cells" = "#CDA646",   
                               "Progenitor cells" = "brown",  
                               "Unknwon" = "#000000")

COLOR_CODE_FETAHU2023_v2  = c("NB (8)" = "#006400",                  # darkgreen
                                               "T (5)" = "#00008B",                   # darkblue variants
                                               "T (6)" = "#000099",
                                               "T (9)" = "#0000B8",
                                               "T (18)" = "#0011AA",
                                               "NK (4)" = "#b4d1ee",
                                               "M (1)" = "#FF0080",                   # hot pink variants
                                               "M (2)" = "#E60073",
                                               "M (10)" = "#FF3399",
                                               "M (15)" = "#CC0066",
                                               "B (19)" = "#800080",                  # purple variants
                                               "B (3)" = "#9932CC",
                                               "B (7)" = "#A020F0",
                                               "B (11)" = "#BA55D3",
                                               "B (16)" = "#9400D3",
                                               "pDC (12)" = "#FAA21B",
                                               "E (13)" = "#CDA646",
                                               "SC (14)" = "#A52A2A",                 # brown variants
                                               "SC (17)" = "#8B4513",
                                               "SC (20)" = "#5C4033",
                                               "other (21)" = "#000000")              # unknown


COLOR_CODE_FETAHU2023_META_v2 = c("NB (8)" = "palegreen4",
                               "T-cells" = "darkred",   
                               "NK-cells" = "#0B5394",
                               "Myeloid cells" = "#39BEB1",
                               "B-cells" = "#fc58a9",   
                               "Dendritic cells" = "#15acbc",   
                               "Erythroid cells" = "#113556",   
                               "Progenitor cells" = "#573807",  
                               "Unknwon" = "#000000")


COLOR_CODE <- c("I"   = "#7EA3AC",
                "III" = "#4f4e92",
                "II"  = "#A7563C",
                "IV"  = "#FAA21B")

COLOR_CODE_v2 <- c("hNK_Bm1"  = "#6376AE",
                   "hNK_Bm2"  = "#679F8F",
                   "hNK_Bm3"  = "#BF4E4E",
                   "hNK_Bm4"  = "#CDA646")

COLOR_CODE_MiloR <- c("hNK_Bm1"  = "#6376AE",
                   "hNK_Bm2"  = "#679F8F",
                   "hNK_Bm3"  = "#BF4E4E",
                   "hNK_Bm4"  = "#CDA646",
                   "Mixed" = "#520c61")

COLOR_CODE_CCC <- c("hNK_Bm1"  = "#6376AE",
                   "hNK_Bm2"  = "#679F8F",
                   "hNK_Bm3"  = "#BF4E4E",
                   "hNK_Bm4"  = "#CDA646",
                   "NB (8)"   = "#520c61")

COLOR_CODE_CCC_2 <- c("hNK_Bm1"  = "#6376AE",
                      "hNK_Bm2"  = "#679F8F",
                      "hNK_Bm3"  = "#BF4E4E",
                      "hNK_Bm4"  = "#CDA646",
                      "NB"   = "#520c61")

COLOR_CODE_v3 <- c("I"   = "#7EA3AC",
                   "III"  = "#4f4e92",
                   "II" = "#A7563C",
                   "IV"  = "#FAA21B",
                    "V"  = "#575C55")


COLOR_CODE_v4 <- c("hNK_Bm1"  = "#6376AE",
                   "hNK_Bm2"  = "#679F8F",
                   "hNK_Bm3"  = "#BF4E4E",
                   "hNK_Bm4"  = "#CDA646",
                   "M (1)"    = "#520c61",
                   "M (2)"    = "#871F78",
                   "M (10)"   = "#ac87c1",
                   "M (15)"   = "#cd8ec0",
                   "pDC (12)" = "#efa26e")

COLOR_CODE_GROUP <- c("control"   = "#7EA3AC",
                      "ATRXdel"   = "#4f4e92",
                      "MYCNamp"       = "#A7563C",
                      "ATRXwtMYCNwt"  = "#FAA21B")

COLOR_CODE_GROUP_v1 <- c("control"   = "#7EA3AC",
                      "ATRXdel"   = "#4f4e92",
                      "MYCNamp"       = "#A7563C",
                      "ATRXwtMYCNwt"  = "#FAA21B")

COLOR_CODE_GROUP_v2 <- c("control"   = "#7EA3AC",
                         "ATRXdel"   = "#4f4e92",
                         "MYCNamp"       = "#A7563C",
                         "ATRXwtMYCNwt"  = "#FAA21B",
                         "adult" = "#575C55")

COLOR_CODE_GROUP_v3 <- c("ATRXdel"   = "#4f4e92",
                         "MYCNamp"       = "#A7563C",
                         "ATRXwtMYCNwt"  = "#FAA21B")

COLOR_CODE_GROUP_v4 <- c("C"   = "#7EA3AC",
                         "A"   = "#4f4e92",
                         "M"       = "#A7563C",
                         "S"  = "#FAA21B")

COLOR_CODE_CLUST_v1 <- c("1"   = "#6376AE",
                         "2"  = "#679F8F",
                         "3" = "#BF4E4E",
                         "4"  = "#CDA646")

COLOR_CODE_CLUST_v2 <- c("0"   = "#6376AE",
                         "1"  = "#679F8F",
                         "2" = "#BF4E4E",
                         "3"  = "#CDA646")

COLOR_CODE_CLUST_v3 <- c("0"   = "#264653",
                         "1"  = "#287271",
                         "2" = "#2a9d8f",
                         "3"  = "#4f4e92",
                         "4" = "#e9c46a",
                         "5" = "#f4a261",
                         "6" = "#e76f51")

COLOR_CODE_DONORS <- c("donor1"   = "#264653",
                       "donor2"  = "#287271",
                       "donor3" = "#2a9d8f",
                       "donor4"  = "#4f4e92",
                       "donor5" = "#e9c46a",
                       "donor6" = "#f4a261",
                       "donor7" = "#e76f51",
                       "donor8" = "#ce9c85")

COLOR_CODE_DOTPLOT <- c("#264653",
                        "#287271",
                        "#2a9d8f",
                        "#4f4e92",
                        "#e9c46a",
                        "#f4a261",
                        "#e76f51")

COLOR_CODE_PATIENTS <-c("2005_1702" = "#b440a3",
                        "2006_2684" = "#ff91ab",
                        "2014_0102" = "#79c220",
                        "2016_1853" = "#f1e899",
                        "2016_2950" = "#2a1a54",
                        "2016_3924" = "#20798b",
                        "2016_4503" = "#6c77c1",
                        "2018_1404" = "#ac87c1",
                        "2018_1625" = "#4a6e6e",
                        "2018_4252" = "#babd2f",
                        "2018_6056" = "#b14d2b",
                        "2019_2495" = "pink",
                        "2019_5022" = "#846b29",
                        "2019_5754" = "#ce9c85",
                        "2020_1288" = "#1c3c41",
                        "2020_1667" = "#6fb1aa")

COLOR_CODE_PATIENTS_v3 <-c("ATRXdel_01" = "#b440a3",
                           "ATRXdel_02" = "#f76cc6",
                           "ATRXdel_03" = "pink",
                           "ctrl_01" = "#4a9e48",
                           "ctrl_02" = "#247a4d",
                           "ctrl_03" = "#30c67c",
                           "ctrl_04" = "#099773",
                           "ctrl_05" = "#2e703b",
                           "MYCNA_01" = "#6c77c1",
                           "MYCNA_02" = "#ac87c1",
                           "MYCNA_03" = "#c471f2",
                           "MYCNA_04" = "#f44369",
                           "sporadic_01" = "#f9a87b",
                           "sporadic_02" = "#f97d5b",
                           "sporadic_03" = "#b14d2b",
                           "sporadic_05" = "#d6ba55")

COLOR_CODE_PATIENTS_v4 <-c("donor1" = "#4a9e48",
                           "donor2" = "#247a4d",
                           "donor3" = "#30c67c",
                           "donor4" = "#099773",
                           "donor5" = "#2e703b")

COLOR_CODE_PATIENTS_v5 <-c("donor1" = "#b440a3",
                           "donor2" = "#f76cc6")

COLOR_CODE_PATIENTS_v6 <-c("ATRXdel_01" = "#b440a3",
                           "ATRXdel_02" = "#f76cc6",
                           "ATRXdel_03" = "pink")

COLOR_CODE_PATIENTS_v7 <-c("2014_0102" = "#4a9e48",
                           "2020_1288" = "#247a4d",
                           "2018_4252" = "#30c67c",
                           "2016_1853" = "#099773",
                           "2016_2950" = "#2e703b")

PATIENT_GROUPS <- c("2014_0102" = "I",
                    "2016_1853" = "I",
                    "2016_2950" = "I",
                    "2018_4252" = "I",
                    "2020_1288" = "I",
                    "2005_1702" = "III",
                    "2016_3924" = "III",
                    "2019_2495" = "III",
                    "2016_4503" = "II",
                    "2018_1404" = "II",
                    "2019_5022" = "II",
                    "2019_5754" = "II",
                    "2006_2684" = "IV",
                    "2020_1667" = "IV",
                    "2018_1625" = "IV",
                    "2018_6056" = "IV")

PATIENT_GROUPS_v2 <- c("2014_0102" = "ctrl_01",
                       "2016_1853" = "ctrl_02",
                       "2016_2950" = "ctrl_03",
                       "2018_4252" = "ctrl_04",
                       "2020_1288" = "ctrl_05",
                       "2005_1702" = "ATRXdel_01", 
                       "2016_3924" = "ATRXdel_02",
                       "2019_2495" = "ATRXdel_03",
                       "2016_4503" = "MYCNA_01",
                    "2018_1404" = "MYCNA_02",
                    "2019_5022" = "MYCNA_03",
                    "2019_5754" = "MYCNA_04",
                    "2006_2684" = "sporadic_01",
                    "2018_1625" = "sporadic_02",
                    "2018_6056" = "sporadic_03",
                    "2020_1667" = "sporadic_05")

COLOR_CODE_PATIENTS_GROUP <-c("C1" = "#7EA3AC",
                              "C2" = "#7EA3AC",
                              "C3" = "#7EA3AC",
                              "C4" = "#7EA3AC",
                              "C5" = "#7EA3AC",
                              "A1" = "#4f4e92",
                              "A2" = "#4f4e92",
                              "A3" = "#4f4e92",
                              "M1" = "#A7563C",
                              "M2" = "#A7563C",
                              "M3" = "#A7563C",
                              "M4" = "#A7563C",
                              "S1" = "#FAA21B",
                              "S2" = "#FAA21B",
                              "S3" = "#FAA21B",
                              "S5" = "#FAA21B")

COLOR_CODE_PATIENTS_GROUP_v2 <- c("2014_0102" = "#7EA3AC",
                       "2016_1853" = "#7EA3AC",
                       "2016_2950" = "#7EA3AC",
                       "2018_4252" = "#7EA3AC",
                       "2020_1288" = "#7EA3AC",
                       "2005_1702" = "#4f4e92", 
                       "2016_3924" = "#4f4e92",
                       "2019_2495" = "#4f4e92",
                       "2016_4503" = "#A7563C",
                       "2018_1404" = "#A7563C",
                       "2019_5022" = "#A7563C",
                       "2019_5754" = "#A7563C",
                       "2006_2684" = "#FAA21B",
                       "2018_1625" = "#FAA21B",
                       "2018_6056" = "#FAA21B",
                       "2020_1667" = "#FAA21B")


PATIENT_COLORS <- c("ctrl_01" = "#7EA3AC",
                    "ctrl_02" = "#7EA3AC",
                    "ctrl_03" = "#7EA3AC",
                    "ctrl_04" = "#7EA3AC",
                    "ctrl_05" = "#7EA3AC",
                    "ATRXdel_01" = "#4f4e92",
                   "ATRXdel_02" = "#4f4e92",
                   "ATRXdel_03" = "#4f4e92",
                   "MYCNA_01" = "#A7563C",
                   "MYCNA_02" = "#A7563C",
                   "MYCNA_03" = "#A7563C",
                   "MYCNA_04" = "#A7563C",
                   "sporadic_01" = "#FAA21B",
                   "sporadic_02" = "#FAA21B",
                   "sporadic_03" = "#FAA21B",
                   "sporadic_05" = "#FAA21B")

PATIENT_COLORS_v2 <- c( "C1" = "#4081a5",
                    "C2" = "#64a0c7",
                    "C3" = "#79a4c6",
                    "C4" = "#7ea3ac",
                    "C5" = "#5d8197",
                    "A1" = "#9481b2",
                    "A2" = "#737fb4",
                    "A3" = "#6869a7",
                    "M1" = "#c06448",
                    "M2" = "#ad5a40",
                    "M3" = "#965553",
                    "M4" = "#8f3c3c",
                    "S1" = "#fcb11a",
                    "S2" = "#fbaa1f",
                    "S3" = "#f8971f",
                    "S5" = "#f78d21")

COLOR_CODE_UMAP <- c("NB (8)" = "#D3D3D3",
                          "T (5)" = "#D3D3D3",   
                          "T (6)" = "#D3D3D3", 
                          "T (9)" = "#D3D3D3",
                          "T (18)" = "#D3D3D3",
                          "NK (4)" = "#b4d1ee",
                          "M (1)" = "#D3D3D3", 
                          "M (2)" = "#D3D3D3",    
                          "M (10)" = "#D3D3D3",
                          "M (15)"= "#D3D3D3",
                          "B (19)" = "#D3D3D3",  
                          "B (3)" = "#D3D3D3",   
                          "B (7)" = "#D3D3D3",   
                          "B (11)" = "#D3D3D3",   
                          "B (16)" = "#D3D3D3",   
                          "pDC (12)" = "#D3D3D3",   
                          "E (13)" = "#D3D3D3",   
                          "SC (14)" = "#D3D3D3",  
                          "SC (17)" = "#D3D3D3",
                          "SC (20)" = "#D3D3D3",  
                          "other (21)" = "#D3D3D3")

COLOR_CODE_UMAP_v2 <- c("NB (8)" = "darkgreen",
                     "T (5)" = "red",   
                     "T (6)" = "darkred", 
                     "T (9)" = "red3",
                     "T (18)" = "red4",
                     "NK (4)" = "#b4d1ee",
                     "M (1)" = "blue", 
                     "M (2)" = "blue3",    
                     "M (10)" = "darkblue",
                     "M (15)"= "blue4",
                     "B (19)" = "purple",  
                     "B (3)" = "purple2",   
                     "B (7)" = "purple3",   
                     "B (11)" = "purple4",   
                     "B (16)" = "purple",   
                     "pDC (12)" = "orange",   
                     "E (13)" = "darkorange",   
                     "SC (14)" = "gray",  
                     "SC (17)" = "darkgray",
                     "SC (20)" = "lightgray",  
                     "other (21)" = "black")


COLOR_CODE_UMAP_v3 <- c("NB (8)" = "darkgreen",
                        "T (5)" = "darkblue",   
                        "T (6)" = "darkblue", 
                        "T (9)" = "darkblue",
                        "T (18)" = "darkblue",
                        "NK (4)" = "#b4d1ee",
                        "M (1)" = "#FF0080", 
                        "M (2)" = "#FF0080",    
                        "M (10)" = "#FF0080",
                        "M (15)"= "#FF0080",
                        "B (19)" = "purple",  
                        "B (3)" = "purple",   
                        "B (7)" = "purple",   
                        "B (11)" = "purple",   
                        "B (16)" = "purple",   
                        "pDC (12)" = "#FAA21B",   
                        "E (13)" = "#CDA646",   
                        "SC (14)" = "brown",  
                        "SC (17)" = "brown",
                        "SC (20)" = "brown",  
                        "other (21)" = "black")

PATIENT_ORDER <- c("C1", "C2", "C3", "C4", "C5",
                   "A1", "A2", "A3",
                   "M1", "M2", "M3", "M4",
                   "S1", "S2", "S3", "S5")

GROUP_ORDER <- c("control",
                 "ATRXdel",
                 "MYCNamp",
                 "ATRXwtMYCNwt")

CLUSTER_ORDER <- c("hNK_Bm1",
                   "hNK_Bm2",
                   "hNK_Bm3",
                   "hNK_Bm4")

CLUSTER_ORDER_v2 <- c("Mixed",
                      "hNK_Bm4",
                      "hNK_Bm3",
                      "hNK_Bm2",
                      "hNK_Bm1")


GROUP_CLUSTER_ORDER <- c("hNK_Bm1_control", "hNK_Bm1_ATRXdel", "hNK_Bm1_MYCNamp", "hNK_Bm1_ATRXwtMYCNwt",
                         "hNK_Bm2_control", "hNK_Bm2_ATRXdel", "hNK_Bm2_MYCNamp", "hNK_Bm2_ATRXwtMYCNwt", 
                         "hNK_Bm3_control", "hNK_Bm3_ATRXdel", "hNK_Bm3_MYCNamp", "hNK_Bm3_ATRXwtMYCNwt", 
                         "hNK_Bm4_control", "hNK_Bm4_ATRXdel", "hNK_Bm4_MYCNamp", "hNK_Bm4_ATRXwtMYCNwt")

CLUSTER_GROUP_ORDER <- c("control_hNK_Bm1", "control_hNK_Bm2", "control_hNK_Bm3", "control_hNK_Bm4",
                         "ATRXdel_hNK_Bm1", "ATRXdel_hNK_Bm2", "ATRXdel_hNK_Bm3", "ATRXdel_hNK_Bm4",
                         "MYCNamp_hNK_Bm1", "MYCNamp_hNK_Bm2", "MYCNamp_hNK_Bm3", "MYCNamp_hNK_Bm4",
                         "ATRXwtMYCNwt_hNK_Bm1","ATRXwtMYCNwt_hNK_Bm2", "ATRXwtMYCNwt_hNK_Bm3", "ATRXwtMYCNwt_hNK_Bm4")


GROUP_COLORS <- c("C"   = "#7EA3AC",
                  "A"   = "#4f4e92",
                  "M"   = "#A7563C",
                  "S"   = "#FAA21B")

PATIENT_COLORS <-c("C1" = "#7EA3AC",
                   "C2" = "#7EA3AC",
                   "C3" = "#7EA3AC",
                   "C4" = "#7EA3AC",
                   "C5" = "#7EA3AC",
                   "A1" = "#4f4e92",
                   "A2" = "#4f4e92",
                   "A3" = "#4f4e92",
                   "M1" = "#A7563C",
                   "M2" = "#A7563C",
                   "M3" = "#A7563C",
                   "M4" = "#A7563C",
                   "S1" = "#FAA21B",
                   "S2" = "#FAA21B",
                   "S3" = "#FAA21B",
                   "S5" = "#FAA21B")

# Define the patient rename mapping
PATIENTS_RENAME <- c("C1" = "ctrl_01",
                     "C2" = "ctrl_02",
                     "C3" = "ctrl_03",
                     "C4" = "ctrl_04",
                     "C5" = "ctrl_05",
                     "A1" = "ATRXdel_01",
                     "A2" = "ATRXdel_02",
                     "S4" = "ATRXdel_03",
                     "M1" = "MYCNA_01",
                     "M2" = "MYCNA_02",
                     "M3" = "MYCNA_03",
                     "M4" = "MYCNA_04",
                     "S1" = "sporadic_01",
                     "S2" = "sporadic_02",
                     "S3" = "sporadic_03",
                     "S5" = "sporadic_05")
