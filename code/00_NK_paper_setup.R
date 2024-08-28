# gene lists for NK project

# Yang et al., 2019
library(readxl)
yang2019_clust_mark <- 
  read_excel("/home/rstudio/mnt_resources/NK_cells_gene_sets/Yang2019_supp_Table_1.xlsx")

yang2019_clust_mark <- yang2019_clust_mark$gene
length(yang2019_clust_mark)


## ___NK gene lists___


adaptive_NK_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/adaptive_NK_marker_genes.txt", 
                              header = FALSE)
adaptive_NK_genes <- adaptive_NK_genes$V1
conventional_NK_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/conventional_NK_marker_genes.txt", 
                                  header = FALSE)
conventional_NK_genes <- conventional_NK_genes$V1
mature_NK_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/mature_NK_marker_genes.txt", 
                            header = FALSE)
mature_NK_genes <- mature_NK_genes$V1
NKP_genes <- read.csv("/home/rstudio/mnt_resources/NK_cells_gene_sets/NKP_marker_genes.txt", 
                      header = FALSE)
NKP_genes <- NKP_genes$V1


COLOR_CODE <- c("I"   = "#161B33",
                "II"  = "#419D78",
                "III" = "#AF1B3F",
                "IV"  = "#E88873")

COLOR_CODE_GROUP <- c("control"   = "#161B33",
                "MNA"  = "#419D78",
                "ATRXdel" = "#AF1B3F",
                "sporadic"  = "#E88873")


COLOR_CODE_v2 <- c("hNK_Bm1"   = "#161B33",
                "hNK_Bm2"  = "#419D78",
                "hNK_Bm3" = "#AF1B3F",
                "hNK_Bm4"  = "#E88873")

COLOR_CODE_CLUST_v1 <- c("1"   = "#161B33",
                         "2"  = "#419D78",
                         "3" = "#AF1B3F",
                         "4"  = "#E88873")

COLOR_CODE_CLUST_v2 <- c("0"   = "#161B33",
                         "1"  = "#419D78",
                         "2" = "#AF1B3F",
                         "3"  = "#E88873")


COLOR_CODE_CLUST_v3 <- c("0"   = "#264653",
                         "1"  = "#287271",
                         "2" = "#2a9d8f",
                         "3"  = "#8ab17d",
                         "4" = "#e9c46a",
                         "5" = "#f4a261",
                         "6" = "#e76f51")



COLOR_CODE_DONORS <- c("donor1"   = "#264653",
                       "donor2"  = "#287271",
                       "donor3" = "#2a9d8f",
                       "donor4"  = "#8ab17d",
                       "donor5" = "#e9c46a",
                       "donor6" = "#f4a261",
                       "donor7" = "#e76f51",
                       "donor8" = "#ce9c85")

COLOR_CODE_DOTPLOT <- c("#264653",
                        "#287271",
                        "#2a9d8f",
                        "#8ab17d",
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
                        "2019_2495" = "#c58941",
                        "2019_5022" = "#846b29",
                        "2019_5754" = "#ce9c85",
                        "2020_1288" = "#1c3c41",
                        "2020_1667" = "#6fb1aa")


COLOR_CODE_PATIENTS_v2 <-c("patient 01" = "#b440a3",
                        "patient 02" = "#ff91ab",
                        "patient 03" = "#79c220",
                        "patient 04" = "#f1e899",
                        "patient 05" = "#2a1a54",
                        "patient 06" = "#20798b",
                        "patient 07" = "#6c77c1",
                        "patient 08" = "#ac87c1",
                        "patient 09" = "#4a6e6e",
                        "patient 10" = "#babd2f",
                        "patient 11" = "#b14d2b",
                        "patient 12" = "#c58941",
                        "patient 13" = "#846b29",
                        "patient 14" = "#ce9c85",
                        "patient 15" = "#1c3c41",
                        "patient 16" = "#6fb1aa")

COLOR_CODE_PATIENTS_v3 <-c("ATRXmut_01" = "#b440a3",
                           "ATRXmut_02" = "#f76cc6",
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
                           "sporadic_04" = "#c58941",
                           "sporadic_05" = "#d6ba55")

COLOR_CODE_PATIENTS_v4 <-c("donor1" = "#4a9e48",
                           "donor2" = "#247a4d",
                           "donor3" = "#30c67c",
                           "donor4" = "#099773",
                           "donor5" = "#2e703b")

COLOR_CODE_PATIENTS_v5 <-c("donor1" = "#b440a3",
                           "donor2" = "#f76cc6")

COLOR_CODE_PATIENTS_v6 <-c("ATRXmut_01" = "#b440a3",
                           "ATRXmut_02" = "#f76cc6")


COLOR_CODE_PATIENTS_v7 <-c("2014_0102" = "#4a9e48",
                           "2020_1288" = "#247a4d",
                           "2018_4252" = "#30c67c",
                           "2016_1853" = "#099773",
                           "2016_2950" = "#2e703b")


PATIENT_GROUPS <- c("2005_1702" = "III", 
                    "2006_2684" = "IV",
                    "2014_0102" = "I",
                    "2016_1853" = "I",
                    "2016_2950" = "I",
                    "2016_3924" = "III",
                    "2016_4503" = "II",
                    "2018_1404" = "II",
                    "2018_1625" = "IV",
                    "2018_4252" = "I",
                    "2018_6056" = "IV",
                    "2019_2495" = "IV",
                    "2019_5022" = "II",
                    "2019_5754" = "II",
                    "2020_1288" = "I",
                    "2020_1667" = "IV")


PATIENT_GROUPS_v2 <- c("2005_1702" = "ATRXmut_01", 
                    "2006_2684" = "sporadic_01",
                    "2014_0102" = "ctrl_01",
                    "2016_1853" = "ctrl_02",
                    "2016_2950" = "ctrl_03",
                    "2016_3924" = "ATRXmut_02",
                    "2016_4503" = "MYCNA_01",
                    "2018_1404" = "MYCNA_02",
                    "2018_1625" = "sporadic_02",
                    "2018_4252" = "ctrl_04",
                    "2018_6056" = "sporadic_03",
                    "2019_2495" = "sporadic_04",
                    "2019_5022" = "MYCNA_03",
                    "2019_5754" = "MYCNA_04",
                    "2020_1288" = "ctrl_05",
                    "2020_1667" = "sporadic_05")
