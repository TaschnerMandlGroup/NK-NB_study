DEG_hNK_Bm2_control_vs_ATRXdel <- read.delim("~/mnt_out/NK_project/tables/DEG_hNK-Bm2_control_vs_ATRXdel_2025-04-24.csv")
DEG_hNK_Bm3_control_vs_ATRXdel <- read.delim("~/mnt_out/NK_project/tables/DEG_hNK-Bm3_control_vs_ATRXdel_2025-04-24.csv")

NK2_signif <- DEG_hNK_Bm2_control_vs_ATRXdel[DEG_hNK_Bm2_control_vs_ATRXdel$padj < 0.05, ]
nrow(NK2_signif)
NK2_pos <- NK2_signif[NK2_signif$log_fc > 1, ]
nrow(NK2_pos)
NK2_neg <- NK2_signif[NK2_signif$log_fc < -1, ]
nrow(NK2_neg)


NK3_signif <- DEG_hNK_Bm3_control_vs_ATRXdel[DEG_hNK_Bm3_control_vs_ATRXdel$padj < 0.05, ]
nrow(NK3_signif)
NK3_pos <- NK3_signif[NK3_signif$log_fc > 1, ]
nrow(NK3_pos)
NK3_neg <- NK3_signif[NK3_signif$log_fc < -1, ]
nrow(NK3_neg)
