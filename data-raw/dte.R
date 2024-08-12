## code to prepare `dte` dataset goes here
ribo <- read.table("data-raw/All.Riboseq.CDS.count.txt", header = T)
rna <- read.table("data-raw/All.RNAseq.CDS.count.txt", header = T)

dte <- merge(ribo, rna)

colnames(dte) <- c(
    "geneid",
    paste0("ribo", rep("_ctrl_dpi5", 4), "_", 1:4),
    paste0("ribo", rep("_trt_dpi5", 5), "_", 1:5),
    paste0("ribo", rep("_ctrl_dpi8", 4), "_", 1:4),
    paste0("ribo", rep("_trt_dpi8", 5), "_", 1:5),
    paste0("rna", rep("_ctrl_dpi5", 4), "_", 1:4),
    paste0("rna", rep("_trt_dpi5", 5), "_", 1:5),
    paste0("rna", rep("_ctrl_dpi8", 4), "_", 1:4),
    paste0("rna", rep("_trt_dpi8", 5), "_", 1:5)
)

usethis::use_data(dte, overwrite = TRUE)

