## use this dataset
# [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50499]
mov10 <- readr::read_tsv("data-raw/GSE50499_raw_counts_GRCh38.p13_NCBI.tsv")

#
cpm_filter <- rowMeans(mov10[, 2:9] / (colSums(mov10) / 1e6)) > 2
mov10 <- mov10[cpm_filter, ]

colnames(mov10) <- c(
  "geneid",
  "low_1",
  "low_2",
  "high_1",
  "high_2",
  "high_3",
  "control_1",
  "control_2",
  "control_3"
)

usethis::use_data(mov10, overwrite = TRUE)
