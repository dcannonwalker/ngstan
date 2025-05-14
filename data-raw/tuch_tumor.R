# tuch data
# clean up the raw data
tuch_tumor <- readxl::read_xls("data-raw/tuch.xls")
names(tuch_tumor)[1:3] <- unlist(tuch_tumor[2, 1:3])
tuch_tumor <- tuch_tumor[3:nrow(tuch_tumor), 1:9]
for (i in seq(3, 9)) {
  tuch_tumor[[i]] <- as.numeric(tuch_tumor[[i]])
}

names(tuch_tumor) <- snakecase::to_snake_case(
  stringr::str_remove(
    names(tuch_tumor), "\\...\\d"
  ),
  abbreviations = paste0(rep(c(8, 33, 51), each = 2), c("N", "T")),
)

usethis::use_data(tuch_tumor, overwrite = TRUE)
