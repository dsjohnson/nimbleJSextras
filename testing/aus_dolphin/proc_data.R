setwd("~/research/projects/r_packages/nimbleJSextras/testing/aus_dolphin")

library(RMark)

# Read the file
raw_data <- import.chdata("Mark_file_2States.txt", field.names = c("ch", "freq"), header = FALSE)
ch_matrix <- do.call(rbind, strsplit(as.character(raw_data$ch), ""))
ch_matrix <- ifelse(ch_matrix=="A",1,ch_matrix)
ch_matrix <- ifelse(ch_matrix=="B",2,ch_matrix)
ch_matrix <- apply(ch_matrix, 2, as.integer)
saveRDS(ch_matrix, file="aus_dolphin_2state.rds")

raw_data <- import.chdata("Mark_file_3States.txt", field.names = c("ch", "freq"), header = FALSE)
ch_matrix <- do.call(rbind, strsplit(as.character(raw_data$ch), ""))
ch_matrix <- apply(ch_matrix, 2, as.integer)
saveRDS(ch_matrix, file="aus_dolphin_3state.rds")
