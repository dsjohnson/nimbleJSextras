library(RMark)
setwd("~/research/projects/r_packages/nimbleJSextras/testing/scot_dolphin")
# Read the file
raw_data <- import.chdata("juvAd_Robust+design+models.txt", field.names = c("ch", "freq"), header = FALSE)
ch_matrix <- do.call(rbind, strsplit(as.character(raw_data$ch), ""))
ch_matrix <- apply(ch_matrix, 2, as.integer)
intervals <- c(0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0)
p_lengths <- rle(c(intervals, 1))$lengths[rle(c(intervals, 1))$values == 0] + 1
idx <- rep(seq_along(p_lengths), p_lengths)

output <- NULL
for(t in seq_along(p_lengths)){
  tmp <- ch_matrix[,idx==t]
  capt <- rowSums(tmp)
  output <- cbind(output, capt)
}


