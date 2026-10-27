#' @title Collapse all capture histories to only unique histories
#' @param ch_matrix A matrix of capture histories where each row represents a single individual.
#' @description
#' Collapse all individuals to unique histories and number of times each history occurs
#' @export
collapse_ch <- function(ch_matrix){
  ch_strings <- apply(ch_matrix, 1, paste, collapse = ",")
  counts <- table(ch_strings)
  unique_histories_str <- names(counts)
  unique_matrix <- t(sapply(strsplit(unique_histories_str, ","), as.numeric))
  final_counts <- as.vector(counts)
  return(list(unique_ch=unique_matrix, w=final_counts))
}
