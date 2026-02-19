
#' @export
collapse_ch <- function(ch_matrix){
  ch_strings <- apply(ch_matrix, 1, paste, collapse = ",")
  counts <- table(ch_strings)
  unique_histories_str <- names(counts)
  unique_matrix <- t(sapply(strsplit(unique_histories_str, ","), as.numeric))
  final_counts <- as.vector(counts)
  return(list(unique_ch=unique_matrix, w=final_counts))
}
