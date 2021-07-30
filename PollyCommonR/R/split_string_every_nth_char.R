#' split_string_every_nth_char
#'
#' This function splits a string after every nth position
#'
#' @param string A character string or vector of strings
#' @param n The number to split string after every nth position
#' @return A vector of splitted string
#' @examples 
#' split_string_every_nth_char(string)
#' @import stringr
#' @export
split_string_every_nth_char <- function(string = NULL, n = 20){
  message("Split String Every Nth Char Started...")
  require(stringr)
  
  if (identical(string, NULL)) {
    warning("The string is null")
    return(NULL) 
  }
  
  if (!identical(class(string), "character")) {
    warning("The string is not a character")
    return(NULL) 
  }
  
  split_str <- function(string = NULL, n = 20){
    string_splits_count <- round(nchar(string)/n) + 1
    start_vec <- 0:string_splits_count * n + 1
    stop_vec <- 1:(string_splits_count + 1) * n
    string_split <- stringr::str_sub(string = string, start = start_vec, end = stop_vec)
    string_split <- string_split[!is.na(string_split) & string_split != ""]
    return (string_split)
  }
  splitted_string <- sapply(string, split_str, n = n)
  
  message("Split String Every Nth Char Completed...")
  
  return (splitted_string)
}