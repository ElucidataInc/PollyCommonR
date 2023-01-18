#' create_supertest_sets
#'
#' It provides the sets of common ids shared between different sets in upset plot.
#'
#' @param sets_list A list of sets to be compared.
#' @examples 
#' create_supertest_sets(sets_list)
#' @import SuperExactTest stringr
#' @export
create_supertest_sets <- function(sets_list){
  message("Create Supertest Sets Started...")
  
  overall_supertest_summary <- data.frame()
  sets_list <- sets_list[!is.na(names(sets_list))]
  
  names(sets_list) <- stringr::str_trim(names(sets_list))
  if (any(as.logical(lapply(sets_list, function(x) identical(any(stringr::str_detect(x, ",")), TRUE))))){
    message('Some of the elements contain ",", so replacing "," by "." special character.')
    sets_list <- lapply(sets_list, function(x) stringr::str_replace_all(x, pattern = ",", replacement = "."))
  }
  
  super_test = SuperExactTest::supertest(sets_list)
  supertest_summary <- summary(super_test)$Table
  rownames(supertest_summary) <- 1:nrow(supertest_summary)
  
  get_setdiff <- function(degree_index){
    intersect_name <- supertest_summary[degree_index,]$Intersections
    intersect_names_type <<- unique(c(intersect_name, stringr::str_trim(stringr::str_split(intersect_name, pattern = "&")[[1]])))
    others_to_intersect_names <<- supertest_summary[-degree_index,]$Intersections
    others_to_intersect_index <<- as.numeric(rownames(supertest_summary[-degree_index,]))
    others_to_intersect_names_u <- lapply(stringr::str_split(others_to_intersect_names, pattern = "&"), function(x) stringr::str_trim(x))
    other_to_unique_names_index <- others_to_intersect_index[!(as.logical(lapply(others_to_intersect_names_u, function(x) all(x %in% intersect_names_type))))]                                      
    seleceted_to_name <- unique(unlist(lapply(stringr::str_split(supertest_summary[degree_index,]$Elements, pattern = ", "), function(x) stringr::str_trim(x))))
    others_to_name <- unique(unlist(lapply(stringr::str_split(supertest_summary[other_to_unique_names_index,]$Elements, pattern = ", "), function(x) stringr::str_trim(x))))
    
    return (setdiff(seleceted_to_name, others_to_name))
    
  }
  
  single_degree_index <- as.numeric(rownames(supertest_summary))
  unique_to_name_elements <- lapply(single_degree_index, function(x) get_setdiff(x))
  unique_to_name_names <- paste("unique", supertest_summary[single_degree_index,]$Intersections, sep = " ")
  degree_ids <- supertest_summary[single_degree_index,]$Degree
  observed_overlap <- as.numeric(lapply(unique_to_name_elements, function(x) length(x)))
  unique_to_name_elements_chr <- as.character(lapply(unique_to_name_elements, function(x) paste(x, collapse = ", ")))
  unique_to_name_summary <- data.frame(Intersections = unique_to_name_names, Degree = degree_ids, Observed.Overlap = observed_overlap, Elements = unique_to_name_elements_chr, stringsAsFactors = FALSE)
  
  overall_supertest_summary <- rbind(supertest_summary, unique_to_name_summary)
  
  message("Create Supertest Sets Completed...")
  
  return (overall_supertest_summary)
}