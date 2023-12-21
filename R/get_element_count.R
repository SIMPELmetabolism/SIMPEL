#' @title Calculate number of elements in a given formula
#'
#' @description This function is calculating the number of each of the elements present in the formula
#'
#' @param comp_formula is the formula from the annotation file
#'
#' @return This function returns a list of counts of each element present in the formula
#'
#' @importFrom CHNOSZ makeup
#'
#' @export
#'
#' @examples
#' require(dplyr)
#' get_element_count("C5H5N5")

get_element_count <- function(comp_formula){
  comp_formula %>%
    makeup() %>%
    as.list() %>%
    return()
}
