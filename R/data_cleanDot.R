#' @title Extract the first field of labels
#'
#' @description Pull out the first field of string separated by dots ('.')
#'
#' @param x string (1) Input string for first field extraction
#'
#' @returns Returns the first field from a string separated by "."
#'
#' @export

data_cleanDot <- function(x) sapply (strsplit(x , '[.]' ), `[` , 1)
