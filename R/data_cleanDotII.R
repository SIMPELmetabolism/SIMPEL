#' @title Extract the second field of labels
#'
#' @description Pull out the second field of character separated dots ('.')
#'
#' @param x string (1) Input string for second field extraction
#'
#' @returns Returns the second field from a string separated by "."
#'
#' @export

data_cleanDotII <- function(x) sapply (strsplit(x , '[.]' ), `[` , 2)
