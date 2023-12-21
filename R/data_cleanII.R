#' @title Extract the second field of the column labels
#'
#' @description Pull out the second field which corresponds to category of sample of non-stationary labeling experiment
#'
#' @param x string (1) Input string for second field extraction
#'
#' @returns Returns the second field from a string separated by "\_"
#'
#' @export

data_cleanII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)
