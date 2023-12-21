#' @title Extract the third field of the column labels
#'
#' @description Pull out the third field which corresponds to replicate of sample of non-stationary labeling experiment
#'
#' @param x string (1) Input string for third field extraction
#'
#' @returns Returns the third field from a string separated by "\_"
#'
#' @export

data_cleanIII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 3)
