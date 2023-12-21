#' @title Extract the first field of the column labels
#'
#' @description Pull out the first field which corresponds to timepoint of non-stationary labeling experiment
#'
#' @param x string (1) Input string for first field extraction
#'
#' @returns Returns the first field from a string separated by "\_"
#'
#' @export

data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)
