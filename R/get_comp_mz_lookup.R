#' @title Extract information on isotopologues to be searched for
#'
#' @description This function is going to evaluate chemical formulas to determine a set of possible isotopologues
#'
#' @param compound_data is the annotation file
#' @param comp_formula is the formula
#' @param r_time This will be the rt from your annotation file
#' @param ppm is the ppm error to be used to get a range of possible mass
#' @param polarity will be the polarity from both Compound_Data and XCMS_data
#' @param label_scheme defines the labeling scheme of the experiments
#'
#' @return The function returns a list of various information of possible isotopologues including:
#' mass, lower bound and upper bound calculated based on ppm error provided, as well as numbers of
#' isotopes
#'
#' @export

get_comp_mz_lookup <- function(compound_data, comp_formula, r_time, ppm, polarity, label_scheme){

  #### INNITIALIZATION ####
  w_c = 1.003355
  w_n = 0.997035
  w_h = 1.006277

  d = list()
  gln_mass = get_comp_mass(comp_formula, polarity)
  elementInLabel = unlist(strsplit(label_scheme, ""))
  c = unlist(ifelse("C" %in% elementInLabel, list(get_element_count(comp_formula)[['C']]), list(NULL)))
  n = unlist(ifelse("N" %in% elementInLabel, list(get_element_count(comp_formula)[['N']]), list(NULL)))
  h = unlist(ifelse("H" %in% elementInLabel, list(get_element_count(comp_formula)[['H']]), list(NULL)))

  comp_prefix <- compound_data %>%
    dplyr::filter(formula == comp_formula & rt == r_time) %>%
    pull(prefix)

  # create a dataframe containing element info
  elementDf = data.frame(element = c('C','N','H'), w = c(w_c, w_n, w_h),
                         InLabel = as.numeric(c("C","N","H") %in% elementInLabel),#lengths(list(c, n, h))
                         number = unlist(replace(list(c,n,h), unlist(lapply(list(c,n,h), is.null)), 0)))

  # loops
  ind = which(elementDf$InLabel==1)
  isotopeDf = expand.grid(sapply(elementDf$number[ind], function(x) c(0:x), simplify = F))
  colnames(isotopeDf) = elementDf$element[ind]
  tempName = as.matrix(apply(isotopeDf[,1:ncol(isotopeDf), drop=F], 1 ,
                             function(x) paste0(x, elementDf$element[ind], collapse = ''))
  )
  isotopeDf$compound = paste0(comp_prefix, "_",
                              apply(tempName[,,drop=F], 1, paste0, collapse="")
  )
  isotopeDf$weight_compound = gln_mass + apply(mapply(`*`, isotopeDf[,c(1:sum(elementDf$InLabel)), drop=F], elementDf$w[ind]), 1, sum)
  isotopeDf$deviance = ppm*0.000001*isotopeDf$weight_compound
  isotopeDf$upper_bound = isotopeDf$weight_compound + isotopeDf$deviance
  isotopeDf$lower_bound = isotopeDf$weight_compound - isotopeDf$deviance
  isotopeDf$isotopeNumbers = apply(isotopeDf[,c(1:sum(elementDf$InLabel)), drop=F], 1, sum)
  for (i in 1:nrow(isotopeDf)){
    d[[isotopeDf$compound[i]]] = list('lb' = isotopeDf$lower_bound[i],
                                      'wc' = isotopeDf$weight_compound[i],
                                      'ub' = isotopeDf$upper_bound[i],
                                      'carbon' = ifelse(elementDf$element[1]==0, 0,
                                                        ifelse(is.null(isotopeDf$C[i]), 0, isotopeDf$C[i])),
                                      'nitrogen' = ifelse(elementDf$element[2]==0, 0,
                                                          ifelse(is.null(isotopeDf$N[i]), 0, isotopeDf$N[i])),
                                      'hydrogen' = ifelse(elementDf$element[3]==0, 0,
                                                          ifelse(is.null(isotopeDf$H[i]), 0, isotopeDf$H[i])),
                                      'isotope_numbers' =  isotopeDf$isotopeNumbers[i])
  }
  return(d)
}
