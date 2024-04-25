#' @title Calculate m/z for given formulas
#'
#' @description This function is calculating the m/z for a chemical formula provided
#'
#' @param comp_formula Formula of the compound
#' @param polarity Polarity of the compound, can be "pos" or "neg" or \code{NULL}. Default to \code{NULL} for a neutral mass.
#'
#' @return The function returns the m/z value calculated
#'
#' @export

get_comp_mass <- function(comp_formula, polarity = NULL){


  #### INNITIALIZATION ####
  w_c = 1.003355
  w_n = 0.997035
  mass_dict = list('N' = 14.003074,
                   'C' = 12.000000,
                   'O' = 15.994915,
                   'H' = 1.0078250,
                   'S' = 31.972072,
                   'P' = 30.973763,
                   'Cl' = 34.968853,
                   'Fe' = 55.934939,
                   'Cu' = 62.929599,
                   'Ca' = 39.962591,
                   'K' = 38.963708,
                   'Br' = 78.91833,
                   'F' = 18.998405,
                   'D' = 2.014102)


  #mass = mass_dict[['H']]

  if(is.null(polarity)){
    mass = 0
  }else if(polarity == "pos")
  {
    mass = 1.007276
  }else if(polarity == "neg")
  {
    mass = - 1.007276
  }

  element_list = get_element_count(comp_formula)
  for(ele in names(element_list)){
    mass = mass + (element_list[[ele]]*mass_dict[[ele]])
  }
  return(mass)
}
