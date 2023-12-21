#' @title Identify isotopologues from the XCMS data
#'
#' @description This function is going to evaluate all of the XCMS_data to identify isotopologues
#'
#' @param x this is the m/z from your XCMS_data
#' @param y will be RT from your XCMS data
#' @param comp_lookup_table This will be your m/z table that was created based on the annotation file
#' @param r_time This will be the rt from your annotation file
#' @param rt_tolerance will be the retention time window for the upper and lower bounds
#'
#' @return This function returns a list containing compound name, and numbers of carbon, nitrogen, and hydrogen
#' for identified isotopologue
#'
#' @export

get_comp_stage <- function(x, y, comp_lookup_table, r_time, rt_tolerance)
{
  d = list()
  for(k in names(comp_lookup_table)){
    lb = comp_lookup_table[[k]]['lb']
    ub = comp_lookup_table[[k]]['ub']
    c = comp_lookup_table[[k]]['carbon']
    n = comp_lookup_table[[k]]['nitrogen']
    h = comp_lookup_table[[k]]['hydrogen']
    #for a given compound in the lookup table
    #determine if its characteristics fall within the lower and upper
    #bounds of the m/z, rt for features that are actually in the XCMS_data
    if (lb <= x & x <= ub){
      if((r_time - rt_tolerance) <= y & y <= (r_time + rt_tolerance)){

        d = list('compound' = k,
                 'carbon' = c,
                 'nitrogen' = n,
                 'hydrogen' = h)
        return(d)
      }
    }
  }
}
