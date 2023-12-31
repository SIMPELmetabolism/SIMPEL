% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_comp_stage.R
\name{get_comp_stage}
\alias{get_comp_stage}
\title{Identify isotopologues from the XCMS data}
\usage{
get_comp_stage(x, y, comp_lookup_table, r_time, rt_tolerance)
}
\arguments{
\item{x}{this is the m/z from your XCMS_data}

\item{y}{will be RT from your XCMS data}

\item{comp_lookup_table}{This will be your m/z table that was created based on the annotation file}

\item{r_time}{This will be the rt from your annotation file}

\item{rt_tolerance}{will be the retention time window for the upper and lower bounds}
}
\value{
This function returns a list containing compound name, and numbers of carbon, nitrogen, and hydrogen
for identified isotopologue
}
\description{
This function is going to evaluate all of the XCMS_data to identify isotopologues
}
