#' IONiseR: A package for assessing quality of MinION data
#' 
#' IONiseR provides tools for the quality assessment of Oxford Nanopore MinION 
#' data. It extracts summary statistics from a set of fast5 files and can be 
#' used either before or after base calling.  In addition to standard summaries
#' of the read-types produced, it provides a number of plots for visualising 
#' metrics relative to experiment run time or spatially over the surface of a 
#' flowcell.
#'
#' @docType package
#' @name IONiseR
#' 
#' @import ggplot2
#' @import rhdf5
#' 
#' @importFrom methods setClass setGeneric setMethod
#' @importFrom magrittr "%>%"
#' @importFrom ShortRead strand
#' 
#' @importClassesFrom ShortRead ShortReadQ
NULL
