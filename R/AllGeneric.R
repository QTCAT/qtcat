#' @title Extract position information from object.
#' @docType methods
#' @param object R object.
#' @importFrom methods setGeneric
#' @export
setGeneric("getPos", function(object) standardGeneric("getPos"))

#' @title Frequency of alleles in data.
#' @docType methods
#' @param x R object.
#' @importFrom methods setGeneric
#' @export
setGeneric("alleleFreq", function(x) standardGeneric("alleleFreq"))

#' @title Frequency of heterozygosity in data.
#' @docType methods
#' @param x R object.
#' @param dim Dimension over which heterozygosity is estimated.
#' @importFrom methods setGeneric
#' @export
setGeneric("hetFreq", function(x, dim=1L) standardGeneric("hetFreq"))
