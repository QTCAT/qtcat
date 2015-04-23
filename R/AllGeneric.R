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
setGeneric("allele.freq", function(x) standardGeneric("allele.freq"))

#' @title Frequency of heterozygosity in data.
#' @docType methods
#' @param x R object.
#' @param dim Dimension over which heterozygosity is estimated.
#' @importFrom methods setGeneric
#' @export
setGeneric("het.freq", function(x, dim=1L) standardGeneric("het.freq"))
