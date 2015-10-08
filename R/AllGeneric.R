#' @title Extract position information from object.
#'
#' @description Extract position information from object.
#'
#' @param object An object.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("getPos", function(object) standardGeneric("getPos"))


#' @title Frequency of alleles in data.
#'
#' @description Frequency of alleles in data.
#'
#' @param x An object.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("alleleFreq", function(x) standardGeneric("alleleFreq"))


#' @title Frequency of heterozygosity in data.
#'
#' @description Frequency of heterozygosity in data.
#'
#' @param x An object.
#' @param dim Dimension over which heterozygosity is estimated.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("hetFreq", function(x, dim=1L) standardGeneric("hetFreq"))
