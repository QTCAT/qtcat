#' @title Extract position and allele information from object.
#'
#' @description Extract position and allele information from object.
#'
#' @param object An object.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("snpInfo", function(object) standardGeneric("snpInfo"))


#' @title Frequency of alleles in data.
#'
#' @description Frequency of alleles in data.
#'
#' @param x An object.
#' @param maf If true minor allele frequency (default), other ways allele frequency of 
#' 'allele.1'.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("alleleFreq", function(x, maf = TRUE) standardGeneric("alleleFreq"))


#' @title Frequency of heterozygosity in data.
#'
#' @description Frequency of heterozygosity in data.
#'
#' @param x An object.
#' @param dim Dimension over which heterozygosity is estimated.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("hetFreq", function(x, dim = 1L) standardGeneric("hetFreq"))


#' @title Frequency of NAs in data.
#'
#' @description Frequency of NAs in data.
#'
#' @param x An object.
#' @param dim Dimension over which heterozygosity is estimated.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("naFreq", function(x, dim=1L) standardGeneric("naFreq"))
