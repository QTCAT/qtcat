#' @title Extract genomic position and allele information.
#'
#' @description Extract genomic position and allele information from object.
#'
#' @param object an object, for which a corresponding method exists.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("snpInfo", function(object) standardGeneric("snpInfo"))


#' @title  Allele Frequency.
#'
#' @description Frequency of alleles in data set.
#'
#' @param x an object, for which a corresponding method exists.
#' @param maf logical, if true minor allele frequency (default), other ways allele
#' frequency of 'allele.1'.
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("alleleFreq", function(x, maf = TRUE) standardGeneric("alleleFreq"))


#' @title Heterozygosity Frequency.
#'
#' @description Frequency of heterozygosity in data set.
#'
#' @param x an object, for which a corresponding method exists.
#' @param dim dimension over which heterozygosity is estimated. 1 (default) is for rows
#' (individuals), 2 is for columns (SNPs).
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("hetFreq", function(x, dim = 1L) standardGeneric("hetFreq"))


#' @title Missing Data Frequency.
#'
#' @description Frequency of missing Data in data set.
#'
#' @param x an object, for which a corresponding method exists.
#' @param dim dimension over which heterozygosity is estimated. 1 (default) is for rows
#' (individuals), 2 is for columns (SNPs).
#'
#' @importFrom methods setGeneric
#' @export
setGeneric("naFreq", function(x, dim = 1L) standardGeneric("naFreq"))
