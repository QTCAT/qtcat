# define
#' @importFrom methods setClassUnion
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

#' @title Class "snpData" for a SNP matrix (S4 class)
#' @name snpData-class
#' @docType class
#' @description A S4 class for a SNP matrix. Storing SNP information using a 
#' byte-level (raw) storage scheme jointly with positions.
#' @slot snpData A matrix of SNPs stored in type 'raw'.
#' @slot position A matrix with two rows. The first row contains the 
#' chromosomes the second row contains the positions.
#' @slot alleles A matrix of allele labels: A, T, C or, G.
#' @slot dim integer vector with exactly two non-negative values.
#' @slot dimnames list of length two; each component containing NULL or a 
#' character vector length equal the corresponding Dim element.
#' @importFrom methods setClass
#' @export
setClass("snpData",
         slots=c(snpData = "matrix",
                 position = "matrix",
                 alleles = "matrixOrNULL",
                 dim = "integer", 
                 dimnames = "list"))
# TODO: validity check of class 
#    - names individuals and snps unique?
#    - no NAs
#    - ...