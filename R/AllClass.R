
# define
#' @importFrom methods setClassUnion
setClassUnion("matrixOrNULL", members=c("matrix", "NULL"))

#' @title Class "snpData" for a SNP matrix (S4 class)
#' @name snpData-class
#' @docType class
#' @description A S4 class for a SNP matrix. Storing SNP information using a 
#' byte-level (raw) storage scheme jointly with positions.
#' @slot snpData A matrix of SNPs stored in type 'raw'.
#' @slot position A matrix with two rows. The first row contains the 
#' chromosomes the second row contains the positions.
#' @slot alleles A matrix of allele labels: A, T, C or, G.
#' @slot Dim Dimentions of snpData
#' @slot Dimnames Row and column names of snpData
#' @importFrom methods setClass
#' @export
setClass("snpData",
         slots=c(snpData = "matrix",
                 position = "matrix",
                 alleles = "matrixOrNULL",
                 dim = "integer", 
                 dimnames = "list"))
