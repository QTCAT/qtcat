#' @title A S4 class to represent a SNP-matrix
#'
#' @description A S4 class to represent a SNP matrix. Storing SNP information, by using a
#' byte-level (raw) storage scheme, jointly with genomic position and allele information.
#'
#' @slot snpData a matrix of SNPs stored in type 'raw'. 00 is NA, 01 homozygote AA, 02
#' heterozygote AB, and 03 homozygote BB.
#' @slot snpInfo data.frame with four columns. The first col. contains the chromosomes,
#' the second col. the positions, the third col. the first allele and the fourth the second
#' allele.
#' @slot dim an integer vector with exactly two non-negative values.
#' @slot dimnames a list of length two; each component containing NULL or a character vector
#' length equal the corresponding dim element.
#'
#' @importFrom methods setClass
#' @export
setClass("snpMatrix",
         slots = c(snpData = "matrix",
                   snpInfo = "data.frame",
                   dim = "integer",
                   dimnames = "list"))
