# define
#' @importFrom methods setClassUnion
setClassUnion("matrixOrNULL", c("matrix", "NULL"))


#' @title An S4 class to represent a SNP matrix
#'
#' @description An S4 class to represent a SNP matrix. Storing SNP information
#' by using a byte-level (raw) storage scheme jointly with genomic positions
#' and alleles.
#'
#' @slot snpData A matrix of SNPs stored in type 'raw'. 0 is NA,
#' 1 homozygote AA, 2 and 4 heterozygote phased for AB and BA, and
#' 3 hetrozygoe unphased AB/BA.
#' @slot position A matrix with two rows. The first row contains the
#' chromosomes the second row contains the positions.
#' @slot alleles A matrix of allele labels: A, T, C or, G.
#' @slot dim integer vector with exactly two non-negative values.
#' @slot dimnames list of length two; each component containing NULL or a
#' character vector length equal the corresponding dim element.
#'
#' @importFrom methods setClass
#' @export
setClass("snpData",
         slots = c(snpData = "matrix",
                   position = "matrix",
                   alleles = "matrixOrNULL",
                   dim = "integer",
                   dimnames = "list"))
# TODO: check validity of class
#    - dims of slots
#    - ...
