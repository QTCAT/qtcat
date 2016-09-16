#' @title Quantitative Trait Cluster Association Test
#'
#' @description QTCAT tests all SNPs jointly in ther association to the phenotype and at
#' the same time considers correlation between them. Thus, correction for population
#' structure becomes unnecessary, which in many cases results in a power advantages
#' compared to other methods.
#'
#' @author Jonas Klasen
#'
#' @useDynLib qtcat
#' @importFrom Rcpp evalCpp
#'
#' @docType package
#' @name qtcat-package
#' @aliases qtcat
NULL
