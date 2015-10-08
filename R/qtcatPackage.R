#' @title Quantitative Trait Cluster Association Test
#'
#' @description An association mapping method which jointly analyses all SNPs
#' at once and at the same time accounts for the correlation between them. This
#' makes correction for population structure unnecessary and therefore
#' increases power compared to classical methods like the mixed model.
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
