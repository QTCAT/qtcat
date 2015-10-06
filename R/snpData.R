#' @title Read SNP tables to snpData object
#' @description Reads a file in table format and returns a 'snpData' object.
#' @param file The name of the file which the data are to be read from. If it
#' does not contain an absolute path, the file name is relative to the current
#' working directory, \code{getwd()}. Tilde-expansion is performed where
#' supported.
#' @param sep The field separator character. Values on each line of the file
#' are separated by this character.
#' @param quote The set of quoting characters. To disable quoting altogether,
#' use \code{quote = ""}.
#' @param na.strings A character vector of strings which are to be interpreted
#' as NA values. Is not implemented yet.
#' @param nrows Integer, the maximum number of rows to read in.
#' @importFrom methods new
#' @export
read.snpData <- function(file, sep = " ",  quote = "\"",
                         na.strings, nrows = -1L) {
  # checks
  if (!file.exists(file))
    stop("No such file or directory")
  file <- normalizePath(file)
  if (!missing(na.strings))
    stop("'qtcat' dos not allowed missing values in the SNP-matrix")
  testRead <- strsplit(readLines(file, n = 2L), sep)
  if (length(testRead[[1L]]) <= 3L)
    stop("In line one the separator character 'sep' doesn't exist")
  if (length(testRead[[1L]]) != length(testRead[[2L]]))
    stop("Line one and two are of differnt length")
  if (!identical(unique(testRead[[1L]]), testRead[[1L]]))
    stop("Column names are not unique")
  firstCols <- tolower(substring(testRead[[1L]][1L:3L], 1L, 3L))
  if (identical(firstCols, c("nam", "chr", "pos"))) {
    rowNames <- TRUE
    snp1 <- testRead[[2L]][-1L:-3L]
  } else if (identical(firstCols[1L:2L], c("chr", "pos"))) {
    rowNames <- FALSE
    snp1 <- testRead[[2L]][-1L:-2L]
  } else {
    stop("Either the first three columns contain 'names', 'chr', and 'pos' or alternatively the first two columns contain 'chr', and 'pos'")
  }
  if (any(nchar(snp1) != 2L))
    stop("Every position in the SNP-matrix has to be specified by two characters, missing values are not allowed")
  # read data from file
  temp <- read_snpData(file, sep, quote, rowNames, "ZZZ", nrows)
  # make snpData object
  if (identical(temp$lociNames, character(0))) {
    lociNames <- paste0("loci", seq_len(ncol(temp$snpData)))
  } else {
    lociNames <- make.unique(temp$lociNames)
  }
  out <- new("snpData",
             snpData = temp$snpData,
             position = temp$position,
             alleles = temp$alleles,
             dim = dim(temp$snpData),
             dimnames = list(temp$indivNames, lociNames))
  out
} # read.snpData

#' @title converts data to  snpData object
#' @param x matrix with indoviduals in rows and SNPs in columns.
#' @param position matrix with chromosome and position in rows and
#' SNPs in columns.
#' @param alleleCoding coding of \code{x} for hom het hom.
#' @param alleles  labels of alleles, for each SNP.
#' @importFrom methods new
#' @export
as.snpData <- function(x, position, alleleCoding = c(-1, 0, 1),
                       alleles = NULL) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (is.null(colnames(x))) {
    stop("Column names are missing for 'x'")
  }
  if (ncol(x) < 2L) {
    stop("'x' has less than two columns")
  }
  if (missing(position)) {
    stop("'position' must be specified")
  } else {
    if (nrow(position) != 2L) {
      stop("'position' must have two rows")
    }
    if (is.null(colnames(position))) {
      stop("Column names are missing for 'position'")
    }
    if (any(colnames(x) != colnames(position))) {
      stop("Column names of 'x' and 'position' differ")
    }
  }
  if (!is.vector(alleleCoding)) {
    stop("'alleleCoding' is not a vector")
  }
  x.allele <- na.omit(unique(c(x)))
  if (!all(c(x.allele %in% alleleCoding, alleleCoding %in% x.allele))) {
    stop("'alleleCoding' do not match to 'x'")
  }
  if (!is.null(alleles)) {
    if (nrow(alleles) != 2L) {
      stop("'alleles' must have two rows")
    }
    if (is.null(colnames(alleles))) {
      stop("Column names are missing for 'alleles'")
    }
    if (any(colnames(x) != colnames(alleles))) {
      stop("Column names of 'x' and 'alleles' differ")
    }
    attr(alleles, 'dimnames') <- NULL
  }
  if (is.null(rownames(x))) {
    indiv.names <- paste0("indiv", seq_len(nrow(x)))
  } else {
    indiv.names <- rownames(x)
  }
  loci.names <- colnames(x)
  attr(x, 'dimnames') <- NULL
  attr(position, 'dimnames') <- NULL
  nLabels <- length(alleleCoding)
  if (nLabels == 2L) {
    newLabels <- as.raw(c(1, 5))
  } else if (nLabels == 3L) {
    newLabels <- as.raw(c(1, 3, 5))
  } else if (nLabels == 4L) {
    newLabels <- as.raw(c(1, 2, 4, 5))
  }
  y <- matrix(raw(0), nrow(x), ncol(x))
  for (i in 1:nLabels) {
    y[alleleCoding[i] == x] <- newLabels[i]
  }
  out <- new("snpData",
             snpData = y,
             position = position,
             alleles = alleles,
             dim = dim(y),
             dimnames = list(indiv.names, loci.names))
  out
}

#' @title Sub snpData
#' @docType methods
#' @param x snpData object.
#' @param i Indices specifying elements to extract or replace. Indices are
#' numeric or character vectors.
#' @param j indices specifying elements to extract or replace. Indices are
#' numeric or character vectors.
#' @param ... Not implemented.
#' @param drop Not implemented.
#' @importFrom methods setMethod signature new
#' @export
setMethod("[", signature(x = "snpData", i = "ANY", j = "ANY", drop = "missing"),
          function(x, i, j, ..., drop) {
            if (!missing(i)) {
              if (is.character(i)) {
                i <- match(i, rownames(x))
              }
            }
            if (!missing(j)) {
              if (is.character(j)) {
                j <- match(j, colnames(x))
              }
            }
            snpData <- x@snpData[i, j, drop = FALSE]
            out <- new("snpData",
                       snpData = snpData,
                       position = x@position[, j, drop = FALSE],
                       alleles = if (is.null(x@alleles)) NULL
                                 else x@alleles[, j, drop = FALSE],
                       dim = dim(snpData),
                       dimnames = list(rownames(x)[i], colnames(x)[j]))
            out
          }
) # `[`

#' @title snpData as matrix
#' @docType methods
#' @param x snpData object.
#' @param inx1 ...
#' @param inx2 ...
#' @param ... Not implemented.
#' @importFrom methods setMethod signature
#' @export
setMethod("as.matrix", signature(x = "snpData"),
          function(x, inx1 = 1:ncol(x), inx2 = NULL, ...) {
            stopifnot(is(x, "snpData"))
            if (is.null(inx2)) {
              out <- design(x@snpData, inx1 - 1, inx1 - 1)
              colnames(out) <- colnames(x)[inx1]
            } else {
              out <- design(x@snpData, inx1 - 1, inx2 - 1)
              nam.dat <- data.frame(colnames(x)[inx1], colnames(x)[inx2],
                                    stringsAsFactors = FALSE)
              colnames(out) <- apply(nam.dat, 1, function(x) {
                if (x[1] == x[2]) {
                  return(x[1])
                } else {
                  return(paste(x, collapse = ":"))
                }
              }) # apply nam.dat
            } # if else is.null(inx2)
            rownames(out) <- rownames(x)
            out
          }
) # as.matrix

#' @title Get position from snpData
#' @docType methods
#' @param object snpData object.
#' @importFrom methods setMethod signature
#' @export
setMethod("getPos", signature(object = "snpData"),
          function(object) {
            out <- object@position
            if (!is.null(out)) {
              colnames(out) <- colnames(object)
              rownames(out) <- c("chr", "pos")
            } else {
              cat("No position information available")
            }
            out
          }
) # getPos

#' @title Allele frequency
#' @docType methods
#' @param x snpData object.
#' @importFrom methods setMethod signature
#' @export
setMethod("alleleFreq", signature(x = "snpData"),
          function(x) {
            out <- freqs2(x@snpData)[[1]]
            names(out) <- colnames(x)
            out
          }
) # alleleFreq

#' @title Heterozygosity
#' @docType methods
#' @param x snpData object.
#' @param dim Interger for dimension.
#' @importFrom methods setMethod signature
#' @export
setMethod("hetFreq", signature(x = "snpData"),
          function(x, dim = c(1, 2)) {
            if (dim[1] == 2)
              return(freqs2(x@snpData)[[2]])
            if (dim[1] == 1)
              return(freq1(x@snpData))
            return(NULL)
          }
) # "hetFreq"
