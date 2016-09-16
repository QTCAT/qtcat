#' @title Read SNP tables as a snpMatrix object
#'
#' @description Reads a file in table format and as a \linkS4class{snpMatrix} object.
#'
#' @param file the name of the file which the data are to be read from. If it does not
#' contain an absolute path, the file name is relative to the current working directory,
#' \code{getwd()}. Tilde-expansion is performed where supported.
#' @param sep a field separator character. Values on each line of the file are separated
#' by this character.
#' @param quote the set of quoting characters. To disable quoting altogether, use
#' \code{quote = ""}.
#' @param na.string a string which is  interpreted as NA value.
#' @param nrows a integer, the maximum number of rows to read.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#'
#' @importFrom methods new
#' @export
read.snpData <- function(file, sep = " ",  quote = "\"",
                         na.string = "NA", nrows = -1L) {
  if (!file.exists(file))
    stop("No such file or directory")
  file <- normalizePath(file)
  testRead <- strsplit(readLines(file, n = 2L), sep)
  if (sep != "")
    testRead <- lapply(testRead, function(x, sep) gsub(sep, '',x), sep = sep)
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
  if (na.string !=  "") {
    na.string <- as.character(na.string)
    na.string[is.na(na.string)] <- "NA"
  }
  temp <- read_snpData(file, sep, quote, rowNames, na.string, nrows)
  if (!length(temp$lociNames)) {
    lociNames <- paste0("loci", seq_len(ncol(temp$snpData)))
  } else {
    lociNames <- make.unique(temp$lociNames)
  }
  chr <- suppressWarnings(as.numeric(temp$chr))
  if (any(is.na(chr)))
    chr <- temp$chr
  out <- new("snpMatrix",
             snpData = temp$snpData,
             snpInfo = data.frame(chr = chr,
                                   pos = temp$pos,
                                   allele = t(temp$alleles),
                                   row.names = lociNames),
             dim = dim(temp$snpData),
             dimnames = list(temp$indivNames, lociNames))
  out
}


#' @title A snpMatrix constructor
#'
#' @description Constructs a \linkS4class{snpMatrix} object from the given data.
#'
#' @param x a matrix with individuals in rows and SNPs in columns.
#' @param chr a vector with chromosoms at which SNPs are located.
#' @param pos a vector with genomic positions at which SNPs are located.
#' @param alleleCoding a coding scheme of \code{x} for hom (AA), het (AB), and hom (BB).
#' @param allele.1 labels of allele one, for each SNP.
#' @param allele.2 labels of allele two, for each SNP.
#'
#' @importFrom methods new
#' @importFrom stats na.omit
#' @export
as.snpMatrix <- function(x, chr, pos, alleleCoding = c(-1, 0, 1),
                       allele.1 = NULL, allele.2 = NULL) {
  if (!is.matrix(x))
    x <- as.matrix(x)
  if (is.null(colnames(x)))
    stop("Column names are missing for 'x'")
  if (ncol(x) < 2L)
    stop("'x' has less than two columns")
  if (missing(chr))
    stop("'chr' must be specified")
  if (missing(pos))
    stop("'pos' must be specified")
  if (!is.vector(alleleCoding))
    stop("'alleleCoding' is not a vector")
  x.allele <- na.omit(unique(c(x)))
  if (!all(c(x.allele %in% alleleCoding, alleleCoding %in% x.allele)))
    stop("'alleleCoding' do not match to 'x'")
  if (is.null(allele.1) || is.null(allele.2)) {
    allele.1 <- rep("A", ncol(x))
    allele.2 <- rep("B", ncol(x))
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
    newLabels <- as.raw(c(1, 3))
  } else if (nLabels == 3L) {
    newLabels <- as.raw(c(1, 2, 3))
  }
  y <- matrix(raw(0), nrow(x), ncol(x))
  for (i in 1:nLabels)
    y[alleleCoding[i] == x] <- newLabels[i]
  out <- new("snpMatrix",
             snpData = y,
             snpInfo = data.frame(chr = chr,
                                  pos = pos,
                                  allele.1 = allele.1,
                                  allele.2 = allele.2,
                                  row.names = loci.names),
             dim = dim(y),
             dimnames = list(indiv.names, loci.names))
  out
}


#' @title Subsetting snpMatrix
#'
#' @description Subsetting an object of class \linkS4class{snpMatrix}.
#'
#' @param x  an object of class \linkS4class{snpMatrix}.
#' @param i a indices specifying elements to extract or replace. Indices are booleans,
#' numeric or character vectors.
#' @param j indices specifying elements to extract or replace. Indices are booleans,
#' numeric or character vectors.
#' @param ... not implemented.
#' @param drop not implemented.
#'
#' @importFrom methods setMethod signature new
#' @export
setMethod("[", signature(x = "snpMatrix", i = "ANY", j = "ANY", drop = "missing"),
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
            out <- new("snpMatrix",
                       snpData = snpData,
                       snpInfo = x@snpInfo[j, ],
                       dim = dim(snpData),
                       dimnames = list(rownames(x)[i], colnames(x)[j]))
            out
          }
)


#' @title As matrix method for snpMatrix
#'
#' @description As matrix method for an object of class \linkS4class{snpMatrix}.
#'
#' @param x an object of class \linkS4class{snpMatrix}.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' snpmat <- as.matrix(snp)
#'
#' @importFrom methods setMethod signature
#' @export
setMethod("as.matrix", signature(x = "snpMatrix"),
          function(x) {
            out <- design(x@snpData)
            dimnames(out) <- dimnames(x)
            out
          }
)


#' @title Get position from snpMatrix
#'
#' @description Genetic position info from an object of class
#'  \linkS4class{snpMatrix}.
#'
#' @param object an object of class \linkS4class{snpMatrix}.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' info <- snpInfo(snp)
#'
#' @importFrom methods setMethod signature
#' @export
setMethod("snpInfo", signature(object = "snpMatrix"),
          function(object) {
            out <- object@snpInfo
            if (is.null(out)) {
              cat("No position information available")
            }
            out
          }
)


#' @title Allele frequency
#'
#' @description Allele frequency an object of class \linkS4class{snpMatrix}.
#'
#' @param x an object of class \linkS4class{snpMatrix}.
#' @param maf logical, if true minor allele frequency (default), other ways allele
#' frequency of 'allele.1'.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' af <- alleleFreq(snp)
#'
#' @importFrom methods setMethod signature
#' @export
setMethod("alleleFreq", signature(x = "snpMatrix"),
          function(x, maf = TRUE) {
            out <- afreq(x@snpData, maf)
            names(out) <- colnames(x)
            out
          }
)


#' @title Heterozygosity
#'
#' @description Heterozygosity an object of class \linkS4class{snpMatrix}.
#'
#' @param x an object of class \linkS4class{snpMatrix}.
#' @param dim a integer for dimension. 1 (default) is for rows (individuals), 2 is for
#' columns (SNPs).
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' hf1 <- hetFreq(snp, 1)
#' hf2 <- hetFreq(snp, 2)
#'
#' @importFrom methods setMethod signature
#' @export
setMethod("hetFreq", signature(x = "snpMatrix"),
          function(x, dim = 1) {
            if (dim != 1L && dim != 2L)
              stop("'dim' must be '1' or '2'")
            out <- hetfreq(x@snpData, dim)
            if (dim == 1L)
              names(out) <- rownames(x)
            else if (dim == 2L)
              names(out) <- colnames(x)
            return(out)
          }
)


#' @title NA frequency
#'
#' @description NA frequency in an object of class \linkS4class{snpMatrix}.
#'
#' @param x an object of class \linkS4class{snpMatrix}.
#' @param dim a integer for dimension. 1 (default) is for rows (individuals), 2 is for
#' columns (SNPs).
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' na1 <- naFreq(snp, 1)
#' na2 <- naFreq(snp, 2)
#'
#' @importFrom methods setMethod signature
#' @export
setMethod("naFreq", signature(x = "snpMatrix"),
          function(x, dim = 1) {
            if (dim != 1L && dim != 2L)
              stop("'dim' must be '1' or '2'")
            out <- nafreq(x@snpData, dim)
            if (dim == 1L)
              names(out) <- rownames(x)
            else if (dim == 2L)
              names(out) <- colnames(x)
            return(out)
          }
)
