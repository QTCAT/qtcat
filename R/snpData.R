#' @title Reads a file in table format and creates a 'snpData' object from it.
#' @param file the name of the file which the data are to be read from. If it 
#' does not contain an absolute path, the file name is relative to the current 
#' working directory, \code{getwd()}. Tilde-expansion is performed where 
#' supported.
#' @param sep the field separator character. Values on each line of the file 
#' are separated by this character.
#' @param quote the set of quoting characters. To disable quoting altogether, 
#' use \code{quote = ""}.
#' @param na.strings a character vector of strings which are to be interpreted 
#' as NA values. Is not implemented yet.
#' @param nrows integer: the maximum number of rows to read in.
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
    stop("Line one and tow are of differnt length")
  if (identical(tolower(testRead[[1L]][1L:3L]), c("nam", "chr", "pos"))) {
    rowNames <- TRUE
    snp1 <- testRead[[2L]][-1L:-3L]
  } else if (identical(tolower(testRead[[1L]][1L:2L]), c("chr", "pos"))){
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
    temp$lociNames <- paste0("loci", seq_len(ncol(temp$snpData)))
  }
  out <- new("snpData",
             snpData=temp$snpData,
             alleles=temp$alleles,
             position=temp$position,
             dim=dim(temp$snpData),
             dimnames=list(temp$indivNames, temp$lociNames))
  out
} # read.snpData
