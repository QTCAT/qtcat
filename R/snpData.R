
#' @title read snp data as snpData object
#' @param file the name of the file which the data are to be read from.
#' @param sep the field separator character. Values on each line of the file 
#' are separated by this character.
#' @param na.strings a character vector of strings which are to be interpreted 
#' as NA values. Is not implemented yet.
#' @param nrows integer: the maximum number of rows to read in.
#' @importFrom methods new
#' @export
read.snpData <- function(file, sep=" ",  na.strings, nrows) {
  if (!missing(na.strings))
    stop("Missing data are not allowed in qtcat")
  if (missing(nrows))
    nrows <- -1L
  test <- read.table(file, header=TRUE, sep=sep, nrows=10, 
                     stringsAsFactors=FALSE)
  nam <- tolower(substring(colnames(test[1L:3L]), 1L, 3L))
  if (identical(nam[1L:2L], c("chr", "pos"))){
    rowNames <- FALSE
  } else if (identical(nam, c("nam", "chr", "pos"))) {
    rowNames <- TRUE
  } else {
    stop("The first two columns have to be either 'chr', and 'pos' or the first three columns have to be 'names', 'chr', and 'pos'")
  }
  if (any(nchar(unique(unlist(test[, -1:-3]))) != 2L))
      stop("allele labels have not two characters")
  # read data
  temp <- read_snpData(file, sep, '"', rowNames, "N", nrows)
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
