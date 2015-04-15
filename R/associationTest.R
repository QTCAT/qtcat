
#' geno data
#' @param x snpData
#' @param qtcatClust clusteing
#' @param height height
#' @param max.height max.height
#' @importFrom hit hierarchy
#' @export
qtcatGeno <- function(x, qtcatClust, height, max.height) {
  stopifnot(is(x, "snpData"))
  stopifnot(is(qtcatClust, "qtcatClust"))
  if (is.null(names <- qtcatClust$medoids)) {
    names <- names(qtcatClust)
  }
  x <- as.matrix(x[, colnames(x) %in% names])
  hier <- hierarchy(qtcatClust$dendrogram, height, max.height, colnames(x))
  out <- list(x = x,
              hierarchy = hier,
              clusters = qtcatClust$clusters,
              medoids = qtcatClust$medoids)
  class(out) <- "qtcatGeno"
  out
}

#' peno of xxx
#' @param x data.frame fisrt column individual names secend column pheno  
#' additional columns aditional variables.
#' @export
qtcatPheno <- function(x) {
  if (any(is.na(x))) 
    stop("Missing values are not allowed")
  if(!identical(substring(tolower(colnames(x)[1:2]), 1, 4), c("name", "phen")))
     stop("first column must be 'names', second column must be 'pheno'")
  if (!is.numeric(x[, 2])) 
    stop("phenotype is not numeric")
  if (ncol(x) > 2L) {
    design <- model.matrix(~ . , data=x[, -1:-2])[, -1L, drop=FALSE]
  } else {
    design <- c()
  }
  out <- list(names = as.character(x[, 1]),
              pheno=x[, 2],
              design=design)
  class(out) <- "qtcatPheno"
  out
}

# #'
# qtcatHit <- function() {
#   
# }

