
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
