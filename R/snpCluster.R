#' @title Correlation based distance among SNPs
#'
#' @description This function computes and returns a distance matrix.
#' Distance is estimated as one minus absolute value of the correlation
#' coefficient 1-abs(cor). The estimation is computed for all pairwise
#' combinations SNPs.
#'
#' @param x An object of class \linkS4class{snpData}.
#' @details See \code{\link[stats]{dist}} for details about the output object.
#' @seealso \code{\link[stats]{dist}}
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' dist <- distCor(snp[, 1:10])
#'
#' @importFrom methods is
#' @export
distCor <- function(x) {
  stopifnot(is(x, "snpData"))
  out <- corDists(x@snpData)
  attr(out,"Labels") <- colnames(x)
  attr(out,"Size") <- ncol(x)
  attr(out,"Diag") <- FALSE
  attr(out,"Upper") <- FALSE
  attr(out,"method") <- "1-abs(cor(x))"
  attr(out,"call") <- match.call()
  class(out) <- "dist"
  out
}


#' @title Perfect simiarity  SNP clusters
#'
#' @description Finds perfect similarity cluster of SNPs. This is specially
#' us full in crossed populations.
#'
#' @param x An object of class \linkS4class{snpData}.
#' @param mc.cores A positive integer for the number of cores for parallel
#' computing. See \code{\link[parallel]{mclapply}} for details.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' ident <- identicals(snp)
#'
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
identicals <- function(x, mc.cores = 1) {
  stopifnot(is(x, "snpData"))
  p <- ncol(x@snpData)
  s <- optimise(function(s, p, m) p * s + p / s * (p / s - 1) / 2 * s / m,
                interval = c(2, p - 1), p = p, m = mc.cores)$minimum
  step <- as.integer(p / (s + 1))
  preclust <- unlist(corPreIdenticals(x@snpData, step), FALSE)
  kidenticals <- mclapply(preclust, function(i, x) corIdenticals(x, i),
                         x = x@snpData, mc.cores = mc.cores)
  identclust <- joinCorIdenticals(ncol(x@snpData), preclust, kidenticals)
  clust <- identclust[[1L]]
  names(clust) <- colnames(x)
  out <- list(clusters = clust,
              medoids = colnames(x)[identclust[[2L]] + 1])
  class(out) <- "identicals"
  out
}


#' @title K-medoids clustering of SNPs using randomized search
#'
#' @description Partitioning (clustering) into k clusters "around medoids" by
#' randomized search. 1-abs(cor) is used as distance among SNPs.
#'
#' @param x An object of class \linkS4class{snpData}.
#' @param k A positive integer specifying the number of clusters, greater than
#' one and less than the number of SNPs.
#' @param maxNeigbours A positive integer specifying the maximum number of
#' randomized searches.
#' @param nLocal A positive integer specifying the number of optimisation runs.
#' @param mc.cores A positive integer for the number of cores for parallel
#' computing. See \code{\link[parallel]{mclapply}} for details.
#'
#' @details The K-medoids clustering is implemented as clustering large
#' applications based upon randomized search (CLARANS) algorithm
#' (Ng and Han 2002). CLARANS is a modification of the partitioning around
#' medoids (PAM) algorithm \code{\link[cluster]{pam}}. Where the PAM
#' algorithm is estimating all distances between SNPs and the
#' respective medoids SNPs, CLARANS is searching a random subset of the SNPs.
#' This is independently repeated several times and the result which minimises
#' the average distance the most is reported. This produces results close to
#' those of the PAM algorithm (Ng and Han 2002), even though the number of runs
#' and the subset size have to be arbitrarily chosen by the user. The algorithm
#' has two advantages: (i) the number of distance comparisons is dramatically
#' reduced; and (ii) parallelizing is straightforward.
#'
#' @references Ng and J. Han (2002). CLARANS: A method for clustering objects
#' for spatial data mining. \emph{IEEE Transactions on Knowledge and Data
#' Engineering}. \url{http://dx.doi.org/10.1109/TKDE.2002.1033770}).
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- clarans(snp, 3)
#'
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
clarans <- function(x, k, maxNeigbours = 100, nLocal = 10, mc.cores = 1) {
  stopifnot(is(x, "snpData"))
  if (missing(k))
    stop("'k' must be specifid")
  if (k < 2L)
    stop("'k' must be at least two")
  # cluster optimisation by clarans in parallel
  clarans.i <- function(i, x, k, maxNeigbours) {
    # cluster optimisation by clarans
    out <- corClarans(x@snpData, k, maxNeigbours)
    out
  } # clarans.i
  out.nLocal <- mclapply(1L:nLocal, clarans.i,
                            x, k, maxNeigbours,
                            mc.cores = mc.cores)
  opt.func <- function(i, x) {x[[i]][[3L]]}
  all.objectives <- sapply(1:nLocal, opt.func, out.nLocal)
  out.opt <- out.nLocal[[which.min(all.objectives)]]
  clusters <- out.opt[[1L]]
  names(clusters) <- colnames(x)
  medoids <- out.opt[[2L]] + 1
  names(medoids) <- colnames(x)[medoids]
  # output
  out <- list(clusters = clusters,
              medoids = medoids,
              objective = out.opt[[3L]],
              all.objectives = all.objectives)
  class(out) <- "k-medoids"
  out
}


#' @title A three step approximated hierarchical clustering of SNPs
#'
#' @description Hierarchical clustering for big SNP data sets.
#'
#' @param x A object of class \linkS4class{snpData}.
#' @param k A positive integer specifying the number of clusters, less than
#' the number of observations.
#' @param identicals Logical, if zero clustering.
#' @param maxNeigbours Positive integer, specifying the maximum number of
#' randomized searches.
#' @param nLocal Positive integer, specifying the number of optimisation runs.
#' Columns have to be similar to \code{x}.
#' @param method See hclust.
#' @param mc.cores Number of cores for parallel computing. See \code{mclapply}
#' in package parallel for details.
#' @param trace If \code{TRUE} it prints current status of the program.
#' @param ... additional argruments for \code{\link[stats]{hclust}}
#'
#' @seealso clarans
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#'
#' @importFrom parallel mclapply
#' @importFrom stats hclust
#' @importFrom methods is
#' @export
qtcatClust <- function(x, k, identicals = TRUE,
                       maxNeigbours = 100, nLocal = 10,
                       method = "complete", mc.cores = 1, trace = FALSE, ...) {
  if (is.element("fastcluster", rownames(installed.packages()))) {
    hclust <- fastcluster::hclust
  }
  stopifnot(is(x, "snpData"))
  if (identicals) {
    # identicals
    if (trace)
      cat("Step 1: Search for identicals is running\n")
    identicalFit <- identicals(x, mc.cores)
    x <- x[, identicalFit$medoids]
  } else if (trace) {
    cat("Step 1: Search for identicals is switch off\n")
  }
  # CLARANS
  if (missing(k))
    k <- as.integer(ncol(x) / 10000L)
  if (k >= 2L) {
    if (identicals && length(identicalFit$medoids) <= k * 2)
      stop("Number of medoids from pefect correlated clustering is < k * 2")
    if (trace)
      cat("Step 2: CLARANS is running, 'k' is:", k, "\n")
    clarFit <- clarans(x, k, maxNeigbours, nLocal, mc.cores)
    if (trace)
      cat("   objectives:",
          format(clarFit$all.objectives, sort = TRUE, digits = 4), "\n")
    # if cluster < 2 add to a bigger cluster
    clust.inx <- seq_len(k)
    cluster.size <- rep(NA, k)
    for (i in clust.inx)
      cluster.size[i] <- sum(clarFit$clusters == i)
    smallclust  <- which(cluster.size < 2)
    if (length(smallclust)) {
      clust.inx <- clust.inx[-smallclust]
      min.bigclust <- clust.inx[which.min(min(cluster.size[clust.inx]))]
      clarFit$clusters[clarFit$clusters %in% smallclust] <- min.bigclust
    }
    if (max(cluster.size) > 65536L)
      stop("Clusters from CLARANS are to big for hclust, choose larger 'k'")
    # HClust
    if (trace)
      cat("Step 3: HClust is running\n")
    hclust.sub <- function(i, x, clarFit, method, ...) {
      inx.i <- which(clarFit$clusters == i)
      out <- as.dendrogram(hclust(distCor(x[, inx.i]), method, ...))
      out
    } # hclust.sub
    hclustFit <- mclapply(clust.inx, hclust.sub,
                          x, clarFit, method, ...,
                          mc.cores = mc.cores)
    dendro <- do.call(merge, c(hclustFit, height = 1, adjust = "add.max"))
  } else {
    if (ncol(x) > 65536L)
      stop("Data size is to big for hclust, choose larger 'k'")
    # HClust
    if (trace)
      cat("Step 2: CLARANS is switch off\nStep 3: HClust is running\n")
    dendro <- as.dendrogram(hclust(distCor(x), method, ...))
  }
  if (identicals) {
    out <- list(dendrogram = dendro,
                clusters = identicalFit$clusters,
                medoids = identicalFit$medoids)
  } else {
    out <- list(dendrogram = dendro,
                clusters = NULL,
                medoids = NULL)
  }
  class(out) <- "qtcatClust"
  out
}


#' @title Cut a qtcatClust object
#'
#' @description Cut a qtcatClust object at an specific height.
#'
#' @param x qtcatClust object.
#' @param absCor cutting height in absolute value of correlation.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' cclust <- cutClust(clust, .5)
#'
#' @importFrom methods is
#' @export
cutClust <- function(x, absCor) {
  stopifnot(is(x, "qtcatClust"))
  stopifnot(!missing(absCor))
  clust.dend <- cutDend(x$dendrogram, absCor)
  if (!is.null(x$clusters)) {
    inx.names <- match(names(clust.dend), names(x$clusters))
    inx.list <- list()
    clust.list <- list()
    for (i in unique(clust.dend)) {
      inx.i <- inx.names[clust.dend == i]
      inx.list[[i]] <- which(x$clusters %in% x$clusters[inx.i])
      clust.list[[i]] <- rep(i, length(inx.list[[i]]))
    }
    out <- unlist(clust.list)[order(unlist(inx.list))]
    names(out) <- names(x$clusters)
  } else {
    out <- clust.dend
  }
  out
}


#' @title Cut a dendrogram
#'
#' @description Cut a dendrogram at a spcific hight.
#'
#' @param x A dendrogram.
#' @param absCor Cutting height in absolute value of correlation.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' cdend <- cutDend(clust$dendrogram, .5)
#'
#' @importFrom methods is
#' @export
cutDend <- function(x, absCor) {
  stopifnot(is(x, "dendrogram"))
  stopifnot(!missing(absCor))
  h <- 1 - absCor
  if (h >= attr(x, "height")) {
    names.clust <- labels(x)
    cluster <- rep(1, length(names.clust))
    names(cluster) <- names.clust
  } else {
    cut.x <- cut(x, h = h)$lower
    clust.member <- function(i, x) {
      names.clust <- labels(x[[i]])
      clust <- rep(i, length(names.clust))
      names(clust) <- names.clust
      return(clust)
    }
    cluster <- lapply(1:length(cut.x), clust.member, cut.x)
  }
  unlist(cluster)
}


#' @title Find the medoid
#'
#' @description Find the medoid of each cluster.
#'
#' @param x An object of class \linkS4class{snpData}.
#' @param clusters Vector of cluster groups \code{cutClust}.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' cclust <- cutClust(clust, .5)
#' mclust <- medoids(snp, cclust)
#'
#' @importFrom methods is
#' @export
medoids <- function(x, clusters) {
  stopifnot(is(x, "snpData"))
  medoids <- corMedoids(x@snpData, clusters)
  out <- rbind(clusters = clusters, medoids = medoids)
  class(out) <- "medoids"
  out
}
