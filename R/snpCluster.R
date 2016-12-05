#' @title Correlation based distance between SNPs
#'
#' @description This function computes a distance matrix. The distance is estimated
#' as one minus the absolute value of the correlation coefficient \code{1 - abs(cor)}.
#'
#' @param snp an object of class \linkS4class{snpMatrix}.
#'
#' @details See \code{\link[stats]{dist}} for details about the output object.
#' @seealso \code{\link[stats]{dist}}
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#'
#' dist <- distCor(snp[, 1:10])
#'
#' @importFrom methods is
#' @export
distCor <- function(snp) {
  stopifnot(is(snp, "snpMatrix"))
  out <- corDists(snp@snpData)
  attr(out,"Labels") <- colnames(snp)
  attr(out,"Size") <- ncol(snp)
  attr(out,"Diag") <- FALSE
  attr(out,"Upper") <- FALSE
  attr(out,"method") <- "1-abs(cor(snp))"
  attr(out,"call") <- match.call()
  class(out) <- "dist"
  out
}


#' @title Perfect simiarity clusters of SNP
#'
#' @description Finds perfect similarity cluster of SNPs. This is specially usfull in
#' artificial crossing populations.
#'
#' @param snp an object of class \linkS4class{snpMatrix}.
#' @param mc.cores a positive integer for the number of cores for parallel computing. See
#' \code{\link[parallel]{mclapply}} for details.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#'
#' ident <- identicals(snp)
#'
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom stats optimise
#' @export
identicals <- function(snp, mc.cores = 1) {
  stopifnot(is(snp, "snpMatrix"))
  p <- ncol(snp@snpData)
  s <- optimise(function(s, p, m) p * s + p / s * (p / s - 1) / 2 * s / m,
                interval = c(2, p - 1), p = p, m = mc.cores)$minimum
  step <- as.integer(p / (s + 1))
  preclust <- unlist(corPreIdenticals(snp@snpData, step), FALSE)
  kidenticals <- mclapply(preclust, function(i, snp) corIdenticals(snp, i),
                         snp = snp@snpData, mc.cores = mc.cores)
  identclust <- joinCorIdenticals(p, preclust, kidenticals)
  clust <- identclust[[1L]]
  names(clust) <- colnames(snp)
  out <- list(clusters = clust,
              medoids = colnames(snp)[identclust[[2L]] + 1])
  class(out) <- "identicals"
  out
}


#' @title K-medoids clustering of SNPs using randomized search
#'
#' @description Partitioning (clustering) into k clusters "around medoids" by randomized
#' search. \code{1-abs(cor)} is used as distance between SNPs.
#'
#' @param snp an object of class \linkS4class{snpMatrix}.
#' @param k a positive integer specifying the number of clusters, has to be greater than
#' one and less than the number of SNPs.
#' @param maxNeigbours a positive integer specifying the maximum number of randomized
#' searches.
#' @param nLocal a positive integer specifying the number of optimisation runs.
#' @param mc.cores a positive integer for the number of cores for parallel computing. See
#' \code{\link[parallel]{mclapply}} for details.
#'
#' @details The K-medoids clustering is implemented as clustering large applications based
#' upon randomized search (CLARANS) algorithm (Ng and Han 2002). CLARANS is a modification
#' of the partitioning around medoids (PAM) algorithm \code{\link[cluster]{pam}}. Where the
#' PAM algorithm is estimating all distances between SNPs and the respective medoids,
#' CLARANS is searching a random subset of the SNPs. This is independently repeated several
#' times and the result which minimises the average distance the most is reported. This
#' produces results close to those of the PAM algorithm (Ng and Han 2002), though the
#' number of runs and the subset size have to be arbitrarily chosen by the user. The
#' algorithm has two advantages: (i) the number of distance comparisons is dramatically
#' reduced; and (ii) parallelizing is straightforward.
#'
#' @references Ng and J. Han (2002). CLARANS: A method for clustering objects for spatial
#' data mining. \emph{IEEE Transactions on Knowledge and Data Engineering}.
#' \url{http://dx.doi.org/10.1109/TKDE.2002.1033770}).
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#'
#' clust <- clarans(snp, 3)
#'
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
clarans <- function(snp, k, maxNeigbours = 100, nLocal = 10, mc.cores = 1) {
  stopifnot(is(snp, "snpMatrix"))
  if (missing(k))
    stop("'k' must be specifid")
  if (k < 2L)
    stop("'k' must be at least two")
  # cluster optimisation by clarans in parallel
  clarans.i <- function(i, snp, k, maxNeigbours) {
    # cluster optimisation by clarans
    out <- corClarans(snp@snpData, k, maxNeigbours)
    out
  }
  out.nLocal <- mclapply(1L:nLocal, clarans.i,
                         snp, k, maxNeigbours,
                         mc.cores = mc.cores)
  opt.func <- function(i, snp) {snp[[i]][[3L]]}
  all.objectives <- sapply(1:nLocal, opt.func, out.nLocal)
  out.opt <- out.nLocal[[which.min(all.objectives)]]
  clusters <- out.opt[[1L]]
  names(clusters) <- colnames(snp)
  medoids <- out.opt[[2L]] + 1
  names(medoids) <- colnames(snp)[medoids]
  # output
  out <- list(clusters = clusters,
              medoids = medoids,
              objective = out.opt[[3L]],
              all.objectives = all.objectives)
  class(out) <- "k-medoids"
  out
}


#' @title Hierarchical clustering for big SNP data sets.
#'
#' @description A three step approximated hierarchical clustering of SNPs suitable to
#' large data sets.
#'
#' @param snp an object of class \linkS4class{snpMatrix}.
#' @param k a positive integer specifying the number of clusters, less than the number of
#' observations.
#' @param identicals logical, if zero clustering.
#' @param maxNeigbours a positive integer, specifying the maximum number of randomized
#' searches.
#' @param nLocal a positive integer, specifying the number of optimisation runs. Columns
#' have to be similar to \code{snp}.
#' @param method see hclust.
#' @param mc.cores a number of cores for parallel computing. See \code{mclapply} in package
#' parallel for details.
#' @param trace logical, if \code{TRUE} it prints current status of the program.
#' @param ... additional argruments for \code{\link[fastcluster]{hclust}}
#'
#' @seealso clarans
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#'
#' clust <- qtcatClust(snp)
#'
#' @importFrom fastcluster hclust
#' @importFrom stats as.dendrogram
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom utils installed.packages
#' @export
qtcatClust <- function(snp, k, identicals = TRUE, maxNeigbours = 100, nLocal = 10,
                       method = "complete", mc.cores = 1, trace = FALSE, ...) {
  stopifnot(is(snp, "snpMatrix"))
  if (identicals) {
    # identicals
    if (trace)
      cat("Step 1: Search for identicals is running\n")
    identicalFit <- identicals(snp, mc.cores)
    snp <- snp[, identicalFit$medoids]
  } else if (trace) {
    cat("Step 1: Search for identicals is switch off\n")
  }
  # CLARANS
  if (missing(k))
    k <- as.integer(ncol(snp) / 10000L)
  if (k >= 2L) {
    if (identicals && length(identicalFit$medoids) <= k * 2)
      stop("Number of medoids from pefect correlated clustering is < k * 2")
    if (trace)
      cat("Step 2: CLARANS is running, 'k' is:", k, "\n")
    clarFit <- clarans(snp, k, maxNeigbours, nLocal, mc.cores)
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
      stop(paste("Clusters from CLARANS are to big for hclust, choose larger 'k'. The current 'k' is",
                 k))
    # HClust
    if (trace)
      cat("Step 3: HClust is running\n")
    hclust.sub <- function(i, snp, clarFit, method, ...) {
      inx.i <- which(clarFit$clusters == i)
      out <- as.dendrogram(hclust(distCor(snp[, inx.i]), method, ...))
      out
    }
    hclustFit <- mclapply(clust.inx, hclust.sub,
                          snp, clarFit, method, ...,
                          mc.cores = mc.cores)
    dendro <- do.call(merge, c(hclustFit, height = 1, adjust = "add.max"))
  } else {
    if (ncol(snp) > 65536L)
      stop("Data size is to big for hclust, choose larger 'k'")
    # HClust
    if (trace)
      cat("Step 2: CLARANS is switch off\nStep 3: HClust is running\n")
    dendro <- as.dendrogram(hclust(distCor(snp), method, ...))
  }
  if (identicals) {
    out <- list(dendrogram = dendro,
                clusters = identicalFit$clusters,
                medoids = identicalFit$medoids)
  } else {
    medos <- labels(dendro)
    clust <- 1:length(medos)
    names(clust) <- medos
    out <- list(dendrogram = dendro,
                clusters = clust,
                medoids = medos)
  }
  class(out) <- "qtcatClust"
  out
}


#' @title Cut a qtcatClust object
#'
#' @description Cut a qtcatClust object at an specific height.
#'
#' @param snp an object of class \linkS4class{snpMatrix}.
#' @param snpClust an object of class \code{\link{qtcatClust}}.
#' @param absCor a cutting height in absolute value of correlation.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#'
#' cclust <- cutClust(snp, clust, .5)
#'
#' @importFrom methods is
#' @importFrom stats na.omit
#' @export
cutClust <- function(snp, snpClust, absCor = 1) {
  stopifnot(is(snp, "snpMatrix"))
  stopifnot(is(snpClust, "qtcatClust"))
  stopifnot(!missing(absCor))
  dend <- snpClust$dendrogram
  clust <- snpClust$clusters
  if ((1 - absCor) >= attr(dend, "height")) {
    stop("'absCor' outside of range")
  } else {
    cut.dend  <- cut(dend, h = 1 - absCor)
    clust.member <- function(i, dlist, clust) {
      namesclust <- names(clust)[clust %in% unique(clust[labels(dlist[[i]])])]
      clusti <- rep(i, length(namesclust))
      names(clusti) <- namesclust
      return(clusti)
    }
    clust <- unlist(lapply(1:length(cut.dend$lower), clust.member, cut.dend$lower, clust))
    medo <- names(clust)[corMedoids(snp[, names(clust)]@snpData, clust)]
    dend <- rename.leafs(cut.dend$upper, medo)
    clust <- clust[na.omit(match(colnames(snp), names(clust)))]
  }
  out <- list(dendrogram = dend,
              clusters = clust,
              medoids = medo)
  class(out) <- "qtcatClust"
  out
}


#' @title Rename dendrogram leafs
#'
#' @description Rename dendrogram leafs.
#'
#' @param dend a dendrogram.
#' @param labels a vector of new names.
#'
#' @importFrom stats dendrapply is.leaf
#' @keywords internal
rename.leafs <- function(dend, labels){
  dendlabel <- function(n) {
    if(is.leaf(n)) {
      i <<- i + 1L
      attr(n, "label") <- labels[i]
    }
    n
  }
  i <- 0L
  dendrapply(dend, dendlabel)
}
