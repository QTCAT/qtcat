#' geno data
#' @param x snpData
#' @param clusters clusteing
#' @param height height
#' @param max.height max.height
#' @importFrom hit hierarchy
#' @export
qtcatGeno <- function(x, clusters, height, max.height=.3) {
  stopifnot(is(x, "snpData"))
  stopifnot(is(clusters, "qtcatClust"))
  if (is.null(names <- clusters$medoids)) {
    names <- names(clusters)
  }
  # TODO: chack allele.freq
  desMat <- as.matrix(x[, colnames(x) %in% names])
  hier <- hierarchy(clusters$dendrogram, height, max.height, colnames(x))
  out <- list(x = desMat,
              hierarchy = hier,
              clusters = clusters$clusters,
              medoids = clusters$medoids,
              positions = getPos(x))
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

#' hit
#' @param pheno qtcatPheno
#' @param geno qtcatGeno
#' @param B number of sample-splits.
#' @param p.samp1 fraction of data used for the LASSO. The ANOVA uses 
#' \code{1 - p.samp1}.
#' @param gamma vector of gamma-values.
#' @param max.p.esti maximum alpha level. All p-values above this value are set 
#' to one. Small \code{max.p.esti} values reduce computing time.
#' @param mc.cores number of cores for parallelising. Theoretical maximum is 
#' 'B'. For details see \code{\link[parallel]{mclapply}}.
#' @param trace if TRUE it prints current status of the program.
#' @param ... additional arguments for \code{\link[glmnet]{cv.glmnet}}.
#' @export
qtcatHit <- function(pheno, geno, 
                     B=50, p.samp1=0.5, gamma=seq(0.05, 0.99, by=0.01), 
                     max.p.esti=1, mc.cores=1L, trace=FALSE, ...) {
  stopifnot(is(pheno, "qtcatPheno"))
  stopifnot(is(geno, "qtcatGeno"))
  id <- intersect(pheno$names, rownames(geno$x))
  if (!length(id))
    stop("The ID intersect of 'pheno' and 'geno' is emty")
  if (length(id.uniqueGeno <- setdiff(rownames(geno$x), id)))
    warning("The following individuals are part of 'geno' but not of 'pheno':\n", 
            paste(id.uniqueGeno, collapse=" "))
  if (length(id.uniquePheno <- setdiff(pheno$names, id)))
    warning("The following individuals are part of 'pheno' but not of 'geno':\n", 
            paste(id.uniquePheno, collapse=" "))
  phenoInx <- match(pheno$names, id, nomatch=0)
  genoInx <- match(pheno$names[phenoInx], rownames(geno$x))
  x <- cbind(geno$x[genoInx, ], pheno$design[phenoInx, ])
  y <- pheno$pheno[phenoInx]
  fitHit <- hit(x, y, geno$hierarchy, 
                B, p.samp1, gamma, max.p.esti, mc.cores, trace, ...)
  out <- c(fitHit,
           geno[3:5])
  class(out) <- c("hit", "qtcatHit")
  out
}

#' @title Include zero clusters x
#' @description Include zero clusters
#' @param x hit object.
#' @param clusters object of class qtcatClust.
#' @param alpha alpha level.
#' @param max.height max. height to consider.
#' @importFrom hit hit
#' @export
qtcatSigClust <- function(x, clusters, alpha = 0.05, max.height) {
  stopifnot(is(clusters, "qtcatClust"))
  y <- summary(x, alpha, max.height)
  signames <- rownames(y)
  sigclust <- match(signames, clusters$medoids)
  sigClust <- matrix(0, length(clusters$clusters), 3L)
  sigClust[, 2] <- 1
  for (i in seq_along(signames)) {
    inx <- which(clusters$clusters == sigclust[i])
    sigClust[inx, ] <- matrix(unlist(y[signames[i], ]), length(inx), 3, byrow=TRUE)
  }
  rownames(sigClust) <- names(clusters$clusters)
  colnames(sigClust) <- colnames(y)
  sigClust <- as.data.frame(sigClust[sigClust[, 1] != 0, ,drop=FALSE])
  sigClust <- cbind(t(getPos(x)[, rownames(sigClust), drop=FALSE]), sigClust)
  sigClust
} # sigCluster
