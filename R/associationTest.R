#' @title A genotype object constructor
#'
#' @description Constructs an S3 object containing SNP matrix and SNP
#' hierarchy. This is needed for \code{\link{qtcatHit}} as input.
#'
#' @param x An object of S4 class \linkS4class{snpData}.
#' @param clusters An object of class \code{\link{qtcatClust}}.
#' @param absCor Vector of absolute value of correlations considered in the
#' hierarchy.
#' @param min.absCor Minimum absolute value of correlation considered. A value
#' in the range from 0 to 1.
#'
#' @examples
#' # file containing example data for SNP data
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' # Construct geotype object
#' geno <- qtcatGeno(snp, clust)
#'
#' @importFrom hit as.hierarchy
#' @importFrom methods is
#' @export
qtcatGeno <- function(x, clusters, absCor, min.absCor=.7) {
  stopifnot(is(x, "snpData"))
  stopifnot(is(clusters, "qtcatClust"))
  if (!setequal(names(clusters$clusters), colnames(x)))
    stop("SNP names of 'x' and 'clusters' differ")
  if (is.null(names <- clusters$medoids))
    names <- names(clusters)
  # TODO: chack alleleFreq?
  desMat <- as.matrix(x[, colnames(x) %in% names])
  if (missing(absCor))
    hier <- as.hierarchy(clusters$dendrogram, max.height = 1 - min.absCor,
                      names = colnames(desMat))
  else
    hier <- as.hierarchy(clusters$dendrogram, 1 - absCor, 1 - min.absCor,
                      colnames(desMat))
  out <- list(x = desMat,
              hierarchy = hier,
              clusters = clusters$clusters,
              medoids = clusters$medoids,
              positions = getPos(x))
  class(out) <- "qtcatGeno"
  out
}


#' @title A phenotype object constructor
#'
#' @description Constructs an S3 object containing phenotype and if
#' additional covariats exist a design matrix of those. This is needed for
#' \code{\link{qtcatHit}} as input.
#'
#' @param x A data frame, containing in the first column individual names, in
#' the second column phenotypic observations and eventuality additional
#' variables in additional columns.
#'
#' @examples
#' # file containing example data for a phenotype.
#' pfile <- system.file("extdata/phenodata.csv", package = "qtcat")
#' pdat <- read.csv(pfile, header = TRUE)
#' # Construct phenotype object
#' pheno <- qtcatPheno(pdat)
#'
#' @export
qtcatPheno <- function(x) {
  if (any(is.na(x)))
    stop("Missing values are not allowed")
  if (!identical(substring(tolower(colnames(x)[1L:2L]), 1L, 4L),
                c("name", "phen")))
     stop("first column must be 'names', second column must be 'pheno'")
  if (!is.numeric(x[, 2L]))
    stop("phenotype is not numeric")
  if (ncol(x) > 2L) {
    xdat <- x[, -(1L:2L), drop = FALSE]
    design <- model.matrix(~., data = xdat)[, -1L, drop = FALSE]
  } else {
    design <- matrix(nrow = nrow(x), ncol = 0L)
  }
  out <- list(names = as.character(x[, 1L]),
              pheno = x[, 2L],
              design = design)
  class(out) <- "qtcatPheno"
  out
}


#' @title Fitting Hierarchical Inference Testing
#'
#' @description Hierarchical inference testing for phenotype-SNP association.
#'
#' @param pheno An object of class \code{\link{qtcatPheno}}.
#' @param geno An object of class \code{\link{qtcatGeno}}.
#' @param B A integer indicating the number of sample-splits.
#' @param p.samp1 A value specifying the fraction of data used for the LASSO
#' sample-split. The ANOVA sample-split is \code{1 - p.samp1}.
#' @param lambda.opt Criterion for optimum selection of cross-validated in the
#' LASSO data selection step. Either 'lambda.1se' (default) or 'lambda.min'.
#' See \code{\link[glmnet]{cv.glmnet}} for more details.
#' @param nfolds Number of folds for the LASSO cross-validated (default is 5),
#' see \code{\link[glmnet]{cv.glmnet}} for details.
#' @param gamma Vector of gamma-values used in significant estimation.
#' @param max.p.esti Maximum for computed p-values. All p-values above this value are set
#' to one. Small \code{max.p.esti} values reduce computing time.
#' @param mc.cores Number of cores for parallelising. Theoretical maximum is
#' \code{'B'}. For details see \code{\link[parallel]{mclapply}}.
#' @param trace If \code{TRUE} it prints current status of the program.
#' @param ... additional arguments for \code{\link[glmnet]{cv.glmnet}}.
#'
#' @examples
#' # files containing example data for SNP data and the phenotype
#' pfile <- system.file("extdata/phenodata.csv", package = "qtcat")
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' pdat <- read.csv(pfile, header = TRUE)
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' geno <- qtcatGeno(snp, clust)
#' pheno <- qtcatPheno(pdat)
#' # fitting HIT
#' fitted <- qtcatHit(pheno, geno)
#'
#' @importFrom methods is
#' @export
qtcatHit <- function(pheno, geno, B = 50, p.samp1 = 0.5,
                     lambda.opt = c("lambda.1se", "lambda.min"), nfolds = 5,
                     gamma = seq(.05, .99, by = .01),
                     max.p.esti = 1, mc.cores = 1L, trace = FALSE, ...) {
  stopifnot(is(pheno, "qtcatPheno"))
  stopifnot(is(geno, "qtcatGeno"))
  id <- intersect(pheno$names, rownames(geno$x))
  if (!length(id))
    stop("The ID intersect of 'pheno' and 'geno' is emty")
  if (length(id.uniqueGeno <- setdiff(rownames(geno$x), id)))
    cat("The following individuals are part of 'geno' but not of 'pheno':\n",
            paste(id.uniqueGeno, collapse = " "), "\n")
  if (length(id.uniquePheno <- setdiff(pheno$names, id)))
    cat("The following individuals are part of 'pheno' but not of 'geno':\n",
            paste(id.uniquePheno, collapse = " "), "\n")
  phenoInx <- which(pheno$names %in% id)
  genoInx <- match(pheno$names[phenoInx], rownames(geno$x))
  if (ncol(pheno$design) == 0L)
    x <- geno$x[genoInx, ]
  else
    x <- cbind(geno$x[genoInx, ], pheno$design[phenoInx, ])
  y <- pheno$pheno[phenoInx]
  fitHit <- hit(x, y, geno$hierarchy,
                B, p.samp1, lambda.opt, nfolds, gamma, max.p.esti, mc.cores,
                trace, standardize = FALSE)
  out <- c(fitHit,
           geno[3L:5L])
  class(out) <- c("hit", "qtcatHit")
  out
}


#' @title Summarize results of Hierarchical Inference Test
#'
#' @description Summarizing the QTCs (significant cluster of SNPs) and their
#' position at the genome.
#'
#' @param x An object of class \code{\link{qtcatHit}}.
#' @param alpha An alpha level for significance estimation.
#' @param min.absCor Minimum absolute value of correlation considered.
#'
#' @examples
#' # files containing example data for SNP data and the phenotype
#' pfile <- system.file("extdata/phenodata.csv", package = "qtcat")
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' pdat <- read.csv(pfile, header = TRUE)
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' geno <- qtcatGeno(snp, clust)
#' pheno <- qtcatPheno(pdat)
#' fitted <- qtcatHit(pheno, geno)
#' # Summarizing the QTCs
#' qtcatQtc(fitted)
#'
#' @importFrom hit hit
#' @importFrom methods is
#' @export
qtcatQtc <- function(x, alpha = 0.05, min.absCor) {
  stopifnot(is(x, "qtcatHit"))
  if (missing(min.absCor))
    y <- summary(x, alpha)
  else
    y <- summary(x, alpha, 1 - min.absCor)
  signames <- rownames(y)
  sigclust <- match(signames, x$medoids)
  sigClust <- matrix(0, length(x$clusters), 3L)
  sigClust[, 2L] <- 1
  for (i in seq_along(signames)) {
    inx <- which(x$clusters == sigclust[i])
    sigClust[inx, ] <- matrix(unlist(y[signames[i], ]),
                              length(inx), 3L, byrow = TRUE)
  }
  rownames(sigClust) <- names(x$clusters)
  sigClust[, 2L] <- 1 - sigClust[, 2L]
  colnames(sigClust) <- c(colnames(y)[1L], "absCor", colnames(y)[3L])
  sigClust <- as.data.frame(sigClust[sigClust[, 1L] != 0, ,drop = FALSE])
  out <- cbind(t(x$positions[, rownames(sigClust), drop = FALSE]),
                    sigClust)
  out
}


#' @title Fitting a Linear Model to medoids
#'
#' @description Linear model between phenotype and medoids of QTCs (significant
#' SNP clusters).
#'
#' @param x An object of class \code{\link{qtcatHit}}.
#' @param pheno An object of class \code{\link{qtcatPheno}}.
#' @param geno An object of class \code{\link{qtcatGeno}}.
#' @param alpha An alpha level for significance estimation.
#' @param min.absCor Minimum absolute value of correlation considered.
#'
#'@examples
#' # files containing example data for SNP data and the phenotype
#' pfile <- system.file("extdata/phenodata.csv", package = "qtcat")
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' pdat <- read.csv(pfile, header = TRUE)
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' geno <- qtcatGeno(snp, clust)
#' pheno <- qtcatPheno(pdat)
#' fitted <- qtcatHit(pheno, geno)
#' # fitting a LM to the phenotype and QTC medoids
#' lmfitted <- medoidsLm(fitted, pheno, geno)
#'
#' @importFrom hit hit
#' @importFrom methods is
#' @export
medoidsLm <- function(x, pheno, geno, alpha = 0.05, min.absCor) {
  stopifnot(is(x, "qtcatHit"))
  stopifnot(is(pheno, "qtcatPheno"))
  stopifnot(is(geno, "qtcatGeno"))
  id <- intersect(pheno$names, rownames(geno$x))
  if (!length(id))
    stop("The ID intersect of 'pheno' and 'geno' is emty")
  if (length(id.uniqueGeno <- setdiff(rownames(geno$x), id)))
    cat("The following individuals are part of 'geno' but not of 'pheno':\n",
        paste(id.uniqueGeno, collapse = " "), "\n")
  if (length(id.uniquePheno <- setdiff(pheno$names, id)))
    cat("The following individuals are part of 'pheno' but not of 'geno':\n",
        paste(id.uniquePheno, collapse = " "), "\n")
  sigClust <- summary(x, alpha, min.absCor)
  clusters <- split(rownames(sigClust), sigClust$clusters)
  medoids <- lapply(clusters, function(names, geno) {
    if (length(names) > 1L)
      return(names[which.max(abs(cor(geno$x[, names])))])
    else
      return(names)
  }, geno = geno)
  xg <- geno$x[, colnames(geno$x) %in% unlist(medoids)]
  phenoInx <- which(pheno$names %in% id)
  genoInx <- match(pheno$names[phenoInx], rownames(xg))
  rownames(xg) <- NULL
  if (is.null(pheno$design))
    dat <- data.frame(y = pheno$pheno[phenoInx], xg[genoInx, ])
  else
    dat <- data.frame(y = pheno$pheno[phenoInx],
                      pheno$design[phenoInx, ],
                      xg[genoInx, ])
    lm(y ~ ., data = dat)
}
