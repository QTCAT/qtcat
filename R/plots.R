#' @title Plot resulting QTCs of the Hierarchical Inference Test
#'
#' @description Plot the QTCs (significant cluster of SNPs) at their
#' position at the genome.
#'
#' @param x An object of class \code{\link{qtcatHit}}.
#' @param alpha An alpha level for significance estimation.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param col.axis Colors for axis line, tick marks, and title respectively.
#' @param ... Other graphical parameters may also be passed as arguments to this function.
#'
#' @examples
#' # If you want to run the examples, use:
#' # example(plotQtc, run.dontrun = TRUE)
#' \dontrun{
#' # files containing example data for SNP data and the phenotype
#' pfile <- system.file("extdata/phenodata.csv", package = "qtcat")
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' pdat <- read.csv(pfile, header = TRUE)
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' geno <- qtcatGeno(snp, clust)
#' pheno <- qtcatPheno(names = pdat[, 1],
#'                     pheno = pdat[, 2],
#'                     covariates = model.matrix(~ pdat[, 3]))
#' fitted <- qtcatHit(pheno, geno)
#'
#' # Plot the QTCs (loci37, loci260, and loci367 are causal)
#' plotQtc(fitted)
#' }
#'
#' @importFrom graphics plot axis mtext
#' @export
plotQtc <- function(x, alpha = 0.05, xlab = "Chromosomes",
                    ylab = expression(-log[10](italic(p))), col.axis = NULL, ...) {
  stopifnot(is(x, "qtcatHit"))
  # make positions linear with gaps between chr's
  pos <- t(x$positions)
  chrminmax <- vapply(split(pos[, 2L], pos[, 1L]), function(x) c(min(x), max(x)), c(1, 2))
  chrgap <- sum(chrminmax[2L, ]) * .01
  chrsize <- cumsum(c(0, chrminmax[2L, -ncol(chrminmax)])) -
    cumsum(chrminmax[1L, ]) +
    cumsum(c(0, rep(chrgap, ncol(chrminmax) - 1L)))
  chr <- unique(pos[, 1L])
  for (i in seq_along(chr)) {
    inx <- which(pos[, 1L] == chr[i])
    pos[inx, 2L] <- pos[inx, 2L] + chrsize[i]
  }
  chrstartend <- vapply(split(pos[, 2L], pos[, 1L]), function(x) c(min(x), max(x)), c(1, 2))
  xlim <- c(chrstartend[1L, 1L] - chrgap, chrstartend[2L, ncol(chrstartend)] + chrgap)
  # get QTCs from result
  qtc <- qtcatQtc(x, alpha = alpha, min.absCor = .01)
  for (i in seq_along(chr)) {
    inx2 <- which(qtc[, 1L] == chr[i])
    qtc[inx2, 2L] <- qtc[inx2, 2L] + chrsize[i]
  }
  qtc$pValues[qtc$pValues] <- qtc$pValues[qtc$pValues] + 1e-308
  qtc$log.pValues <- -log10(qtc$pValues)
  # plot
  plot(qtc$pos, qtc$log.pValues,  xlim = xlim, ylim = c(0, max(9, qtc$log.pValues) + 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  # x
  for (i in seq_along(chr))
    axis(1, labels = FALSE, at = c(chrstartend[1, i], chrstartend[2, i]), col = col.axis)
  axis(1, at = colMeans(chrstartend), labels = chr, col = NA, col.axis = col.axis)
  mtext(xlab, 1, 2.5, col = col.axis)
  # y
  axis(2, col.axis = col.axis, col = col.axis)
  mtext(expression(-log[10](italic(p))), 2, 2.5, col = col.axis)
}

#' @title Plot markers selection frequencies of the Hierarchical Inference Test
#'
#' @description Plot markers selection frequencies at their
#' position at the genome.
#'
#' @param x An object of class \code{\link{qtcatHit}}.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param col.axis Colors for axis line, tick marks, and title respectively.
#' @param ... Other graphical parameters may also be passed as arguments to this function.
#'
#' @examples
#' # If you want to run the examples, use:
#' # example(plotSelFreq, run.dontrun = TRUE)
#' \dontrun{
#' # files containing example data for SNP data and the phenotype
#' pfile <- system.file("extdata/phenodata.csv", package = "qtcat")
#' gfile <- system.file("extdata/snpdata.csv", package = "qtcat")
#' pdat <- read.csv(pfile, header = TRUE)
#' snp <- read.snpData(gfile, sep = ",")
#' clust <- qtcatClust(snp)
#' geno <- qtcatGeno(snp, clust)
#' pheno <- qtcatPheno(names = pdat[, 1],
#'                     pheno = pdat[, 2],
#'                     covariates = model.matrix(~ pdat[, 3]))
#' fitted <- qtcatHit(pheno, geno)
#'
#' # Plot the selection frequncy of markers (loci37, loci260, and loci367 are causal)
#' plotSelFreq(fitted)
#' }
#'
#' @importFrom graphics plot axis mtext
#' @export
plotSelFreq <- function(x, xlab = "Chromosomes", ylab = "Sel. freq.",
                        col.axis = NULL, ...) {
  stopifnot(is(x, "qtcatHit"))
  # make positions linear with gaps between chr's
  pos <- t(x$positions)
  chrminmax <- vapply(split(pos[, 2L], pos[, 1L]), function(x) c(min(x), max(x)), c(1, 2))
  chrgap <- sum(chrminmax[2L, ]) * .01
  chrsize <- cumsum(c(0, chrminmax[2L, -ncol(chrminmax)])) -
    cumsum(chrminmax[1L, ]) +
    cumsum(c(0, rep(chrgap, ncol(chrminmax) - 1L)))
  chr <- unique(pos[, 1L])
  for (i in seq_along(chr)) {
    inx <- which(pos[, 1L] == chr[i])
    pos[inx, 2L] <- pos[inx, 2L] + chrsize[i]
  }
  chrstartend <- vapply(split(pos[, 2L], pos[, 1L]), function(x) c(min(x), max(x)), c(1, 2))
  xlim <- c(chrstartend[1L, 1L] - chrgap, chrstartend[2L, ncol(chrstartend)] + chrgap)
  # hit lasso selection fequency
  inx <- which(x$selectFreq > 0)
  selfreq <- cbind(pos[names(x$hier), 2L][inx], x$selectFreq[inx])
  # hit lasso selection fequency plot
  plot(selfreq, xlim = xlim, ylim = c(0, 1.1), axes = FALSE, xlab = "", ylab = "", ...)
  # x
  for (i in seq_along(chr))
    axis(1, labels = FALSE, at = c(chrstartend[1, i], chrstartend[2, i]), col = col.axis)
  axis(1, at = colMeans(chrstartend), labels = chr, col = NA, col.axis = col.axis)
  mtext(xlab, 1, 2.5, col = col.axis)
  # y
  axis(2, at = c(0, .5,  1), col = col.axis, col.axis = col.axis)
  mtext(ylab, 2, 2.5, col = col.axis)
}

