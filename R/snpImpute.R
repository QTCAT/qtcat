#' @title Impute missing information at each SNP
#'
#' @description Uses neighboring SNPs in the clustering hierarchy to impute alleles to
#' positions with missing values.
#'
#' @param snp An object of class \linkS4class{snpData}.
#' @param snpClust An object of class \code{\link{qtcatClust}}.
#' @param min.absCor A minimum value of correlation. If missing values still exist if this
#' point in the hierarchy is reached, imputing is done via allele frequencies.
#'
#' @importFrom hit as.hierarchy
#' @importFrom stats reorder
#' @export
imputSnpData <- function(snp, snpClust, min.absCor = .25) {
  stopifnot(is(snp, "snpData"))
  stopifnot(is(snpClust, "qtcatClust"))
  hier <- reorder(as.hierarchy(snpClust$dendrogram), colnames(snp))
  snp <- imputeMenoids(snp, snpClust$clusters, hier, min.absCor)
  # impute non medoid SNPs (if exist)
  nonMedo <- which(!names(snpClust$clusters) %in% snpClust$medoids)
  if (length(nonMedo))
    for (i in nonMedo) {
      m <- snpClust$clusters[i]
      pos <- which(snp@snpData[, i] == is.raw(0))
      snp@snpData[pos, i] <- snp@snpData[pos, m]
    }
  snp
}


#' @title Impute missing information at medoid SNPs
#'
#' @description Uses neighboring SNPs in the clustering hierarchy to impute alleles to
#' positions with missing values at medoid SNPs.
#'
#' @param snp An object of class \linkS4class{snpData}.
#' @param clust A named vector of clusters.
#' @param hier A object of class hierarchy.
#' @param min.absCor A minimum value of correlation. If missing values still exist if this
#' point in the hierarchy is reached, imputing is done via allele frequencies.
#'
#' @keywords internal
imputeMenoids <- function(SNP, clust, hier, min.absCor = .25) {
  # function which is altering the 'SNP' input by imputing info from 'grp' to 'snpOfInt'
  imputSnp <- function(snpOfInt, snpsToComp, unsolved) {
    n <- length(snpsToComp)
    j <- 1
    while (unsolved && j <= n) {
      zeroPos <- which(SNP@snpData[, snpOfInt] == as.raw(0))
      snpInfoPos <- SNP@snpData[zeroPos, snpsToComp[j]]
      SNP@snpData[zeroPos, snpOfInt] <<- snpInfoPos
      unsolved <- any(snpInfoPos == as.raw(0))
      j <- j + 1
    }
    unsolved
  }
  leafs_hier <- which(sapply(hier, function(x) is.null(attr(x, which = "subset"))))
  leafs_hiersnp <- unlist(hier[leafs_hier])
  naSnps <- which(naFreq(SNP, 2) > 0 & colnames(SNP) %in% labels(hier))
  # run thru all SNPs with NAs
  for (i in seq_along(naSnps)) { # i <- 1
    unsolved <- TRUE
    snpOfInt_snp <- snpsNotComp <- naSnps[i]
    # check in clusters of identicals
    snpGrp <- which(clust == clust[snpOfInt_snp])
    snpsToComp <- snpGrp[!(snpGrp %in% snpsNotComp)]
    if (length(snpsToComp))
      unsolved <- imputSnp(snpOfInt_snp, snpsToComp, unsolved)
    if (unsolved) {
      # run thru the heirarchy until NAs of the SNP are filled with  information or the
      # height threshold is reached
      snpOfInt_hier <- leafs_hier[leafs_hiersnp == snpOfInt_snp]
      super <- attr(hier[[snpOfInt_hier]], "superset")
      snpGrp <- hier[[super]]
      h <- attr(snpGrp, "height")
      while (unsolved && h <= min.absCor) {
        snpsNotComp <- c(snpsNotComp, snpsToComp)
        snpsToComp <- snpGrp[!(snpGrp %in% snpsNotComp)]
        if (length(snpsToComp))
          unsolved <- imputSnp(snpOfInt_snp, snpsToComp, unsolved)
        super <- attr(hier[[super]], "superset")
        snpGrp <- hier[[super]]
        h <- attr(snpGrp, "height")
      }
    }
    # if height threshold is reached use alle frequency for imputing
    if (unsolved) {
      zeroPos <- which(SNP@snpData[, snpOfInt_snp] == as.raw(0))
      alleleNo <- table(as.integer(SNP@snpData[, snpOfInt_snp]), exclude = 0L)
      alleles <- as.raw(names(alleleNo))
      prob <- alleleNo / sum(alleleNo)
      SNP@snpData[zeroPos, snpOfInt_snp] <- sample(alleles, length(zeroPos), TRUE, prob)
    }
  }
  SNP
}
