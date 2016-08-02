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
#' @export
imputeSnpData <- function(snp, snpClust, min.absCor = .25) {
  stopifnot(is(snp, "snpData"))
  stopifnot(is(snpClust, "qtcatClust"))
  snpnames <- colnames(snp)
  hier <- as.hierarchy(snpClust$dendrogram, names = snpnames)
  snp <- imputeMedo(snp, snpClust$clusters, hier, min.absCor)
  # impute non medoid SNPs (if exist)
  nonMedo <- which(!(names(snpClust$clusters) %in% snpClust$medoids))
  if (length(nonMedo)) {
    flipAlleles <- as.numeric(alleleFreq(snp, FALSE) <= .5)
    snpList <- list()
    for (i in 1:ncol(snp)) {
      snpList[[i]] <- snp@snpData[, i]
      if (i %in% nonMedo) {
        m <- which(snpClust$medoids[snpClust$clusters[i]] == snpnames)
        js <- which(snpList[[i]] == is.raw(0))
        jAllele <- snp@snpData[js, m]
        if (flipAlleles[i] != flipAlleles[m]) {
          j1 <- which(jAllele == as.raw(1))
          j3 <- which(jAllele == as.raw(3))
          jAllele[j1] <- as.raw(3)
          jAllele[j3] <- as.raw(1)
        }
        snpList[[i]][js] <- jAllele
      }
    }
    snp@snpData <- do.call(cbind, snpList)
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
imputeMedo <- function(snp, clust, hier, min.absCor = .25) {
  snpList <- list()
  max.h <- 1 - min.absCor
  leafs_hier <- which(sapply(hier, function(x) is.null(attr(x, which = "subset"))))
  leafs_hiersnp <- hier[leafs_hier]
  naSnps <- which(naFreq(snp, 2) > 0 & colnames(snp) %in% labels(hier))
  flipAlleles <- as.numeric(alleleFreq(snp, FALSE) <= .5)
  # run thru all SNPs with NAs
  for (inxSnpOfInt in 1:ncol(snp)) { # inxSnpOfInt <- 1
    snpList[[inxSnpOfInt]] <- snp@snpData[, inxSnpOfInt]
    if (inxSnpOfInt %in% naSnps) {
      unsolved <- TRUE
      inxSnpsNotComp <- inxSnpOfInt
      snpOfIntFlip <- flipAlleles[inxSnpOfInt]
      # check in clusters of identicals
      inxSnpGrp <- which(clust == clust[inxSnpOfInt])
      inxSnpsNotComp <- inxSnpGrp[!(inxSnpGrp %in% inxSnpsNotComp)]
      if (length(inxSnpsNotComp)) {
        temp <- imputeSnp(snp, snpList[[inxSnpOfInt]], inxSnpsNotComp,
                          snpOfIntFlip, flipAlleles)
        snpList[[inxSnpOfInt]] <- temp[[1L]]
        unsolved <- temp[[2L]]
      }
      if (unsolved) {
        # run thru the heirarchy until NAs of the SNP are filled with  information or the
        # height threshold is reached
        hierSnpOfInt <- leafs_hier[sapply(leafs_hiersnp, function(x) any(x == inxSnpOfInt))]
        super <- attr(hier[[hierSnpOfInt]], "superset")
        inxSnpGrp <- hier[[super]]
        h <- attr(inxSnpGrp, "height")
        while (unsolved && h <= max.h) {
          inxSnpsNotComp <- c(inxSnpsNotComp, inxSnpsNotComp)
          inxSnpsNotComp <- inxSnpGrp[!(inxSnpGrp %in% inxSnpsNotComp)]
          if (length(inxSnpsNotComp)) {
            temp <- imputeSnp(snp, snpList[[inxSnpOfInt]], inxSnpsNotComp,
                              snpOfIntFlip, flipAlleles)
            snpList[[inxSnpOfInt]] <- temp[[1L]]
            unsolved <- temp[[2L]]
          }
          super <- attr(hier[[super]], "superset")
          inxSnpGrp <- hier[[super]]
          h <- attr(inxSnpGrp, "height")
        }
      }
      # if height threshold is reached use alle frequency for imputing
      if (unsolved) {
        js <- which(snpList[[inxSnpOfInt]] == as.raw(0))
        alleleNo <- table(as.integer(snpList[[inxSnpOfInt]]), exclude = 0L)
        alleles <- as.raw(names(alleleNo))
        prob <- alleleNo / sum(alleleNo)
        snpList[[inxSnpOfInt]][js] <- sample(alleles, length(js), TRUE, prob)
      }
    }
  }
  snp@snpData <- do.call(cbind, snpList)
  snp
}


#' @title Impute missing information at a medoid SNPs from a group of neighbors
#'
#' @description Uses neighboring SNPs in the clustering hierarchy to impute as many as
#' possible alleles to positions with missing values at medoid SNPs.
#'
#' @param snp An object of class \linkS4class{snpData}.
#' @param snpOfInt A vertor of the snp of interest.
#' @param inxSnpsToComp Index of neighbors.
#' @param snpOfIntFlip Flip status of the snp of interest.
#' @param flipAlleles A vertor of telling for each SNP if allele one has allele freq. > 0.5
#' or not.
#' @param min.absCor A minimum value of correlation. If missing values still exist if this
#' point in the hierarchy is reached, imputing is done via allele frequencies.
#'
#' @keywords internal
imputeSnp <- function(snp, snpOfInt, inxSnpsToComp, snpOfIntFlip, flipAlleles) {
  unsolved <- TRUE
  n <- length(inxSnpsToComp)
  i <- 1
  while (unsolved && i <= n) {
    js <- which(snpOfInt == as.raw(0))
    jAllele <- snp@snpData[js, inxSnpsToComp[i]]
    if (snpOfIntFlip != flipAlleles[inxSnpsToComp[i]]) {
      j1 <- which(jAllele == as.raw(1))
      j3 <- which(jAllele == as.raw(3))
      jAllele[j1] <- as.raw(3)
      jAllele[j3] <- as.raw(1)
    }
    snpOfInt[js] <- jAllele
    unsolved <- any(jAllele == as.raw(0))
    i <- i + 1
  }
  list(snpOfInt, unsolved)
}
