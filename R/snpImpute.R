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
imputSnpData <- function(snp, snpClust, min.absCor = .25) {
  stopifnot(is(snp, "snpData"))
  stopifnot(is(snpClust, "qtcatClust"))
  snpnames <- colnames(snp)
  hier <- as.hierarchy(snpClust$dendrogram, names = snpnames)
  flipAlleles <- as.numeric(alleleFreq(snp, FALSE) <= .5)
  snp <- imputeMedo(snp, snpClust$clusters, hier, flipAlleles, min.absCor)
  # impute non medoid SNPs (if exist)
  nonMedo <- which(!(names(snpClust$clusters) %in% snpClust$medoids))
  if (length(nonMedo))
    for (i in nonMedo) {
      m <- which(snpClust$medoids[snpClust$clusters[i]] == snpnames)
      js <- which(snp@snpData[, i] == is.raw(0))
      jAllele <- snp@snpData[js, m]
      if (flipAlleles[i] != flipAlleles[m]) {
        j1 <- which(jAllele == as.raw(1))
        j3 <- which(jAllele == as.raw(3))
        jAllele[j1] <- as.raw(3)
        jAllele[j3] <- as.raw(1)
      }
      snp@snpData[js, i] <- jAllele
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
#' @param flipAlleles A vertor of telling for each SNP if allele one has allele freq. > 0.5 
#' or not.
#' @param min.absCor A minimum value of correlation. If missing values still exist if this
#' point in the hierarchy is reached, imputing is done via allele frequencies.
#'
#' @keywords internal
imputeMedo <- function(SNP, clust, hier, flipAlleles, min.absCor = .25) {
  # function which is altering the 'SNP' input by imputing info from 'grp' to 'snpOfInt'
  imputSnp <- function(snpOfInt, snpsToComp, snpOfIntFlip, flipAlleles, unsolved) {
    n <- length(snpsToComp)
    i <- 1
    while (unsolved && i <= n) {
      js <- which(SNP@snpData[, snpOfInt] == as.raw(0))
      jAllele <- SNP@snpData[js, snpsToComp[i]]
      if (snpOfIntFlip != flipAlleles[snpsToComp[i]]) {
        j1 <- which(jAllele == as.raw(1))
        j3 <- which(jAllele == as.raw(3))
        jAllele[j1] <- as.raw(3)
        jAllele[j3] <- as.raw(1)
      }
      SNP@snpData[js, snpOfInt] <<- jAllele
      unsolved <- any(jAllele == as.raw(0))
      i <- i + 1
    }
    unsolved
  }
  max.h <- 1 - min.absCor
  leafs_hier <- which(sapply(hier, function(x) is.null(attr(x, which = "subset"))))
  leafs_hiersnp <- hier[leafs_hier]
  naSnps <- which(naFreq(SNP, 2) > 0 & colnames(SNP) %in% labels(hier))
  # run thru all SNPs with NAs
  for (i in seq_along(naSnps)) {
    unsolved <- TRUE
    snpOfInt_snp <- snpsNotComp <- naSnps[i]
    snpOfIntFlip <- flipAlleles[snpOfInt_snp]
    # check in clusters of identicals
    snpGrp <- which(clust == clust[snpOfInt_snp])
    snpsToComp <- snpGrp[!(snpGrp %in% snpsNotComp)]
    if (length(snpsToComp))
      unsolved <- imputSnp(snpOfInt_snp, snpsToComp, snpOfIntFlip, flipAlleles, unsolved)
    if (unsolved) {
      # run thru the heirarchy until NAs of the SNP are filled with  information or the
      # height threshold is reached
      snpOfInt_hier <- leafs_hier[sapply(leafs_hiersnp, function(x) any(x == snpOfInt_snp))]
      super <- attr(hier[[snpOfInt_hier]], "superset")
      snpGrp <- hier[[super]]
      h <- attr(snpGrp, "height")
      while (unsolved && h <= max.h) {
        snpsNotComp <- c(snpsNotComp, snpsToComp)
        snpsToComp <- snpGrp[!(snpGrp %in% snpsNotComp)]
        if (length(snpsToComp))
          unsolved <- imputSnp(snpOfInt_snp, snpsToComp, snpOfIntFlip, 
                               flipAlleles, unsolved)
        super <- attr(hier[[super]], "superset")
        snpGrp <- hier[[super]]
        h <- attr(snpGrp, "height")
      }
    }
    # if height threshold is reached use alle frequency for imputing
    if (unsolved) {
      js <- which(SNP@snpData[, snpOfInt_snp] == as.raw(0))
      alleleNo <- table(as.integer(SNP@snpData[, snpOfInt_snp]), exclude = 0L)
      alleles <- as.raw(names(alleleNo))
      prob <- alleleNo / sum(alleleNo)
      SNP@snpData[js, snpOfInt_snp] <- sample(alleles, length(js), TRUE, prob)
    }
  }
  SNP
}
