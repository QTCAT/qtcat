# simulation of phenotype
#
# require(qtcat)
# set.seed(7680)
#
# snp <- read.snpData("./inst/extdata/snpdata.csv", sep = ",")
# snpmat <- as.matrix(snp)
# genotype <- snpmat[,c(37, 260, 367)] %*% c(2, 2, 2)
# pheno <- rbind(genotype, genotype, genotype)
# pheno <- data.frame(name = rownames(pheno),
#                     phenotype = as.vector(pheno[, 1]),
#                     replication = rep(x = c("I", "II", "III"),
#                                       each = nrow(genotype)))
# pheno$phenotype <- pheno$phenotype +
#   rep(x = c(0, 1.5, 2), nrow(genotype)) +
#   rnorm(nrow(pheno), 0, 1)
# write.csv(pheno, "./inst/extdata/phenodata.csv",
#           quote = FALSE, row.names = FALSE)
#
# example("qtcatSigClust")
# phenofile <- system.file("extdata/phenodata.csv", package = "qtcat")
