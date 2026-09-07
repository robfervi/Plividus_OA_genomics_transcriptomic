## ============================================================
## PCA on ddRADseq SNP data
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Run separately on the neutral SNP dataset and on the
## candidate SNPs under pH selection (RDA outliers, z > 2.5),
## to compare clustering patterns between the two datasets.
##
## Population assignment is read from popmap.txt (sample \t population),
## the same file used in the Ho/He/Fis script, instead of a hardcoded
## vector — this avoids silent mismatches if individual order in the
## VCF ever changes.
##
## Input:
##   - populationsNeutro.recode.vcf      (neutral SNPs)
##   - populationsUS_RDA2_5.recode.vcf   (candidate SNPs under selection)
##   - popmap.txt                        (sample, population)
##
## Output:
##   - PCA plots (PC1 vs PC2) for each dataset
## ============================================================

library(vcfR)
library(adegenet)
library(ggplot2)


## ------------------------------------------------------------
## glPcaFast: faster PCA for genlight objects
## (standard replacement for adegenet::glPca on large SNP datasets;
## same algorithm, optimized for speed — see adegenet forum/tutorials)
## ------------------------------------------------------------

glPcaFast <- function(x,
                       center = TRUE,
                       scale = FALSE,
                       nf = NULL,
                       loadings = TRUE,
                       alleleAsUnit = FALSE,
                       returnDotProd = FALSE) {

  if (!inherits(x, "genlight")) stop("x is not a genlight object")

  if (center) {
    vecMeans <- glMean(x, alleleAsUnit = alleleAsUnit)
    if (any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }

  if (scale) {
    vecVar <- glVar(x, alleleAsUnit = alleleAsUnit)
    if (any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }

  # convert to full data, keeping NA handling close to the original
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)

  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = TRUE)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[, 2]]
  } else {
    mx[NAidx] <- 0
  }

  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else FALSE,
              scale = if (scale) vecVar else FALSE)

  # all dot products at once using underlying BLAS
  allProd <- tcrossprod(mx) / nInd(x)  # assume uniform weights

  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric = TRUE, only.values = FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop = FALSE]

  ## number of axes to retain
  if (is.null(nf)) {
    barplot(eigRes$values, main = "Eigenvalues", col = heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }

  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig > 1e-10))

  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x))  # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop = FALSE], 2,
                       sqrt(eigRes$values[1:nf]), FUN = "*")

  ## loadings
  if (loadings) {
    if (scale) vecSd <- sqrt(vecVar)
    res$loadings <- matrix(0, nrow = nLoc(x), ncol = nf)
    myPloidy <- ploidy(x)
    for (k in 1:nInd(x)) {
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if (center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if (scale) temp <- temp / vecSd

      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop = FALSE]
    }

    res$loadings <- res$loadings / nInd(x)
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN = "/")
  }

  ## format output
  colnames(res$scores) <- paste("PC", 1:nf, sep = "")
  if (!is.null(indNames(x))) {
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }

  if (!is.null(res$loadings)) {
    colnames(res$loadings) <- paste("Axis", 1:nf, sep = "")
    if (!is.null(locNames(x)) & !is.null(alleles(x))) {
      rownames(res$loadings) <- paste(locNames(x), alleles(x), sep = ".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }

  if (returnDotProd) {
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }

  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}


## ------------------------------------------------------------
## Population map (sample -> population)
## ------------------------------------------------------------

popmap <- read.table("popmap.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("sample", "population")

# Site colors: V1/V2 = Vent, T1/T2 = Transition, C1/C2 = Ambient
# (see Table 1 in Materials and Methods for full site names/coordinates)
pop_colors <- c(
  "V1" = "darkorchid2", "V2" = "darkorchid2",
  "T1" = "goldenrod1",  "T2" = "goldenrod1",
  "C1" = "deepskyblue", "C2" = "deepskyblue"
)


## ------------------------------------------------------------
## Run PCA for one dataset (VCF) and return scores + plot
## ------------------------------------------------------------

run_pca <- function(vcf_file, popmap, label) {

  snps <- read.vcfR(vcf_file, verbose = FALSE)
  snps_gl <- vcfR2genlight(snps)

  # Assign populations by matching sample name (robust to sample order)
  pop(snps_gl) <- as.factor(popmap$population[match(indNames(snps_gl), popmap$sample)])

  pca <- glPcaFast(snps_gl)
  df <- as.data.frame(pca$scores)
  df$population <- pop(snps_gl)
  eig <- pca$eig

  p <- ggplot(df, aes(PC1, PC2, color = population)) +
    geom_point(size = 2.6, alpha = 1) +
    theme_classic() +
    coord_equal() +
    scale_color_manual(values = pop_colors) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(
      legend.background = element_rect(color = "black"),
      legend.title = element_blank()
    ) +
    xlab(paste0("PC1 (", sprintf("%0.1f%% explained var.", 100 * eig[1] / sum(eig)), ")")) +
    ylab(paste0("PC2 (", sprintf("%0.1f%% explained var.", 100 * eig[2] / sum(eig)), ")")) +
    ggtitle(label)

  list(pca = pca, df = df, plot = p)
}


## ------------------------------------------------------------
## Run for both datasets
## ------------------------------------------------------------

pca_neutral   <- run_pca("populationsNeutro.recode.vcf",    popmap, "Neutral SNPs")
pca_selection <- run_pca("populationsUS_RDA2_5.recode.vcf", popmap, "Candidate SNPs under selection")

pca_neutral$plot
pca_selection$plot
