## ============================================================
## WGCNA: weighted gene co-expression network analysis
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Builds a co-expression network from the RNA-seq counts, correlates
## modules with treatment (CC/CV/VV) and pH, and extracts genes from
## the module of interest for comparison against the DESeq2 DEGs and
## the RDA candidate SNPs under selection.
##
## FIXES vs the original draft:
##   - datTraits now reads all trait columns actually used later
##     (Subtype, InitialpH, FinalpH), not just TreatmentA/B/C - the
##     original only read 3 columns then referenced others that were
##     never loaded, which would have errored.
##   - The module-trait heatmap ("Step 5") now correlates against a
##     proper numeric design matrix (one column per treatment group,
##     built with model.matrix) instead of passing raw text columns
##     straight into cor(), which requires numeric input.
##   - Removed a stray, unfinished `levels` statement and a
##     `groupLabels = colnames(sample)` reference to an undefined
##     "sample" object (replaced with a real, defined label).
##
## Input:
##   - counts_curated_CV4VV5rem.csv  (raw counts, genes x samples)
##   - metadata1.csv                 (Sample, TreatmentA/B/C, Subtype, InitialpH, FinalpH)
##   - merged_genes_VvsA.csv / merged_genes_AVvsA.csv / merged_genes_AVvsV.csv
##       (DESeq2 log2FC per contrast, from 06_DESeq2.R)
##   - Pliv_genes_master_filt.csv    (eggnog annotation, for the GO filtering step)
##
## Output:
##   - step5-Module-trait-relationships.png
##   - magenta_genes.csv, datExpr_magenta.csv, filtered_genes.csv
##   - heatmap_output.png
##   - module_magenta_results.csv
## ============================================================

## ------------------------------------------------------------
## 0. Libraries
## ------------------------------------------------------------

# install.packages("Hmisc", type = "binary")                     # uncomment if needed
# BiocManager::install(c("impute", "preprocessCore"))             # uncomment if needed
# install.packages("WGCNA", dependencies = TRUE)                  # uncomment if needed

library(Hmisc)
library(impute)
library(preprocessCore)
library(WGCNA)
library(dplyr)
library(pheatmap)
options(stringsAsFactors = FALSE)


## ------------------------------------------------------------
## 1. Load expression data and sample traits
## ------------------------------------------------------------

data <- read.csv("counts_curated_CV4VV5rem.csv", row.names = 1, sep = ";")
datExpr <- as.data.frame(t(data))  # WGCNA expects samples in rows, genes in columns
datExpr <- as.data.frame(lapply(datExpr, as.numeric))
rownames(datExpr) <- colnames(data)

metadata <- read.csv("metadata1.csv", sep = ";", header = TRUE)
head(metadata)

# Keep every trait column actually used downstream (treatment groups,
# sample subtype, and pH at start/end of the experimental treatment)
datTraits <- data.frame(
  gsm = metadata$Sample,
  TreatmentA = metadata$TreatmentA,
  TreatmentB = metadata$TreatmentB,
  TreatmentC = metadata$TreatmentC,
  Subtype = metadata$Subtype,
  InitialpH = metadata$InitialpH,
  FinalpH = metadata$FinalpH
)
save(datTraits, file = "datTraits.RData")

# Match datTraits rows to datExpr sample order
traitRows <- match(rownames(datExpr), datTraits$gsm)
datTraits <- datTraits[traitRows, ]
rownames(datTraits) <- datTraits$gsm
head(datTraits)


## ------------------------------------------------------------
## 2. Preprocessing and outlier check
## ------------------------------------------------------------

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5)


## ------------------------------------------------------------
## 3. Soft-thresholding power selection
## ------------------------------------------------------------

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     type = "n", xlab = "Soft threshold (power)", ylab = "Scale-free topology model fit",
     main = "Soft thresholding power")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.9, col = "blue")

softPower <- 6  # pick from the plot above (first power where fit crosses ~0.9)


## ------------------------------------------------------------
## 4. Network construction and module detection
## ------------------------------------------------------------

net <- blockwiseModules(
  datExpr,
  power = softPower,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)
table(net$colors)

moduleColors <- labels2colors(net$colors)
table(moduleColors)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                     "Module colors", dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)


## ------------------------------------------------------------
## 5. Sample dendrogram with a trait heatmap
## ------------------------------------------------------------

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

datExpr_tree <- hclust(dist(datExpr), method = "average")
par(mar = c(0, 5, 2, 0))
plot(datExpr_tree, main = "Sample clustering", sub = "", xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1)

combined_treatments <- factor(paste(datTraits$TreatmentA, datTraits$TreatmentB, datTraits$TreatmentC,
                                     datTraits$InitialpH, datTraits$FinalpH, sep = ";"))
sample_colors <- numbers2colors(as.numeric(combined_treatments),
                                 colors = c("white", "blue", "red", "green", "yellow", "purple", "orange"),
                                 signed = FALSE)

par(mar = c(1, 4, 3, 1), cex = 0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                     groupLabels = "Treatment (combined)",
                     cex.dendroLabels = 0.8, marAll = c(1, 4, 3, 1), cex.rowText = 0.01,
                     main = "Sample dendrogram and trait heatmap")


## ------------------------------------------------------------
## 6. Module-trait relationships
## ------------------------------------------------------------

# Build a numeric design matrix: one 0/1 column per treatment group
# (cor() needs numeric input, not the raw text trait columns)
design <- model.matrix(~ 0 + TreatmentA, data = datTraits)
colnames(design) <- levels(factor(datTraits$TreatmentA))

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, design, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

png("step5-Module-trait-relationships.png", width = 800, height = 1200, res = 120)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()


## ------------------------------------------------------------
## 7. Extract genes from the module of interest
##
## NOTE: "magenta" is the module Robert identified as biologically
## relevant from the heatmap above - if a re-run picks a different
## module as the interesting one, change `module` below accordingly.
## ------------------------------------------------------------

module <- "magenta"
names(moduleColors) <- colnames(datExpr)

magenta_genes <- colnames(datExpr)[moduleColors == module]
head(magenta_genes)

datExpr_magenta <- datExpr[, magenta_genes]
write.csv(magenta_genes, "magenta_genes.csv", row.names = FALSE, quote = FALSE)
write.csv(datExpr_magenta, "datExpr_magenta.csv", row.names = TRUE, quote = FALSE)


## ------------------------------------------------------------
## 8. Filter the module's genes against the functional annotation
## ------------------------------------------------------------

genes_master <- read.csv("Pliv_genes_master_filt.csv", sep = ";", stringsAsFactors = FALSE)
filtered_genes <- genes_master %>% filter(gene_id %in% magenta_genes)
write.csv(filtered_genes, "filtered_genes.csv", row.names = FALSE)


## ------------------------------------------------------------
## 9. Heatmap of the module's expression
## ------------------------------------------------------------

dat <- t(datExpr[, moduleColors == module])
n <- t(scale(t(log(dat + 1))))
n[n > 2] <- 2
n[n < -2] <- -2

ac <- data.frame(Treatment = datTraits$Subtype)
rownames(ac) <- colnames(n)

png("heatmap_output.png", width = 1000, height = 800, res = 120)
pheatmap(n, show_colnames = TRUE, show_rownames = FALSE, annotation_col = ac, cluster_rows = FALSE)
dev.off()


## ------------------------------------------------------------
## 10. Gene significance per contrast, cross-referenced with DESeq2
## ------------------------------------------------------------

res_VvsA  <- read.csv("merged_genes_VvsA.csv",  stringsAsFactors = FALSE)[, c(1, 3)]
res_AVvsA <- read.csv("merged_genes_AVvsA.csv", stringsAsFactors = FALSE)[, c(1, 3)]
res_AVvsV <- read.csv("merged_genes_AVvsV.csv", stringsAsFactors = FALSE)[, c(1, 3)]
colnames(res_VvsA)  <- c("Gene", "log2FC_VvsA")
colnames(res_AVvsA) <- c("Gene", "log2FC_AVvsA")
colnames(res_AVvsV) <- c("Gene", "log2FC_AVvsV")

datExpr_t <- t(datExpr_magenta)  # genes as rows for the correlation step below
metadata <- metadata[match(colnames(datExpr_t), metadata$Sample), ]

# Binary trait vectors, one per comparison (1 = first group, 0 = second group)
trait_VvsC  <- ifelse(metadata$Subtype == "VV", 1, ifelse(metadata$Subtype == "CC", 0, NA))
trait_CVvsC <- ifelse(metadata$Subtype == "CV", 1, ifelse(metadata$Subtype == "CC", 0, NA))
trait_CVvsV <- ifelse(metadata$Subtype == "CV", 1, ifelse(metadata$Subtype == "VV", 0, NA))

genes_module <- magenta_genes[magenta_genes %in% rownames(datExpr_t)]
length(genes_module)

GS_VvsC  <- apply(datExpr_t[genes_module, ], 1, function(g) cor(g, trait_VvsC,  use = "pairwise.complete.obs"))
GS_CVvsC <- apply(datExpr_t[genes_module, ], 1, function(g) cor(g, trait_CVvsC, use = "pairwise.complete.obs"))
GS_CVvsV <- apply(datExpr_t[genes_module, ], 1, function(g) cor(g, trait_CVvsV, use = "pairwise.complete.obs"))

module_results <- data.frame(Gene = genes_module, GS_VvsC, GS_CVvsC, GS_CVvsV)
module_results <- merge(module_results, res_VvsA, by = "Gene", all.x = TRUE)
module_results <- merge(module_results, res_AVvsA, by = "Gene", all.x = TRUE)
module_results <- merge(module_results, res_AVvsV, by = "Gene", all.x = TRUE)

module_results$Regulation_VvsC  <- ifelse(module_results$GS_VvsC > 0, "Up in VV", "Up in CC")
module_results$Regulation_CVvsC <- ifelse(module_results$GS_CVvsC > 0, "Up in CV", "Up in CC")
module_results$Regulation_CVvsV <- ifelse(module_results$GS_CVvsV > 0, "Up in CV", "Up in VV")

module_results <- module_results[order(-abs(module_results$GS_VvsC)), ]
write.csv(module_results, "module_magenta_results.csv", row.names = FALSE)
head(module_results)
