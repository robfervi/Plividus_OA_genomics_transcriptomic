## ============================================================
## Differential gene expression analysis (DESeq2)
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
## Treatment groups: CC = A1A (ambient, chronic), VV = V1V (vent, chronic),
##                    CV = A1V (ambient origin, acute exposure to vent pH)
##
## Three pairwise comparisons (see Materials and Methods):
##   1. VV vs CC  - chronic exposure to natural low pH
##   2. CV vs CC  - acute experimental response to low pH
##   3. CV vs VV  - genotype-of-origin response to experimental low pH
##
## NOTE ON A FIX vs the original script: in the original, "up"/"down"
## regulated gene lists were labelled inconsistently across the three
## comparisons (sign of log2FoldChange did not consistently mean
## "higher in the first group" vs "higher in the second group").
## Here, for contrast(treatment, group1, group2): "up" always means
## higher in group1, "down" always means higher in group2 - check this
## matches the direction you want before citing counts in the paper.
##
## This script only performs differential expression + exports gene
## lists; GO enrichment is done separately in 07_topGO_RNA.R.
## ============================================================

## ------------------------------------------------------------
## 0. Libraries
## ------------------------------------------------------------

# install.packages(c("dplyr", "ggplot2", "pheatmap", "RColorBrewer", "ComplexHeatmap", "circlize", "glmpca"))  # uncomment if needed
# BiocManager::install(c("DESeq2", "apeglm", "ashr", "edgeR"))                                                  # uncomment if needed

library(DESeq2)
library(edgeR)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(readxl)

## ------------------------------------------------------------
## 1. Import counts and build the DESeq2 object
## ------------------------------------------------------------

countData <- read.csv("counts_curated_CV4VV5rem.csv", sep = ";")
rownames(countData) <- countData$Geneid
countData <- countData[, -1]

# Sample metadata: 8 CC, 7 CV, 7 VV (adjust if sample numbers differ)
treatment <- c(rep("CC", 8), rep("CV", 7), rep("VV", 7))
colData <- data.frame(row.names = colnames(countData), treatment = treatment)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~treatment)
dds$treatment <- factor(dds$treatment, levels = c("CC", "CV", "VV"))

# Remove genes with zero counts across all samples
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep, ]


## ------------------------------------------------------------
## 2. Exploratory analysis
## ------------------------------------------------------------

vst_data <- vst(dds, blind = FALSE)

# Sample-distance heatmap
sampleDists <- dist(t(assay(vst_data)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vst_data$treatment
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))

# Standard PCA
plotPCA(vst_data, intgroup = "treatment")

# Generalized PCA (glmpca) on raw counts - useful when variance structure is non-Gaussian
library(glmpca)
gpca <- glmpca(counts(dds), L = 2)
gpca.dat <- gpca$factors
gpca.dat$treatment <- dds$treatment
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = treatment)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


## ------------------------------------------------------------
## 3. Fit the model once
## ------------------------------------------------------------

dds <- DESeq(dds)
resultsNames(dds)

# Consistent treatment color palette used across all comparisons
treatment_colors_all <- c(CC = "#4C72B0", CV = "#F4B86E", VV = "#EC6337")


## ------------------------------------------------------------
## 4. Reusable function: one pairwise contrast, full workflow
## ------------------------------------------------------------

run_contrast <- function(dds, group1, group2, label, outdir = ".") {

  contrast_name <- paste0(group1, "vs", group2)
  message("== Running contrast: ", label, " (", contrast_name, ") ==")

  res <- results(dds, contrast = c("treatment", group1, group2))
  cat("Significant (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")

  write.csv(as.data.frame(res), file.path(outdir, paste0("deseq2_results_", contrast_name, ".csv")))

  ## Volcano plot
  res_df <- as.data.frame(res)
  res_df$significant <- with(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
  print(
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
      geom_point() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
      scale_color_manual(values = c("grey", "lightgreen")) +
      labs(title = paste("Volcano plot -", label), x = "Log2 fold change", y = "-Log10 adjusted p-value") +
      theme_minimal()
  )

  ## MA plot
  plotMA(res, main = paste("MA plot -", label))

  ## Subset dds to the two groups being compared
  sub_dds <- dds[, dds$treatment %in% c(group1, group2)]
  sub_dds$treatment <- factor(sub_dds$treatment, levels = c(group2, group1))
  n_dds <- normTransform(sub_dds, f = log2, pc = 1)

  resSig <- subset(res, padj < 0.05)

  ## Heatmap - top 100 most variable genes
  topVar <- order(rowVars(assay(sub_dds), useNames = TRUE), decreasing = TRUE)[1:min(100, nrow(sub_dds))]
  top_n_dds <- n_dds[topVar, ]
  mat_var <- assay(top_n_dds) - rowMeans(assay(top_n_dds))
  col.anno <- HeatmapAnnotation(treatment = colData(top_n_dds)$treatment,
                                 col = list(treatment = treatment_colors_all[c(group1, group2)]))
  expr_palette <- colorRamp2(c(min(mat_var), 0, max(mat_var)), c("#883BC1", "#F6E5DB", "#A53C2C"))
  draw(Heatmap(mat_var, col = expr_palette, show_row_names = FALSE, show_column_names = FALSE,
               name = "Expression", top_annotation = col.anno))

  ## Heatmap - top 100 most significant genes (if any)
  if (nrow(resSig) > 0) {
    sorted_res <- resSig[order(resSig$padj), ]
    top_sig <- na.omit(rownames(sorted_res)[1:min(100, nrow(sorted_res))])
    mat_sig <- assay(n_dds)[top_sig, ]
    mat_sig <- mat_sig - rowMeans(mat_sig)
    draw(Heatmap(mat_sig, col = expr_palette, show_row_names = FALSE, show_column_names = TRUE,
                 name = "Expression", top_annotation = col.anno, cluster_columns = FALSE))

    ## Top gene, single-gene plot
    topGene <- rownames(resSig)[which.min(resSig$padj)]
    plotCounts(sub_dds, gene = topGene, intgroup = "treatment")
  }

  ## PCA restricted to DEGs
  diff_genes <- rownames(res)[which(res$padj < 0.05)]
  if (length(diff_genes) > 1) {
    dge_dds <- n_dds[diff_genes, ]
    pca <- plotPCA(dge_dds, intgroup = "treatment", returnData = TRUE)
    percentVar <- round(100 * attr(pca, "percentVar"))
    print(
      ggplot(pca, aes(x = PC1, y = PC2, color = treatment)) +
        geom_point(size = 3) +
        labs(x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC2 (", percentVar[2], "%)")) +
        ggtitle(paste("PCA of DEGs -", label)) +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = treatment_colors_all[c(group1, group2)])
    )
  }

  ## Log fold-change shrinkage (ashr supports arbitrary contrasts, unlike apeglm)
  resLFC <- lfcShrink(dds, contrast = c("treatment", group1, group2), type = "ashr")
  plotMA(resLFC, main = paste("MA plot (shrunk LFC) -", label))

  ## Export gene lists for downstream GO enrichment (07_topGO_RNA.R)
  # "up"   = higher in group1 (log2FC > 1)
  # "down" = higher in group2 (log2FC < -1)
  all_genes <- rownames(res)
  deg_genes <- rownames(res)[which(res$padj < 0.05)]
  up_genes  <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 1)]
  down_genes <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < -1)]

  writeLines(all_genes, file.path(outdir, "all_genes.txt"))  # same across contrasts, written once is enough
  writeLines(deg_genes,  file.path(outdir, paste0("DEG_", contrast_name, ".txt")))
  writeLines(up_genes,   file.path(outdir, paste0("DEG_", contrast_name, "_up.txt")))
  writeLines(down_genes, file.path(outdir, paste0("DEG_", contrast_name, "_down.txt")))

  invisible(list(res = res, resLFC = resLFC))
}


## ------------------------------------------------------------
## 5. Run the three comparisons
## ------------------------------------------------------------

res_VVvsCC <- run_contrast(dds, "VV", "CC", "Chronic exposure (V1V vs A1A)")
res_CVvsCC <- run_contrast(dds, "CV", "CC", "Acute exposure (A1V vs A1A)")
res_CVvsVV <- run_contrast(dds, "CV", "VV", "Genotype-of-origin (A1V vs V1V)")
