## ============================================================
## Pairwise FST (Weir & Cockerham 1984) via StAMPP
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Run separately on the neutral SNP dataset and on the candidate
## SNPs under pH selection (RDA outliers, z > 2.5), then combine
## both into a single heatmap (upper triangle = selection,
## lower triangle = neutral) and test whether FST at candidate
## loci correlates with FST at neutral loci across site pairs.
##
## Population assignment is read from popmap.txt, as in the PCA
## and Ho/He/Fis scripts (instead of a hardcoded vector).
##
## Input:
##   - populationsNeutro.recode.vcf      (neutral SNPs)
##   - populationsUS_RDA2_5.recode.vcf   (candidate SNPs under selection)
##   - popmap.txt                        (sample, population)
##
## Output (per dataset, "_neutral" / "_selection" suffix):
##   - pairwise_fst_<dataset>.csv
##   - pairwise_fst_pval_<dataset>.csv
##   - pairwise_fst_bootstraps_<dataset>.csv
##   - pairwise_fst_pval_adj_BH_<dataset>.csv
##
## Combined output:
##   - Heatmap: selection FST (upper triangle) vs neutral FST (lower triangle)
##   - Pearson/Spearman correlation between neutral and selection FST per pair
##
## NOTE: this script was reorganized into functions and could not be
## re-run here to confirm the exact structure (dimnames/orientation)
## of StAMPP's $Pvalues output — check that the upper/lower triangle
## indexing below lines up with your actual population order the
## first time you run it.
## ============================================================

# install.packages(c("StAMPP", "adegenet", "vcfR", "dplyr", "pheatmap", "RColorBrewer", "ggplot2", "reshape2"))  # uncomment if needed

library(StAMPP)
library(adegenet)
library(vcfR)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)


## ------------------------------------------------------------
## Population map (sample -> population)
## ------------------------------------------------------------

popmap <- read.table("popmap.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("sample", "population")


## ------------------------------------------------------------
## Compute pairwise FST (point estimate + bootstrap p-values)
## for one VCF, and write the raw output tables
## ------------------------------------------------------------

compute_fst <- function(vcf_file, popmap, label, nboot = 1000) {

  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  gl <- vcfR2genlight(vcf)
  pop(gl) <- as.factor(popmap$population[match(indNames(gl), popmap$sample)])

  fst_point <- stamppFst(gl, nboots = 1, percent = 95, nclusters = 1)          # FST matrix only
  fst_boot  <- stamppFst(gl, nboots = nboot, percent = 95, nclusters = 4)      # + p-values, bootstraps

  write.csv(fst_boot$Fsts,       paste0("pairwise_fst_", label, ".csv"))
  write.csv(fst_boot$Pvalues,    paste0("pairwise_fst_pval_", label, ".csv"))
  write.csv(fst_boot$Bootstraps, paste0("pairwise_fst_bootstraps_", label, ".csv"))

  list(fst = fst_point, pval = fst_boot$Pvalues)
}

res_neutral   <- compute_fst("populationsNeutro.recode.vcf",    popmap, "neutral")
res_selection <- compute_fst("populationsUS_RDA2_5.recode.vcf", popmap, "selection")


## ------------------------------------------------------------
## BH correction of p-values (lower triangle, excluding diagonal)
## ------------------------------------------------------------

adjust_pvals_bh <- function(pval_matrix) {
  pvals <- as.matrix(pval_matrix)
  lower_vals <- pvals[lower.tri(pvals)]
  lower_vals_adj <- p.adjust(lower_vals, method = "BH")
  pvals_adj <- matrix(NA, nrow(pvals), ncol(pvals), dimnames = dimnames(pvals))
  pvals_adj[lower.tri(pvals_adj)] <- lower_vals_adj
  pvals_adj
}

pval_adj_neutral   <- adjust_pvals_bh(res_neutral$pval)
pval_adj_selection <- adjust_pvals_bh(res_selection$pval)

write.csv(pval_adj_neutral,   "pairwise_fst_pval_adj_BH_neutral.csv")
write.csv(pval_adj_selection, "pairwise_fst_pval_adj_BH_selection.csv")


## ------------------------------------------------------------
## Individual heatmaps (one per dataset)
## ------------------------------------------------------------

plot_fst_heatmap <- function(fst_matrix, title) {
  m <- fst_matrix
  diag(m) <- 0
  m[upper.tri(m)] <- t(m)[upper.tri(m)]  # mirror lower triangle to upper
  pheatmap(m, scale = "row", main = title)
}

plot_fst_heatmap(res_neutral$fst,   "Neutral SNPs")
plot_fst_heatmap(res_selection$fst, "Candidate SNPs under selection")


## ------------------------------------------------------------
## Combined heatmap: selection FST (upper triangle) vs
## neutral FST (lower triangle), built directly from the
## computed matrices (no manually typed values)
## ------------------------------------------------------------

pops <- colnames(res_neutral$fst)  # assumes same population set/order in both datasets
stopifnot(setequal(pops, colnames(res_selection$fst)))

combined_fst <- matrix(NA, nrow = length(pops), ncol = length(pops), dimnames = list(pops, pops))
combined_fst[upper.tri(combined_fst)] <- res_selection$fst[upper.tri(res_selection$fst)]
combined_fst[lower.tri(combined_fst)] <- res_neutral$fst[lower.tri(res_neutral$fst)]

# Significance stars from the BH-adjusted p-values (same triangle placement)
labels_matrix <- matrix("", nrow = length(pops), ncol = length(pops), dimnames = list(pops, pops))
threshold <- 0.05
labels_matrix[upper.tri(labels_matrix)] <- ifelse(
  pval_adj_selection[upper.tri(pval_adj_selection)] < threshold, "*", "")
labels_matrix[lower.tri(labels_matrix)] <- ifelse(
  pval_adj_neutral[lower.tri(pval_adj_neutral)] < threshold, "*", "")

heatmap_long <- melt(combined_fst, na.rm = FALSE)
labels_long  <- melt(labels_matrix, na.rm = FALSE)
colnames(labels_long) <- c("Var1", "Var2", "label")
heatmap_long$label <- labels_long$label

ggplot(heatmap_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 9) +
  scale_fill_gradient2(low = "darkseagreen1", mid = "aquamarine", high = "darkcyan",
                        midpoint = median(combined_fst, na.rm = TRUE), na.value = "white") +
  theme_minimal() +
  labs(x = "", y = "", fill = "Fst",
       title = "Pairwise FST: selection (upper) vs neutral (lower)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_equal()


## ------------------------------------------------------------
## Correlation between neutral and selection FST across site pairs
## ------------------------------------------------------------

vals_selection <- combined_fst[upper.tri(combined_fst)]
vals_neutral   <- combined_fst[lower.tri(combined_fst)]
pairs <- combn(pops, 2, FUN = function(x) paste(x, collapse = "-"))

fst_df <- data.frame(
  pair          = pairs,
  Fst_selection = vals_selection,
  Fst_neutral   = vals_neutral
)

cor_pearson  <- cor.test(fst_df$Fst_selection, fst_df$Fst_neutral, method = "pearson")
cor_spearman <- cor.test(fst_df$Fst_selection, fst_df$Fst_neutral, method = "spearman")
print(cor_pearson)
print(cor_spearman)

ggplot(fst_df, aes(x = Fst_neutral, y = Fst_selection, label = pair)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_text(vjust = -0.7, size = 3) +
  labs(x = "Neutral Fst", y = "Fst under selection",
       title = "Neutral Fst vs Fst under selection, by site pair") +
  theme_minimal()

boxplot(fst_df$Fst_neutral, fst_df$Fst_selection,
        names = c("Neutral", "Under selection"),
        ylab = "Fst", col = c("lightblue", "lightgreen"))

# Pairs with the largest neutral-vs-selection FST difference
fst_df$difference <- abs(fst_df$Fst_selection - fst_df$Fst_neutral)
fst_df[order(-fst_df$difference), ]

diff_threshold <- 0.05
fst_df$highlight <- fst_df$difference > diff_threshold

ggplot(fst_df, aes(x = Fst_neutral, y = Fst_selection, label = pair)) +
  geom_point(aes(color = highlight), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_text(data = subset(fst_df, highlight), vjust = -0.7, size = 3) +
  scale_color_manual(values = c("black", "orange"), guide = "none") +
  labs(x = "Neutral Fst", y = "Fst under selection",
       title = "Neutral vs selection Fst by site pair",
       subtitle = paste("Highlighted pairs: difference >", diff_threshold)) +
  theme_minimal()
