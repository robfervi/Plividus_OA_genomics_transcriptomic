############################################################
## Ho, He and Fis per locus and per population
## Neutral SNPs vs candidate SNPs under pH selection
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
############################################################

############################################################
# 0. LIBRARIES
############################################################
library(vcfR)
library(adegenet)
library(hierfstat)
library(poppr)
library(dplyr)
library(tidyr)
library(ggplot2)

############################################################
# 1. INPUT FILES
############################################################
vcf_neutral_file  <- "neutral.vcf"
vcf_selected_file <- "under_selection.vcf"
popmap_file       <- "popmap.txt"

############################################################
# 2. READ VCF
############################################################
vcf_neutral  <- read.vcfR(vcf_neutral_file)
vcf_selected <- read.vcfR(vcf_selected_file)

############################################################
# 3. CONVERT TO GENIND
############################################################
genind_neutral  <- vcfR2genind(vcf_neutral)
genind_selected <- vcfR2genind(vcf_selected)

############################################################
# 4. READ POPMAP AND ASSIGN POPULATIONS
############################################################
popmap <- read.table(popmap_file, header = FALSE, stringsAsFactors = FALSE)
colnames(popmap) <- c("sample", "population")

pop_neutral  <- popmap$population[match(indNames(genind_neutral), popmap$sample)]
pop_selected <- popmap$population[match(indNames(genind_selected), popmap$sample)]

pop(genind_neutral)  <- as.factor(pop_neutral)
pop(genind_selected) <- as.factor(pop_selected)

############################################################
# 5. CONVERT TO HIERFSTAT
############################################################
hf_neutral  <- genind2hierfstat(genind_neutral)
hf_selected <- genind2hierfstat(genind_selected)

############################################################
# 6. COMPUTE STATISTICS
############################################################
stats_neutral  <- basic.stats(hf_neutral)
stats_selected <- basic.stats(hf_selected)

############################################################
# 7. REMOVE LOCI WITH NA (e.g. dropped automatically if e.g.
#    a locus is monomorphic in one population)
############################################################
valid_loci_neutral  <- colnames(stats_neutral$Hs)[!is.na(colMeans(stats_neutral$Hs))]
valid_loci_selected <- colnames(stats_selected$Hs)[!is.na(colMeans(stats_selected$Hs))]

hf_neutral  <- hf_neutral[, c(TRUE, colnames(hf_neutral)[-1] %in% valid_loci_neutral), drop = FALSE]
hf_selected <- hf_selected[, c(TRUE, colnames(hf_selected)[-1] %in% valid_loci_selected), drop = FALSE]

# recompute stats without the problematic loci
stats_neutral  <- basic.stats(hf_neutral)
stats_selected <- basic.stats(hf_selected)

############################################################
# 8. Fis + p-value FUNCTION, PER POPULATION
############################################################
calc_fis_pvalue_pop <- function(hf, nboot = 1000) {

  loci <- colnames(hf)[-1]
  pops <- unique(hf[, 1])

  results <- expand.grid(
    Locus = loci,
    Population = pops
  )

  results$Fis <- NA
  results$p_value <- NA

  for (i in seq_along(loci)) {
    for (j in seq_along(pops)) {

      pop_i <- pops[j]
      sub <- hf[hf[, 1] == pop_i, c(1, i + 1)]
      sub <- sub[!is.na(sub[, 2]), ]

      if (nrow(sub) < 5) next

      bs <- basic.stats(sub)
      fis_obs <- bs$Fis

      if (is.na(fis_obs)) next

      fis_boot <- replicate(nboot, {
        samp <- sub[sample(1:nrow(sub), replace = TRUE), ]
        bs_boot <- basic.stats(samp)
        bs_boot$Fis
      })

      fis_boot <- fis_boot[!is.na(fis_boot)]
      if (length(fis_boot) == 0) next

      p_val <- mean(abs(fis_boot) >= abs(fis_obs))

      idx <- which(results$Locus == loci[i] & results$Population == pop_i)

      results$Fis[idx] <- fis_obs
      results$p_value[idx] <- p_val
    }
  }

  return(results)
}

############################################################
# 9. COMPUTE Fis + p-values
############################################################
fis_pop_neutral  <- calc_fis_pvalue_pop(hf_neutral)
fis_pop_selected <- calc_fis_pvalue_pop(hf_selected)

############################################################
# 10. PER-LOCUS RESULTS (global)
############################################################
locus_results_neutral <- data.frame(
  Locus = colnames(stats_neutral$Ho),
  Ho = colMeans(stats_neutral$Ho, na.rm = TRUE),
  He = colMeans(stats_neutral$Hs, na.rm = TRUE),
  Fis = colMeans(stats_neutral$Fis, na.rm = TRUE),
  Dataset = "Neutral"
)

locus_results_selected <- data.frame(
  Locus = colnames(stats_selected$Ho),
  Ho = colMeans(stats_selected$Ho, na.rm = TRUE),
  He = colMeans(stats_selected$Hs, na.rm = TRUE),
  Fis = colMeans(stats_selected$Fis, na.rm = TRUE),
  Dataset = "Selected"
)

locus_results <- rbind(locus_results_neutral, locus_results_selected)

############################################################
# 11. FDR CORRECTION FOR PER-POPULATION p-values
############################################################
fis_pop_neutral$p_adj  <- p.adjust(fis_pop_neutral$p_value, method = "fdr")
fis_pop_selected$p_adj <- p.adjust(fis_pop_selected$p_value, method = "fdr")

############################################################
# 12. GLOBAL MEANS
############################################################
results_summary <- data.frame(
  Dataset = c("Neutral", "Selected"),
  Ho = c(mean(stats_neutral$Ho, na.rm = TRUE),
         mean(stats_selected$Ho, na.rm = TRUE)),
  He = c(mean(stats_neutral$Hs, na.rm = TRUE),
         mean(stats_selected$Hs, na.rm = TRUE)),
  Fis = c(mean(stats_neutral$Fis, na.rm = TRUE),
          mean(stats_selected$Fis, na.rm = TRUE))
)

print(results_summary)

############################################################
# 13. EXPORT
############################################################
write.csv(results_summary, "summary_Ho_He_Fis.csv", row.names = FALSE)
write.csv(locus_results, "locus_Ho_He_Fis.csv", row.names = FALSE)
write.csv(fis_pop_neutral, "fis_pop_neutral.csv", row.names = FALSE)
write.csv(fis_pop_selected, "fis_pop_selected.csv", row.names = FALSE)

############################################################
# 14. PLOT
############################################################
results_long <- results_summary %>%
  pivot_longer(cols = c(Ho, He, Fis), names_to = "Metric", values_to = "Value")

ggplot(results_long, aes(x = Dataset, y = Value, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Metric, scales = "free") +
  theme_minimal() +
  theme(text = element_text(size = 14))
