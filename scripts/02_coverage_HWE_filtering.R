## ============================================================
## SNP filtering: coverage (median +/- MAD), linkage disequilibrium,
## and Hardy-Weinberg equilibrium
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Combines the three post-populations filters described in Materials
## and Methods into a single whitelist of loci to keep, consumed by
## step 11 of 10_ddRAD_pipeline.sh.
##
## Input (from 10_ddRAD_pipeline.sh, steps 8-10):
##   - out.ldepth.mean   (VCFtools --site-mean-depth: CHROM, POS, MEAN_DEPTH, VAR_DEPTH)
##   - ld_out.geno.ld    (VCFtools --geno-r2: CHR, POS1, POS2, N_INDV, R^2)
##   - PlinkHW.hwe       (PLINK --hardy: per-site, per-population HWE test)
##
## Output:
##   - snps_to_remove_coverage.txt
##   - snps_to_remove_ld.txt
##   - snps_to_remove_hwe.txt
##   - whitelist_final.txt   (loci kept after all three filters combined)
## ============================================================

library(readr)
library(dplyr)

setwd(".")  # run from POPULATIONS/pass2 (see 10_ddRAD_pipeline.sh)


## ------------------------------------------------------------
## 1. Coverage filter: median +/- k*MAD per locus
## (median - 2.5*MAD to median + 4*MAD, per Materials and Methods)
## ------------------------------------------------------------

depth <- read_tsv("out.ldepth.mean", col_types = cols(
  CHROM = col_character(), POS = col_double(),
  MEAN_DEPTH = col_double(), VAR_DEPTH = col_double()
))

k_lower <- 2.5
k_upper <- 4

m   <- median(depth$MEAN_DEPTH, na.rm = TRUE)
mad_val <- mad(depth$MEAN_DEPTH, constant = 1, na.rm = TRUE)
lower <- m - k_lower * mad_val
upper <- m + k_upper * mad_val

cat("Median coverage:", m, " MAD:", mad_val, "\n")
cat("Keeping loci with", lower, "<= mean depth <=", upper, "\n")

depth_remove <- subset(depth, MEAN_DEPTH < lower | MEAN_DEPTH > upper)
depth_keep   <- subset(depth, MEAN_DEPTH >= lower & MEAN_DEPTH <= upper)
cat("Coverage filter: removing", nrow(depth_remove), "/", nrow(depth), "loci\n")

write.table(paste(depth_remove$CHROM, depth_remove$POS, sep = ":"),
            "snps_to_remove_coverage.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


## ------------------------------------------------------------
## 2. LD filter: for each pair with r^2 above threshold, drop the
## second SNP (greedy pruning - the same logic as PLINK's
## --indep-pairwise, applied here to VCFtools' --geno-r2 output since
## VCFtools itself has no built-in pruning step)
## ------------------------------------------------------------

ld <- read_tsv("ld_out.geno.ld", col_types = cols(
  CHR = col_character(), POS1 = col_double(), POS2 = col_double(),
  N_INDV = col_double(), `R^2` = col_double()
))

r2_threshold <- 0.5
high_ld <- ld %>% filter(`R^2` > r2_threshold) %>% arrange(CHR, POS1, POS2)

ld_remove <- character(0)
for (i in seq_len(nrow(high_ld))) {
  snp1 <- paste(high_ld$CHR[i], high_ld$POS1[i], sep = ":")
  snp2 <- paste(high_ld$CHR[i], high_ld$POS2[i], sep = ":")
  if (!(snp1 %in% ld_remove)) ld_remove <- c(ld_remove, snp2)  # keep the first, drop the second
}
ld_remove <- unique(ld_remove)
cat("LD filter (r^2 >", r2_threshold, "): removing", length(ld_remove), "loci\n")

write.table(ld_remove, "snps_to_remove_ld.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


## ------------------------------------------------------------
## 3. Hardy-Weinberg filter: FDR (Benjamini-Yekutieli)-corrected p-value
## Loci excessively out of HWE (excess observed heterozygosity) in
## a locus's sampling sites are removed.
## ------------------------------------------------------------

hwe <- read.table("PlinkHW.hwe", header = TRUE, stringsAsFactors = FALSE)
hwe <- hwe[order(hwe$P), ]
hwe$FDR <- p.adjust(hwe$P, method = "BY")
write.csv(hwe, "PlinkHWadj.csv", row.names = FALSE, quote = FALSE)

alpha <- 0.05
hwe_remove <- unique(hwe$SNP[hwe$FDR < alpha & hwe$TEST == "ALL"])
cat("HWE filter (FDR/BY <", alpha, "): removing", length(hwe_remove), "loci\n")

write.table(hwe_remove, "snps_to_remove_hwe.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


## ------------------------------------------------------------
## 4. Combine: whitelist = all loci minus the three removal sets
## ------------------------------------------------------------

all_loci <- paste(depth$CHROM, depth$POS, sep = ":")
to_remove <- unique(c(depth_remove %>% mutate(id = paste(CHROM, POS, sep = ":")) %>% pull(id),
                      ld_remove, hwe_remove))

whitelist <- setdiff(all_loci, to_remove)
cat("Final whitelist:", length(whitelist), "/", length(all_loci), "loci retained\n")

write.table(whitelist, "whitelist_final.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
