## ============================================================
## RDA (Redundancy Analysis) for detecting candidate SNPs
## under pH-driven selection
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Input:
##   - populations.snps.vcf   (filtered VCF, output of Stacks/populations)
##   - env_matrix_RDA.csv     (environmental matrix: sample + pH)
##
## Output:
##   - SNPid_chr_pos.txt        Position info (ID, chromosome, bp) for each SNP
##   - geno_matrix_chr.txt      Encoded genotype matrix (0/1/2)
##   - outliers_rda_z2_5.txt    Candidate SNPs under selection (|z| > 2.5)
##
## Method: Forester et al. 2018 (multivariate RDA). Threshold z > 2.5
## (see Materials and Methods).
## ============================================================

# install.packages(c("vegan", "adegenet", "vcfR", "readr"))  # uncomment if needed

library(vegan)
library(adegenet)
library(vcfR)
library(readr)


## ------------------------------------------------------------
## 1. Load VCF and extract SNP position information
## ------------------------------------------------------------

obj.vcfR <- read.vcfR("populations.snps.vcf")

position   <- getPOS(obj.vcfR)    # position in bp
chromosome <- getCHROM(obj.vcfR)  # chromosome / scaffold
id_snp     <- getID(obj.vcfR)     # SNP ID

chr_pos <- as.data.frame(cbind(id_snp, chromosome, position))
chr_pos$position <- as.numeric(as.character(chr_pos$position))

write.table(chr_pos, "SNPid_chr_pos.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## ------------------------------------------------------------
## 2. Build the genotype matrix (0 = homozygous ref,
##    1 = heterozygous, 2 = homozygous alt, 9 = missing)
## ------------------------------------------------------------

geno <- extract.gt(obj.vcfR)

G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

write.table(G, "geno_matrix_chr.txt", sep = "\t", col.names = FALSE, row.names = FALSE)
dim(G)

# Transpose (samples as rows, SNPs as columns) and name columns from chr_pos
gen <- t(G)
colnames(gen) <- paste(chr_pos$chromosome, chr_pos$position, sep = "_")
gen[1:10, 1:10]  # quick check


## ------------------------------------------------------------
## 3. Missing data and imputation
## ------------------------------------------------------------

gen[which(gen == "9")] <- NA
sum(is.na(gen)) / (dim(gen)[1] * dim(gen)[2])  # % missing data

# Impute with the most common genotype at that position
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp))  # sanity check: should be 0


## ------------------------------------------------------------
## 4. Environmental matrix (pH per sample)
## ------------------------------------------------------------

env_matrix <- read_delim(
  "env_matrix_RDA.csv",
  delim = ";",
  escape_double = FALSE,
  col_types = cols(pH = col_number()),
  trim_ws = TRUE
)

info <- env_matrix
head(info)


## ------------------------------------------------------------
## 5. RDA: genotypes ~ pH
## ------------------------------------------------------------

Env.rda <- rda(gen.imp ~ info$pH, scale = TRUE)
Env.rda
plot(Env.rda, scaling = 3)  # circles = samples, red crosses = SNPs, vector = pH along RDA1

load.Env.rda <- scores(Env.rda, choices = c(1), display = "species")
load.Env.rda.pos <- cbind(chr_pos, load.Env.rda)
head(load.Env.rda.pos)


## ------------------------------------------------------------
## 6. Outlier SNP detection (candidates under selection)
## ------------------------------------------------------------

# Threshold: z > 2.5 standard deviations (see Materials and Methods)
# z = 3 would correspond to p = 0.0027; z = 2.5 (used here) corresponds to p = 0.001
z <- 2.5

lim_min <- mean(load.Env.rda.pos$RDA1) - (z * sd(load.Env.rda.pos$RDA1))
lim_max <- mean(load.Env.rda.pos$RDA1) + (z * sd(load.Env.rda.pos$RDA1))

outlier_Env <- load.Env.rda.pos[load.Env.rda.pos$RDA1 >= lim_max | load.Env.rda.pos$RDA1 <= lim_min, ]
outlier_Env

write.table(outlier_Env, "outliers_rda_z2_5.txt", row.names = FALSE, quote = FALSE, sep = "\t")
