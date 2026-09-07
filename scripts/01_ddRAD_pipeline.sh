#!/bin/bash
## ============================================================
## ddRADseq processing pipeline (reference-based) - Stacks + PLINK + VCFtools
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma), 111 individuals,
## 6 sampling sites (A1, A2, T1, T2, V1, V2)
##
## Steps (see Materials and Methods):
##   1. process_radtags   - demultiplex + quality trim
##   2. bwa mem            - align each sample to the reference genome
##   3. samtools            - convert/sort/index BAM, check mapping quality
##   4. gstacks              - assemble loci, call variants, genotype
##   5. populations (pass 1) - initial VCF (-r 0.8, --min-maf 0.03, --write-single-snp)
##   6. PLINK                - remove individuals with >15% missing genotypes
##   7. populations (pass 2) - re-run on the filtered sample set
##   8. VCFtools              - per-locus mean depth (input to 11_coverage_HWE_filtering.R)
##   9. VCFtools              - LD (r^2) between SNP pairs (input to 11_coverage_HWE_filtering.R)
##  10. PLINK                - Hardy-Weinberg exact test p-values (input to 11_coverage_HWE_filtering.R)
##  11. populations (pass 3) - final VCF, using the whitelist produced by the R script
##
## Reference genome: GCA_940671915.1 (Marlétaz et al., 2023)
## ============================================================

set -euo pipefail

## ------------------------------------------------------------
## CONFIG - edit before running
## ------------------------------------------------------------

DEMULTIPLEXED_DIR=./DEMULTIPLEXED
BARCODES_FILE=./SCRIPTS/BARCODES/barcodes_Pliv.txt
REFERENCE=./reference/GCA_940671915.1_reference_genome.fna.gz
ALIGN_OUT=./ALIGN_OUT
GSTACKS_OUT=./GSTACKS_OUT
POPMAP=./SCRIPTS/BARCODES/popmap.txt          # sample <TAB> site (A1/A2/T1/T2/V1/V2)
POPMAP_FILTERED=./SCRIPTS/BARCODES/popmap_filtered.txt  # after step 6, individuals with <=15% missing
POP_OUT_1=./POPULATIONS/pass1
POP_OUT_2=./POPULATIONS/pass2
POP_OUT_FINAL=./POPULATIONS/final
THREADS=20


## ------------------------------------------------------------
## 1. process_radtags - demultiplex + quality trim
## ------------------------------------------------------------

mkdir -p "$DEMULTIPLEXED_DIR"

process_radtags -i gzfastq -P \
  -p ./LIBRARIES/CONCAT_LIBS/ \
  -o "$DEMULTIPLEXED_DIR" \
  -b "$BARCODES_FILE" \
  -c -q -r --inline_index \
  --renz_1 ecoRI --renz_2 mseI \
  --adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG \
  --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  --adapter_mm 2

# Check ./process_radtags.log for retained-read %, barcode/adapter/cutsite drop rates.
# If a large fraction is dropped at "RAD cutsite not found", check whether the reverse
# reads carry an extra leading base before the MseI site - see PIPELINE_STACKS notes
# on --disable_rad_check / trimming the leading base with cutadapt.


## ------------------------------------------------------------
## 2-3. Align each sample to the reference genome (BWA-MEM + Samtools)
## ------------------------------------------------------------

mkdir -p "$ALIGN_OUT" reference

bwa index "$REFERENCE"

for file in $(cat "$DEMULTIPLEXED_DIR/samples.txt"); do
  bwa mem -t $THREADS "$REFERENCE" \
    "$DEMULTIPLEXED_DIR/${file}.1.fq.gz" "$DEMULTIPLEXED_DIR/${file}.2.fq.gz" \
    | samtools view -b -@ $THREADS \
    | samtools sort -@ $THREADS > "$ALIGN_OUT/${file}.bam"
  samtools index "$ALIGN_OUT/${file}.bam"
done

# Mapping quality check (% mapped/properly paired reads per sample)
mkdir -p "$ALIGN_OUT/bam_statistics"
for bam in "$ALIGN_OUT"/*.bam; do
  samtools stats "$bam" > "$ALIGN_OUT/bam_statistics/$(basename "$bam")_stats"
done


## ------------------------------------------------------------
## 4. gstacks - assemble loci, call variants, genotype
## ------------------------------------------------------------

mkdir -p "$GSTACKS_OUT"
gstacks -I "$ALIGN_OUT" -M "$POPMAP" -O "$GSTACKS_OUT" -t $THREADS

# Check effective per-sample coverage in gstacks.log.distribs:
#   stacks-dist-extract "$GSTACKS_OUT/gstacks.log.distribs" effective_coverages_per_sample
# Note: --rm-pcr-duplicates is NOT used for ddRAD (double-digest data) - the two
# cut sites break the duplicate-detection algorithm (see Stacks documentation).


## ------------------------------------------------------------
## 5. populations (pass 1) - initial SNP set for QC/filtering
## ------------------------------------------------------------

mkdir -p "$POP_OUT_1"
populations -P "$GSTACKS_OUT" -M "$POPMAP" -O "$POP_OUT_1" \
  -r 0.8 --min-maf 0.03 --write-single-snp \
  --vcf --genepop --structure --plink


## ------------------------------------------------------------
## 6. PLINK - remove individuals with >15% missing genotypes
## ------------------------------------------------------------

cd "$POP_OUT_1"
plink --vcf populations.snps.vcf --allow-extra-chr --double-id --mind 0.15 \
  --make-just-fam --out missingness_filter
# missingness_filter.fam lists the individuals PLINK kept (--mind 0.15 removes
# individuals with >15% missing genotypes). Build the filtered popmap from it:
awk '{print $1}' missingness_filter.fam > kept_individuals.txt
cd -

grep -Ff "$POP_OUT_1/kept_individuals.txt" "$POPMAP" > "$POPMAP_FILTERED"


## ------------------------------------------------------------
## 7. populations (pass 2) - re-run on the filtered sample set
## ------------------------------------------------------------

mkdir -p "$POP_OUT_2"
populations -P "$GSTACKS_OUT" -M "$POPMAP_FILTERED" -O "$POP_OUT_2" \
  -r 0.8 --min-maf 0.03 --write-single-snp \
  --vcf --genepop --structure --plink


## ------------------------------------------------------------
## 8. VCFtools - per-locus mean depth (-> 11_coverage_HWE_filtering.R)
## ------------------------------------------------------------

cd "$POP_OUT_2"
vcftools --vcf populations.snps.vcf --site-mean-depth --out out
# produces out.ldepth.mean (CHROM, POS, MEAN_DEPTH, VAR_DEPTH)


## ------------------------------------------------------------
## 9. VCFtools - pairwise LD (r^2) between SNPs (-> 11_coverage_HWE_filtering.R)
## ------------------------------------------------------------

vcftools --vcf populations.snps.vcf --geno-r2 --ld-window-bp 50000 --out ld_out
# produces ld_out.geno.ld (CHR, POS1, POS2, N_INDV, R^2)
cd -


## ------------------------------------------------------------
## 10. PLINK - Hardy-Weinberg exact test p-values (-> 11_coverage_HWE_filtering.R)
## ------------------------------------------------------------

cd "$POP_OUT_2"
plink --vcf populations.snps.vcf --allow-extra-chr --double-id --hardy --out PlinkHW
# produces PlinkHW.hwe (with per-site, per-population HWE p-values)
cd -

echo "Now run 11_coverage_HWE_filtering.R on:"
echo "  - $POP_OUT_2/out.ldepth.mean"
echo "  - $POP_OUT_2/ld_out.geno.ld"
echo "  - $POP_OUT_2/PlinkHW.hwe"
echo "It writes a combined whitelist of loci to keep: whitelist_final.txt"


## ------------------------------------------------------------
## 11. populations (final) - using the whitelist from the R filtering script
## ------------------------------------------------------------

mkdir -p "$POP_OUT_FINAL"
populations -P "$GSTACKS_OUT" -M "$POPMAP_FILTERED" -O "$POP_OUT_FINAL" \
  --whitelist "$POP_OUT_2/whitelist_final.txt" \
  --vcf --genepop --structure --plink
