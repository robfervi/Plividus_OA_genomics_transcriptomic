Analysis code for: Natural pH gradients uncover genomic and transcriptomic bases of ocean acidification adaptation in the marine calcifier Paracentrotus lividus

Fernández-Vilert R, Arranz V, Martín-Huete M, Oleksiak MF, Crawford DL, González-Delgado S, Hernández JC, Pérez-Portela R.

Correspondence: r.fvilert@gmail.com

Overview
This repository contains the analysis code used in the study:

Population genomics
ddRADseq-based analyses were used to investigate population structure, detect candidate SNPs associated with pH-driven selection, and functionally annotate genes containing candidate variants.

Transcriptomics
RNA-seq analyses were conducted to investigate differential gene expression and gene co-expression patterns in individuals from ambient (A1) and vent (V1) populations under chronic (natural) and acute (experimental) low-pH exposure.

Full methodological details are provided in the Materials and Methods section of the manuscript.

Scripts are numbered according to the order in which they should be run within each analysis section. The header of each script specifies the required input files and generated outputs.

Population genomics
Processing and filtering
01_ddRAD_pipeline.sh
ddRADseq processing pipeline:
process_radtags → BWA-MEM alignment to the reference genome → gstacks → populations → PLINK missing-data filtering → VCFtools/PLINK filtering and preparation for downstream analyses

02_coverage_HWE_filtering.R
SNP filtering based on:
Sequencing depth using median ± MAD thresholds
Linkage disequilibrium using VCFtools (r²)
Hardy–Weinberg equilibrium using FDR-corrected tests
Produces the final locus whitelist used for downstream analyses.

Analyses
01_RDA_selection.R
Redundancy analysis (RDA) to identify candidate SNPs associated with pH variation. Candidate loci were identified using a |z|-score threshold of 2.5.

02_PCA.R
Principal component analysis (PCA) of the neutral and candidate-SNP datasets.

03_FST.R
Pairwise F<sub>ST</sub> estimation among sampling sites using StAMPP, comparing neutral and candidate SNP datasets.

04_HoHeFis.R
Estimation of observed heterozygosity (H<sub>O</sub>), expected heterozygosity (H<sub>E</sub>), and the inbreeding coefficient (F<sub>IS</sub>) per sampling site and locus.

05_STRUCTURE_pipeline.sh
Bayesian population clustering using STRUCTURE:
K = 1–6
10 replicates per K
Parallelized execution

06_topGO_genpop.R
Gene Ontology (GO) enrichment analysis of genes containing candidate SNPs associated with selection.

Transcriptomics
RNA-seq preprocessing
The RNA-seq preprocessing pipeline, including quality control, trimming, contaminant filtering, alignment, and read quantification, is not included in this repository.
See the Data availability section below for access to the raw sequencing data and preprocessing pipeline.

Analyses
01_DESeq2.R
Differential gene expression analysis using DESeq2.

Three pairwise comparisons were performed:
Chronic exposure: V1V vs. A1A
Acute exposure: A1V vs. A1A
Genotype of origin: A1V vs. V1V

02_topGO_RNA.R
GO enrichment analysis of up- and down-regulated differentially expressed genes (DEGs) from each comparison.

03_WGCNA.R
Weighted gene co-expression network analysis (WGCNA), including module detection and module–trait correlation analyses.

Data availability
Population genomics
ddRADseq sequencing data: NCBI BioProject PRJNA1466738
Transcriptomics
RNA-seq sequencing data: NCBI BioProject PRJNA1466604
Reference genome
Paracentrotus lividus reference genome: GCA_940671915.1
Marlétaz et al. (2023)
RNA-seq preprocessing pipeline

The RNA-seq preprocessing pipeline includes:
FastQC
SortMeRNA
Trimmomatic
Kraken2
HISAT2
FeatureCounts
Repository: vanearranz/Arbacia_lixula_transcriptomics

Software
Analyses were performed using:
General software
R v4.4.3
Linux high-performance computing cluster
Stacks v2.59
BWA-MEM
Samtools
PLINK v1.90b6.24
VCFtools v0.1.16
STRUCTURE v2.3.4
Key R packages
vegan
adegenet
vcfR
StAMPP
hierfstat
DESeq2
WGCNA
topGO
clusterProfiler

Repository maintenance
This repository is maintained by Robert Fernández-Vilert.
Contact: r.fvilert@gmail.com

Citation
If you use the code or data provided in this repository, please cite the associated manuscript:

Fernández-Vilert R, Arranz V, Martín-Huete M, Oleksiak MF, Crawford DL, González-Delgado S, Hernández JC, Pérez-Portela R. Natural pH gradients uncover genomic and transcriptomic bases of ocean acidification adaptation in the marine calcifier Paracentrotus lividus.

Citation information will be updated upon publication.
