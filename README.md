Natural pH gradients uncover genomic and transcriptomic bases of ocean acidification adaptation in the marine calcifier Paracentrotus lividus

Analysis code for the manuscript:

Fernández-Vilert R, Arranz V, Martín-Huete M, Oleksiak MF, Crawford DL, González-Delgado S, Hernández JC, Pérez-Portela R. Natural pH gradients uncover genomic and transcriptomic bases of ocean acidification adaptation in the marine calcifier Paracentrotus lividus.

Correspondence: r.fvilert@gmail.com

Abstract

Ocean acidification (OA) is reshaping marine environments at an unprecedented pace, yet whether and how populations can evolve tolerance — and through which molecular mechanisms — remains a central open question for predicting species persistence under future climate scenarios. Natural CO₂ vent systems, where pH is chronically reduced by volcanic emissions, provide unique opportunities to study adaptive divergence under conditions that mirror future oceans. Here, we integrate population genomics, transcriptomics, and morphological data to dissect the molecular bases of adaptation to OA in Paracentrotus lividus, a marine calcifier inhabiting the Fuencaliente CO₂ vent system (La Palma, Canary Islands). Using ddRADseq on 95 individuals sampled along a natural pH gradient (7.5–8.1), we identified 27,542 SNPs, including 232 candidate loci under pH-driven selection. Despite negligible differentiation at neutral loci — confirming high gene flow — loci under selection revealed strong divergence across the gradient, with enrichment in genes related to ion transport, biomineralization, and acid–base regulation. Transcriptomic analyses revealed pronounced constitutive divergence between ambient (pH 8.1) and vent (pH 7.5) populations, consistent with long-term molecular reprogramming of metabolic pathways. In contrast, acute experimental exposure of ambient individuals to low pH triggered a modest stress response centered on protein stabilization that did not converge with the constitutive transcriptional profile of vent-adapted individuals, underscoring the distinction between transient plasticity and evolved tolerance. Together, our results show that adaptive molecular divergence can emerge at fine spatial scales despite extensive gene flow, and that resilience to OA is driven mainly by evolved molecular reprogramming rather than short-term plasticity.

Keywords: Echinodermata, Sea urchin, Evolutionary genomics, Gene expression, Natural selection, CO2 seeps

Overview

This repository contains the analysis code used in the study:

Population genomics — ddRADseq-based detection of candidate SNPs under pH-driven selection, population structure, and functional annotation of candidate genes. 111 individuals were sampled across 6 sites (A1, A2, T1, T2, V1, V2) spanning ambient, transition, and vent pH conditions; 95 individuals and 27,542 SNPs (232 under candidate selection) passed filtering.
Transcriptomics — RNA-seq differential expression and co-expression network analysis of individuals from ambient (A1) and vent (V1) sites, under chronic (natural) and acute (experimental) low-pH exposure.

Full methodological details are in the Materials and Methods of the manuscript.

Repository structure

Scripts are numbered in the order they're meant to run. Each one's header comment lists its exact inputs and outputs.

Population genomics
Script	Description
10_ddRAD_pipeline.sh	ddRADseq processing: process_radtags → BWA-MEM alignment to the reference genome → gstacks → populations → PLINK (missing-data filter) → VCFtools/PLINK (coverage, LD, HWE inputs)
11_coverage_HWE_filtering.R	SNP filtering: sequencing depth (median ± MAD), linkage disequilibrium (VCFtools r²), Hardy–Weinberg equilibrium (FDR-corrected) — produces the final locus whitelist
01_RDA_selection.R	Redundancy analysis (RDA) to detect candidate SNPs under pH selection (
02_PCA.R	Principal component analysis (neutral and candidate-SNP datasets)
03_FST.R	Pairwise F<sub>ST</sub> between sampling sites (StAMPP), neutral vs. candidate SNPs
04_HoHeFis.R	Observed/expected heterozygosity and F<sub>IS</sub> per site and locus
05_STRUCTURE_pipeline.sh	Bayesian clustering with STRUCTURE (K = 1–6, 10 replicates/K, parallelized)
08_topGO_genpop.R	GO enrichment of genes containing candidate SNPs under selection
Transcriptomics
Script	Description
06_DESeq2.R	Differential gene expression (DESeq2), three pairwise comparisons: chronic exposure (V1V vs. A1A), acute exposure (A1V vs. A1A), genotype-of-origin (A1V vs. V1V)
07_topGO_RNA.R	GO enrichment of the up/down-regulated DEGs from each comparison
09_WGCNA.R	Weighted gene co-expression network analysis, module–trait correlations

RNA-seq preprocessing (quality control, trimming, contaminant filtering, alignment, read quantification) is not included here — see Data availability below.

Data availability
ddRADseq raw reads (BAM): BioProject PRJNA1466738
RNA-seq raw reads: BioProject PRJNA1466604
Reference genome: GCA_940671915.1 (Marlétaz et al. 2023)
RNA-seq preprocessing pipeline (FastQC, SortMeRNA, Trimmomatic, Kraken2, HISAT2, FeatureCounts): vanearranz/Arbacia_lixula_transcriptomics
Software

Analyses were run in R (v4.4.3) and on a Linux cluster running Stacks v2.59, BWA-MEM, Samtools, PLINK v1.90b6.24, VCFtools v0.1.16, and STRUCTURE v2.3.4. Key R packages: vegan, adegenet, vcfR, StAMPP, hierfstat, DESeq2, WGCNA, topGO, clusterProfiler.

Authors

Robert Fernández-Vilert¹˒², Vanessa Arranz¹˒², Marta Martín-Huete¹˒², Marjorie F. Oleksiak³, Douglas L. Crawford³, Sara González-Delgado¹˒², José Carlos Hernández⁴, Rocío Pérez-Portela¹˒² (senior author)

¹ Dept. de Biologia Evolutiva, Ecologia i Ciències Ambientals, Universitat de Barcelona, Spain ² Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona, Spain ³ Marine Biology and Ecology, Rosenstiel School of Marine, Atmospheric and Earth Science, University of Miami, USA ⁴ Dept. Biología Animal, Edafología y Geología, Universidad de La Laguna, Tenerife, Spain

Repository maintained by Robert Fernández-Vilert (r.fvilert@gmail.com).
