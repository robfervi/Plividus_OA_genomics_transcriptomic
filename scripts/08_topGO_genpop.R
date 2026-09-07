## ============================================================
## GO enrichment (topGO) for candidate genes under pH selection
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Genes containing candidate SNPs under selection (from the RDA,
## 01_RDA_selection.R) as the test set, against all annotated genes
## in the reference genome as background.
##
## Reporting follows Materials and Methods: raw elim-Fisher p-values
## (< 0.05), not FDR-adjusted (see 07_topGO_RNA.R header for why).
##
## Input:
##   - genes_ID.txt      (candidate genes containing SNPs under selection)
##   - all_genes.txt     (background: all annotated genes)
##   - association2.txt  (gene -> GO annotation, gene<TAB>GO;GO;GO...)
##
## Output:
##   - TopGO_enrichment_results_all_genpop.csv
##   - GO_enrichment_SIGresults_elim_top15_genpop.csv
##   - bubble_plot_gogenpop.pdf
## ============================================================

# BiocManager::install(c("topGO", "GO.db"))  # uncomment if needed

library(topGO)
library(GO.db)
library(dplyr)
library(ggplot2)

## ------------------------------------------------------------
## 1. Gene lists and GO annotation
## ------------------------------------------------------------

degs_genpop <- read.table("genes_ID.txt", header = FALSE, stringsAsFactors = FALSE)$V1
all_genes   <- read.table("all_genes.txt", header = FALSE, stringsAsFactors = FALSE)$V1

gene_list <- factor(as.integer(all_genes %in% degs_genpop))
names(gene_list) <- all_genes

annotations <- read.table("association2.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
colnames(annotations) <- c("Gene", "GO")
annotations <- annotations[annotations$GO != "-", ]
annotations_expanded <- annotations %>% tidyr::separate_rows(GO, sep = ";")
gene2GO <- split(annotations_expanded$GO, annotations_expanded$Gene)


## ------------------------------------------------------------
## 2. Run all 4 Fisher-test algorithms for one ontology
## ------------------------------------------------------------

run_topGO_ontology <- function(gene_list, ontology, category_name) {
  GOdata <- new("topGOdata", description = "GO analysis", ontology = ontology,
                allGenes = gene_list, geneSelectionFun = function(x) x == 1,
                annot = annFUN.gene2GO, gene2GO = gene2GO)

  allGO <- usedGO(GOdata)
  resultFisher      <- runTest(GOdata, algorithm = "classic",     statistic = "fisher")
  result_elim       <- runTest(GOdata, algorithm = "elim",        statistic = "fisher")
  result_weight     <- runTest(GOdata, algorithm = "weight",      statistic = "fisher")
  result_parentchild <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")

  GenTable(GOdata,
           classicFisher = resultFisher,
           elimFisher = result_elim,
           weightFisher = result_weight,
           parentchildFisher = result_parentchild,
           orderBy = "classicFisher",
           topNodes = length(allGO)) %>%
    mutate(Category = category_name,
           across(c(classicFisher, elimFisher, weightFisher, parentchildFisher), as.numeric))
}

table_BP <- run_topGO_ontology(gene_list, "BP", "Biological Process")
table_MF <- run_topGO_ontology(gene_list, "MF", "Molecular Function")
table_CC <- run_topGO_ontology(gene_list, "CC", "Cellular Component")

combined_results <- bind_rows(table_BP, table_MF, table_CC)
write.csv(combined_results, "TopGO_enrichment_results_all_genpop.csv", row.names = FALSE)


## ------------------------------------------------------------
## 3. elim-Fisher significant terms, top 15 per category, with full names
## ------------------------------------------------------------

sig_elim <- combined_results %>%
  filter(elimFisher < 0.05) %>%
  group_by(Category) %>%
  slice_min(elimFisher, n = 15) %>%
  ungroup()
sig_elim$Term <- Term(GOTERM[sig_elim$GO.ID])

write.csv(sig_elim, "GO_enrichment_SIGresults_elim_top15_genpop.csv", row.names = FALSE)


## ------------------------------------------------------------
## 4. Bubble plot
## ------------------------------------------------------------

go_data_filtered <- sig_elim %>%
  mutate(Category = factor(Category, levels = c("Biological Process", "Molecular Function", "Cellular Component"))) %>%
  arrange(Category, elimFisher) %>%
  mutate(Term_ordered = factor(Term, levels = rev(Term)))

p <- ggplot(go_data_filtered, aes(y = Term_ordered, x = elimFisher, size = Significant, fill = Category)) +
  geom_point(shape = 21, color = "black", stroke = 0.7) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "azure4", linewidth = 0.5) +
  scale_fill_manual(values = c("Biological Process" = "coral2",
                                "Molecular Function" = "cornflowerblue",
                                "Cellular Component" = "#2ca02c")) +
  scale_x_continuous(breaks = c(0.01, 0.03, 0.05), limits = c(0, 0.05)) +
  labs(title = "Enriched GO terms - candidate genes under pH selection",
       x = "P-value (elimFisher)", y = "", size = "Significant genes", fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

ggsave("bubble_plot_gogenpop.pdf", plot = p, width = 18, height = 15, units = "cm", bg = "transparent")
