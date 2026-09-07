## ============================================================
## GO enrichment (topGO) for the transcriptomic DEGs
##
## P. lividus - Fuencaliente CO2 vent gradient (La Palma)
##
## Runs GO enrichment (BP/MF/CC, four Fisher-test algorithms) on the
## up- and down-regulated gene sets from each of the three DESeq2
## comparisons (06_DESeq2.R), then builds a combined 6-panel bubble
## plot of the elim-Fisher significant terms.
##
## Reporting follows the Materials and Methods: results are reported
## as raw elim-Fisher p-values (< 0.05), not FDR-adjusted, because the
## elim algorithm's term-by-term dependence violates the independence
## assumption behind standard multiple-testing correction.
##
## Input (from 06_DESeq2.R):
##   - all_genes.txt
##   - DEG_<contrast>_up.txt / DEG_<contrast>_down.txt   for VVvsCC, CVvsCC, CVvsVV
##   - association2.txt   (gene -> GO annotation, gene<TAB>GO;GO;GO...)
##
## Output:
##   - GO_enrichment_results_all_<contrast>_<direction>.csv       (all terms, all methods)
##   - GO_enrichment_SIGresults_elim_top15_<contrast>_<direction>.csv
##   - GO_BubblePlots_6panel_elim.pdf
## ============================================================

# BiocManager::install(c("topGO", "clusterProfiler", "GO.db"))  # uncomment if needed

library(dplyr)
library(tidyr)
library(ggplot2)
library(topGO)
library(GO.db)
library(patchwork)

## ------------------------------------------------------------
## 1. Load the gene -> GO annotation once (shared across contrasts)
## ------------------------------------------------------------

annotations <- read.table("association2.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
colnames(annotations) <- c("Gene", "GO")
annotations <- annotations[annotations$GO != "-", ]
annotations_expanded <- annotations %>% separate_rows(GO, sep = ";")
gene2GO <- split(annotations_expanded$GO, annotations_expanded$Gene)

all_genes <- read.table("all_genes.txt", header = FALSE, stringsAsFactors = FALSE)$V1


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


## ------------------------------------------------------------
## 3. Full enrichment for one DEG set (up or down of one contrast)
## ------------------------------------------------------------

run_topGO_analysis <- function(degs_file, label) {

  degs <- read.table(degs_file, header = FALSE, stringsAsFactors = FALSE)$V1
  gene_list <- factor(as.integer(all_genes %in% degs))
  names(gene_list) <- all_genes

  table_BP <- run_topGO_ontology(gene_list, "BP", "Biological Process")
  table_MF <- run_topGO_ontology(gene_list, "MF", "Molecular Function")
  table_CC <- run_topGO_ontology(gene_list, "CC", "Cellular Component")

  combined <- bind_rows(table_BP, table_MF, table_CC)
  write.csv(combined, paste0("GO_enrichment_results_all_", label, ".csv"), row.names = FALSE)

  # elim-Fisher significant terms, top 15 per category, with full GO term names
  sig_elim <- combined %>%
    filter(elimFisher < 0.05) %>%
    group_by(Category) %>%
    slice_min(elimFisher, n = 15) %>%
    ungroup()
  sig_elim$Term <- Term(GOTERM[sig_elim$GO.ID])

  write.csv(sig_elim, paste0("GO_enrichment_SIGresults_elim_top15_", label, ".csv"), row.names = FALSE)
  sig_elim
}


## ------------------------------------------------------------
## 4. Run for the 3 comparisons x 2 directions (6 gene sets)
## ------------------------------------------------------------

contrasts <- c("VVvsCC", "CVvsCC", "CVvsVV")
directions <- c("up", "down")

go_results <- list()
for (contrast in contrasts) {
  for (direction in directions) {
    label <- paste0(contrast, "_", direction)
    degs_file <- paste0("DEG_", contrast, "_", direction, ".txt")
    go_results[[label]] <- run_topGO_analysis(degs_file, label)
  }
}


## ------------------------------------------------------------
## 5. Combined 6-panel bubble plot (elim-Fisher significant terms)
## ------------------------------------------------------------

prep_go <- function(df) {
  df %>%
    mutate(Category = factor(Category, levels = c("Biological Process", "Molecular Function", "Cellular Component"))) %>%
    arrange(Category, elimFisher) %>%
    mutate(Term_ordered = factor(Term, levels = rev(Term)))
}

bubble_plot_elim <- function(df, x_limit = 0.05, y_right = FALSE, legend_pos = "right") {
  ggplot(df, aes(y = Term_ordered, x = elimFisher, size = Significant, fill = Category)) +
    geom_point(shape = 21, color = "black", stroke = 0.7) +
    geom_vline(xintercept = 0.01, linetype = "dashed", color = "azure4", linewidth = 0.5) +
    scale_fill_manual(values = c("Biological Process" = "coral2",
                                  "Molecular Function" = "cornflowerblue",
                                  "Cellular Component" = "#2ca02c")) +
    scale_x_continuous(breaks = c(0.01, 0.03, 0.05), limits = c(0, x_limit)) +
    scale_y_discrete(position = ifelse(y_right, "right", "left")) +
    labs(x = "", y = "", size = "Significant genes", fill = "Category") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 18),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.position = legend_pos)
}

panels <- lapply(names(go_results), function(label) prep_go(go_results[[label]]))
names(panels) <- names(go_results)

rows <- lapply(contrasts, function(contrast) {
  p_up   <- bubble_plot_elim(panels[[paste0(contrast, "_up")]],   y_right = FALSE, legend_pos = "left")
  p_down <- bubble_plot_elim(panels[[paste0(contrast, "_down")]], y_right = TRUE,  legend_pos = "right")
  p_up + p_down
})

final_plot <- (rows[[1]] / rows[[2]] / rows[[3]]) + plot_layout(heights = c(5, 2, 3))
final_plot <- final_plot & theme(
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background  = element_rect(fill = "transparent", colour = NA)
)
final_plot

ggsave("GO_BubblePlots_6panel_elim.pdf", plot = final_plot, device = "pdf",
       width = 45, height = 55, units = "cm", bg = "transparent")
