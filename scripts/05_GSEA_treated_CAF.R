# Gene set enrichment analysis: ENG-positive vs ENG-negative CAFs (Treated)
# Gene sets are obtained with msigdbr plus the archived EXTRACELLULAR_MATRIX_PART set.

dir.create("results", showWarnings = FALSE)

library(presto)
library(Seurat)
library(msigdbr)
library(fgsea)
library(dplyr)
library(jsonlite)
library(ggplot2)
library(tibble)





##################
### load files ###
##################
## Processed CAFs --------------------------------------------------------------
treated_CAF.obj <- readRDS("data/processed_treated_CAF.rds")
head(treated_CAF.obj);dim(treated_CAF.obj) #[1] 19008  9967

DimPlot(treated_CAF.obj)

treated_CAF.obj$ENG_status <- ifelse(
  FetchData(treated_CAF.obj, vars = "ENG")[,"ENG"] > 0,
  "ENG_pos",
  "ENG_neg"
)

table(treated_CAF.obj$celltype, treated_CAF.obj$ENG_status)
#         ENG_neg ENG_pos
# iCAF      1926     638
# myCAF     2272     521
# Others    3944     666



####################
# Prepare gene set #
####################
## Prepare gene set ------------------------------------------------------------
msigdbr_collections()

# Hallmark
H.df <- msigdbr(species = "Homo sapiens", collection = "H") 

# C2-KEGG
C2_KEGG.df <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") 

# C2-Reactome
RT.df <- msigdbr(species = "Homo sapiens", collection = "C2",  subcollection = "CP:REACTOME") 

# C5-BP
C5_BP.df <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") 

# C5-CC 
C5_CC.df <- msigdbr(species = "Homo sapiens", collection = "C5",  subcollection = "GO:CC") 


## Target (Interesting) gene sets ----------------------------------------------
# ECM
H_tgf_beta_gs <- H.df %>% dplyr::filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>% pull(gene_symbol)
C5_ecm_assembly_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY") %>% pull(gene_symbol)
C5_ecm_organize_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION") %>% pull(gene_symbol)

# immunity regulation
C5_t_migration_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_T_CELL_MIGRATION") %>% pull(gene_symbol)
C5_t_activation_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_REGULATION_OF_T_CELL_ACTIVATION") %>% pull(gene_symbol)
C5_immresponse_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE") %>% pull(gene_symbol)
RT_immregulate_gs <- RT.df %>% dplyr::filter(gs_name == "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL") %>% pull(gene_symbol)

C5_NK_chemo_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_NATURAL_KILLER_CELL_CHEMOTAXIS") %>% pull(gene_symbol)
KEGG_cytokine_inter_gs <- C2_KEGG.df %>% dplyr::filter(gs_name == "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION") %>% pull(gene_symbol)

C5_NK_active_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_NATURAL_KILLER_CELL_ACTIVATION") %>% pull(gene_symbol)
C5_immune_effect_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS") %>% pull(gene_symbol)

C5_lyp_chemo_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_LYMPHOCYTE_CHEMOTAXIS") %>% pull(gene_symbol)
C5_lyp_migration_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_REGULATION_OF_LYMPHOCYTE_MIGRATION") %>% pull(gene_symbol)


# ECM part (from JSON)
json_url <- "gene_sets/EXTRACELLULAR_MATRIX_PART.v2026.1.Hs.json"
parsed <- fromJSON(json_url)
ecm_gs <- parsed$EXTRACELLULAR_MATRIX_PART$geneSymbols


## Combine all gene sets of interest into a named list ------------------------
# fgsea() requires a named list as the 'pathways' argument
# Passing a plain vector directly will cause an error
target_genesets <- list(
  "HALLMARK_TGF_BETA_SIGNALING" = H_tgf_beta_gs,
  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY" = C5_ecm_assembly_gs,
  "GOBP_POS_REG_ECM_ORGANIZATION" = C5_ecm_organize_gs,
  "EXTRACELLULAR_MATRIX_PART" = ecm_gs,
  
  "GOBP_T_CELL_MIGRATION" = C5_t_migration_gs,
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION" = C5_t_activation_gs,
  "GOBP_POS_REG_IMMUNE_RESPONSE" = C5_immresponse_gs,
  "GOBP_NATURAL_KILLER_CELL_ACTIVATION" = C5_NK_active_gs,
  "GOBP_LYMPHOCYTE_CHEMOTAXIS" = C5_lyp_chemo_gs,
  "GOBP_REGULATION_OF_LYMPHOCYTE_MIGRATION" = C5_lyp_migration_gs
)

sapply(target_genesets, length)  # Check gene count per gene set



#####################
### GSEA analysis ###
#####################
## Build ranked gene list (AUC-based) ------------------------------------------
X_matrix <- GetAssayData(treated_CAF.obj, layer = "data")   # Seurat v5: use layer=
y <- treated_CAF.obj$ENG_status

ENG.genes <- wilcoxauc(
  X = X_matrix,
  y = y,
  groups_use = c("ENG_pos", "ENG_neg")
)
head(ENG.genes);dim(ENG.genes)

# Sort by AUC (descending) from ENG_pos group and convert to named vector
set.seed(42)

ENGpos.genes <- ENG.genes %>%
  filter(group == "ENG_pos") %>%
  mutate(rank_score = sign(logFC) * auc) %>%
  arrange(desc(rank_score)) %>%
  select(feature, rank_score)

ranks <- deframe(ENGpos.genes)
range(ranks)  # -1.0 ~ +1.0

# ENGpos.genes <- ENG.genes %>%
#   filter(group == "ENG_pos") %>%
#   mutate(rank_score = auc - 0.5) %>%  # centered ranking: spans -0.5 to +0.5
#   arrange(desc(rank_score)) %>%
#   select(feature, rank_score)
# 
# ranks <- deframe(ENGpos.genes)

cat(">>> Rank score range (should span negative to positive):", range(ranks), "\n")
cat(">>> Genes with positive rank (ENG_pos enriched):", sum(ranks > 0), "\n")
cat(">>> Genes with negative rank (ENG_neg enriched):", sum(ranks < 0), "\n")


# Add tiny random noise to break ties (does not meaningfully change results)
ranks_jittered <- ranks + runif(length(ranks), min = -1e-7, max = 1e-7)
ranks_jittered


## run GSEA on all target gene sets at once ------------------------------------
set.seed(42)
fgseaRes_target <- fgsea(
  pathways    = target_genesets,
  stats       = ranks_jittered,
  scoreType   = "std",     # pos/neg direction
  nPermSimple = 1000,      # 'nperm' is deprecated; use 'nPermSimple'
  minSize     = 5,         # exclude gene sets smaller than this
  maxSize     = 1000        # exclude gene sets larger than this
)


## Summarize results -----------------------------------------------------------
fgseaResTidy <- fgseaRes_target %>%
  as_tibble() %>%
  arrange(desc(NES))

cat("\n>>> GSEA results (sorted by NES descending):\n")
fgseaResTidy %>%
  select(pathway, size, NES, pval, padj) %>%
  print(n = Inf)

# Save results 
write.csv(
  fgseaResTidy %>% select(pathway, size, NES, pval, padj),
  file = "results/Treated_ENG_GSEA_target_genesets.csv",
  row.names = FALSE, quote = F
)



############
### Plot ###
############
## NES based bar plot ----------------------------------------------------------
NEG_plot <- fgseaResTidy %>%
  mutate(
    pathway   = gsub("_", " ", pathway),
    direction = ifelse(NES > 0, "ENG_pos enriched", "ENG_neg enriched"),
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  ) %>%
  ggplot(aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(
    label = sig_label,
    y     = ifelse(NES > 0, NES + 0.05, NES - 0.05),
    hjust = ifelse(NES > 0, 0, 1)
  ), size = 5) +
  coord_flip() +
  scale_fill_manual(values = c(
    "ENG_pos enriched" = "#D94F3D",
    "ENG_neg enriched" = "#4F87C5"
  )) +
  labs(
    title = "GSEA: CD105+ vs CD105- CAF",
    x     = NULL,
    y     = "Normalized Enrichment Score (NES)",
    fill  = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

NEG_plot

pdf('results/Treated_target_gs_gsea_barPlot.pdf',
    width = 7.5, height = 4)
NEG_plot
dev.off()


## enrichment plot -------------------------------------------------------------
# enrichment plot 
enrich_plots <- lapply(names(target_genesets), function(pw) {
  nes  <- round(fgseaRes_target[fgseaRes_target$pathway == pw, "NES"], 2)
  padj <- round(fgseaRes_target[fgseaRes_target$pathway == pw, "padj"], 4)
  
  plotEnrichment(
    pathway = target_genesets[[pw]],
    stats   = ranks_jittered
  ) +
    labs(
      title    = gsub("_", " ", pw),
      subtitle = paste0("NES = ", nes, "  padj = ", padj)
    ) +
    theme_bw(base_size = 10)
})

enrich_plots

pdf('results/Treated_target_gs_enrichmentPlot.pdf',
    width = 5, height = 3)
enrich_plots
dev.off()

