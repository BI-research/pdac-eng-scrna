dir.create("results", showWarnings = FALSE)

library(presto)
library(Seurat)
library(msigdbr)
library(fgsea)
library(dplyr)
library(jsonlite)
library(ggplot2)
library(tibble)
library(edgeR)



##################
### load files ###
##################
## Processed CAFs --------------------------------------------------------------
naive_CAF.obj <- readRDS("data/processed_naive_CAF.rds")
head(naive_CAF.obj);dim(naive_CAF.obj) #[1] 18115  5446

DimPlot(naive_CAF.obj)

naive_CAF.obj$ENG_status <- ifelse(
  FetchData(naive_CAF.obj, vars = "ENG")[,"ENG"] > 0,
  "ENG_pos",
  "ENG_neg"
)

table(naive_CAF.obj$celltype, naive_CAF.obj$ENG_status)
#         ENG_neg ENG_pos
# iCAF       442     291
# myCAF     2214     867
# Others    1337     295



####################
# Prepare gene set #
####################
## Prepare gene set ------------------------------------------------------------

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
C5_NK_cytotoxi_gs <- C5_BP.df %>% dplyr::filter(gs_name == "GOBP_REGULATION_OF_NK_CELL_MEDIATED_CYTOTOXICITY") %>% pull(gene_symbol)
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
  # ECM / Stromal niche
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
X_matrix <- GetAssayData(naive_CAF.obj, layer = "data")   # Seurat v5: use layer=
y <- naive_CAF.obj$ENG_status


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
  scoreType   = "pos",     # ENG_pos enrichment direction
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
  file = "results/Naive_ENG_GSEA_target_genesets.csv",
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

pdf('results/Naive_target_gs_gsea_barPlot.pdf',
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

pdf('results/Naive_target_gs_enrichmentPlot.pdf',
    width = 5, height = 3)
enrich_plots
dev.off()



#########################################
### Donor-level pseudobulk validation ###
#########################################
table(naive_CAF.obj$orig.ident, naive_CAF.obj$ENG_status)


## Build donor x ENG_status pseudobulk counts ----------------------------------
naive_CAF.obj <- JoinLayers(naive_CAF.obj, assay = "RNA")
counts <- GetAssayData(naive_CAF.obj, assay = "RNA", layer = "counts")

meta   <- naive_CAF.obj@meta.data
meta$pb_sample <- paste(meta$orig.ident, meta$ENG_status, sep = "__")
head(meta);dim(meta)

mm <- model.matrix(~ 0 + pb_sample, data = meta)
colnames(mm) <- sub("^pb_sample", "", colnames(mm))
pb_counts <- counts %*% mm
head(pb_counts)

pb_meta <- data.frame(sample = colnames(pb_counts)) %>%
  mutate(
    donor      = sub("__(ENG_neg|ENG_pos)$", "", sample),
    ENG_status = sub("^.*__", "", sample)
  )

# keep only donors with both ENG_pos and ENG_neg pseudobulk samples
paired_donors <- pb_meta %>%
  count(donor, ENG_status) %>%
  count(donor) %>%
  filter(n == 2) %>%
  pull(donor)

length(paired_donors)

pb_meta   <- pb_meta %>% filter(donor %in% paired_donors)
pb_counts <- pb_counts[, pb_meta$sample]
stopifnot(all(colnames(pb_counts) == pb_meta$sample))

table(pb_meta$donor, pb_meta$ENG_status)


## edgeR pseudobulk DE, donor as blocking factor -------------------------------
pb_meta$donor <- factor(pb_meta$donor)
pb_meta$ENG_status <- factor(pb_meta$ENG_status, levels = c("ENG_neg", "ENG_pos"))
dim(pb_counts)

y <- DGEList(counts = pb_counts)
y <- y[filterByExpr(y, group = pb_meta$ENG_status), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ donor + ENG_status, data = pb_meta)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "ENG_statusENG_pos")

edgeR_res <- topTags(qlf, n = Inf, sort.by = "none")$table


## Re-rank genes and re-run fgsea on target_genesets ---------------------------
ranks_pb <- with(edgeR_res, sign(logFC) * sqrt(pmax(F, 0)))
names(ranks_pb) <- rownames(edgeR_res)
ranks_pb <- sort(ranks_pb[is.finite(ranks_pb)], decreasing = TRUE)

# Hallmark gene set (genome-wide check, not just target_genesets)
H_list <- split(H.df$gene_symbol, H.df$gs_name)

set.seed(42)
fgseaRes_pb <- fgsea(
  pathways    = H_list,
  stats       = ranks_pb,
  scoreType   = "std",
  nPermSimple = 1000,
  minSize     = 5,
  maxSize     = 1000
) %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaRes_pb %>% select(pathway, size, NES, pval, padj) %>% print(n = Inf)

NEG_plot_pseudo <- fgseaRes_pb %>%
  filter(NES > 0) %>%
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
    y     = NES + 0.05,
    hjust = 0
  ), size = 5) +
  coord_flip() +
  scale_fill_manual(values = c(
    "ENG_pos enriched" = "#D94F3D"
  )) +
  labs(
    title = "GSEA: CD105+ vs CD105- CAF (donor pseudobulk, Hallmark)",
    x     = NULL,
    y     = "Normalized Enrichment Score (NES)",
    fill  = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

NEG_plot_pseudo

# Re-run fgsea on the same target gene sets used in the cell-level analysis
set.seed(42)
fgseaRes_pb_target <- fgsea(
  pathways    = target_genesets,
  stats       = ranks_pb,
  scoreType   = "std",
  nPermSimple = 1000,
  minSize     = 5,
  maxSize     = 1000
) %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaRes_pb_target %>% select(pathway, size, NES, pval, padj) %>% print(n = Inf)

enrich_plots_pseudo <- lapply(names(target_genesets), function(pw) {
  nes  <- round(fgseaRes_pb_target[fgseaRes_pb_target$pathway == pw, "NES"], 2)
  padj <- round(fgseaRes_pb_target[fgseaRes_pb_target$pathway == pw, "padj"], 4)

  plotEnrichment(
    pathway = target_genesets[[pw]],
    stats   = ranks_pb
  ) +
    labs(
      title    = gsub("_", " ", pw),
      subtitle = paste0("NES = ", nes, "  padj = ", padj)
    ) +
    theme_bw(base_size = 10)
})

pdf('results/Naive_ENG_donor_pseudobulk_HallMark_GSEAbar.pdf',
    width = 6, height = 8)
NEG_plot_pseudo
dev.off()

pdf('results/Naive_ENG_donor_pseudobulk_enrichmentPlot.pdf',
    width = 4, height = 2.5)
enrich_plots_pseudo
dev.off()

write.csv(
  fgseaRes_pb %>% select(pathway, size, NES, pval, padj),
  file = "results/Naive_ENG_donor_pseudobulk_GSEA.csv",
  row.names = FALSE, quote = FALSE
)

