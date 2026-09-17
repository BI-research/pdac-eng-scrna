dir.create("results", showWarnings = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(1234)




##################
### load files ###
##################
## Processed CAFs --------------------------------------------------------------
naive_CAF.obj <- readRDS("data/processed_naive_CAF.rds")
head(naive_CAF.obj);dim(naive_CAF.obj) #[1] 18115  5446

#DimPlot(naive_CAF.obj)

naive_CAF.obj$ENG_status <- ifelse(
  FetchData(naive_CAF.obj, vars = "ENG")[,"ENG"] > 0,
  "ENG_pos",
  "ENG_neg"
)

table( naive_CAF.obj$celltype, naive_CAF.obj$ENG_status)
#         ENG_neg ENG_pos
# iCAF       442     291
# myCAF     2214     867
# Others    1337     295



###########
### DEG ###
###########
Idents(naive_CAF.obj) <- "ENG_status"


## DEG analysis ----------------------------------------------------------------
naive_CAF.obj <- PrepSCTFindMarkers(naive_CAF.obj)

PosENG_UP_AllCAF_naive <- FindMarkers(naive_CAF.obj, 
                                      ident.1 = "ENG_pos", 
                                      ident.2 = "ENG_neg", 
                                      only.pos = TRUE)
NegENG_UP_AllCAF_naive <- FindMarkers(naive_CAF.obj, 
                                      ident.1 = "ENG_neg", 
                                      ident.2 = "ENG_pos", 
                                      only.pos = TRUE)

dim(PosENG_UP_AllCAF_naive) # 2936 genes
dim(NegENG_UP_AllCAF_naive) # 4948 genes


# filtering with log2FC & q-value
PosENG_UP_AllCAF_naive_list <- PosENG_UP_AllCAF_naive[PosENG_UP_AllCAF_naive$avg_log2FC >= 0.25 & PosENG_UP_AllCAF_naive$p_val_adj < 0.01,]
PosENG_UP_AllCAF_naive_list$UP_in <- "ENG_pos"
head(PosENG_UP_AllCAF_naive_list);dim(PosENG_UP_AllCAF_naive_list) # 352 genes

NegENG_UP_AllCAF_naive_list <- NegENG_UP_AllCAF_naive[NegENG_UP_AllCAF_naive$avg_log2FC >= 0.25 & NegENG_UP_AllCAF_naive$p_val_adj < 0.01,]
NegENG_UP_AllCAF_naive_list$UP_in <- "ENG_neg"
head(NegENG_UP_AllCAF_naive_list);dim(NegENG_UP_AllCAF_naive_list) # 224 genes


# FN1, COL1A1, COL1A2, COL2A1, COL3A1, ACTA2 (a-SMA), VIM (Vinentin), TGFB1
target_genes <- c('ENG','FN1','COL1A1','COL1A2','COL2A1','COL3A1','ACTA2','VIM','TGFB1')
ECM_genes <- c("AGRN","AMTN","ANG","APLP1","CAV3","COL10A1","COL11A1","COL13A1",
               "COL15A1","COL16A1","COL18A1","COL1A2","COL3A1","COL4A2","COL4A3",
               "COL4A4","COL4A5","COL5A1","COL5A2","COL5A3","COL6A3","COL7A1",
               "COL8A1","COL9A1","COL9A2","COL9A3","COLQ","DMD","DST","EFEMP2",
               "ERBIN","FBN1","LAMA2","LAMA3","LAMA4","LAMB1","LAMB2","LAMC1",
               "LUM","MAGEE1","MUC5AC","NID2","ODAM","SGCA","SGCB","SGCD","SGCE",
               "SGCG","SMC3","SNTB1","SNTB2","SNTG1","SNTG2","SSPN","TINAG","TNXB","USH2A")
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/EXTRACELLULAR_MATRIX_PART.html


unique(c(target_genes, ECM_genes)) %>% length() # 64 ECM-related genes 

ECM_genes_PosENG_UP_AllCAF_naive <- PosENG_UP_AllCAF_naive[rownames(PosENG_UP_AllCAF_naive) %in% unique(c(target_genes, ECM_genes)), ]
ECM_genes_NegENG_UP_AllCAF_naive <- NegENG_UP_AllCAF_naive[rownames(NegENG_UP_AllCAF_naive) %in% unique(c(target_genes, ECM_genes)), ]

ECM_genes_PosENG_UP_AllCAF_naive; dim(ECM_genes_PosENG_UP_AllCAF_naive) # [1] 35  5
ECM_genes_NegENG_UP_AllCAF_naive; dim(ECM_genes_NegENG_UP_AllCAF_naive) # [1] 10  5

# In significant DEG list
ECM_genes_PosENG_UP_AllCAF_naive_list <- PosENG_UP_AllCAF_naive_list[rownames(PosENG_UP_AllCAF_naive_list) %in% unique(c(target_genes, ECM_genes)), ]
ECM_genes_NegENG_UP_AllCAF_naive_list <- NegENG_UP_AllCAF_naive_list[rownames(NegENG_UP_AllCAF_naive_list) %in% unique(c(target_genes, ECM_genes)), ]

ECM_genes_PosENG_UP_AllCAF_naive_list; dim(ECM_genes_PosENG_UP_AllCAF_naive_list) # [1] 20  5
ECM_genes_NegENG_UP_AllCAF_naive_list; dim(ECM_genes_NegENG_UP_AllCAF_naive_list) # [1] 2  5


## Visualization with VolcanoPlot ----------------------------------------------

# FC negative value form
temp_NegENG_UP_AllCAF_naive <- NegENG_UP_AllCAF_naive
temp_NegENG_UP_AllCAF_naive$avg_log2FC <- -temp_NegENG_UP_AllCAF_naive$avg_log2FC

# merge ENG positive/negative UP All DEGs
Naive_AllCAF_DEGs.df <- rbind(PosENG_UP_AllCAF_naive, temp_NegENG_UP_AllCAF_naive)
head(Naive_AllCAF_DEGs.df)

# DEGs classification
res_Naive_AllCAF_DEGs.df <- na.omit(Naive_AllCAF_DEGs.df %>%
                                      mutate(Expression = case_when(avg_log2FC >= 0.25 & p_val_adj < 0.01 ~ 'ENG_pos',
                                                                    avg_log2FC <= -0.25 & p_val_adj < 0.01 ~ 'ENG_neg',
                                                                    TRUE ~ "Unchanged")))

res_Naive_AllCAF_DEGs.df
res_Naive_AllCAF_DEGs.df %>% group_by(Expression) %>% summarise(n=n())
 # Expression     n
 # ENG_neg      224
 # ENG_pos      352
 # Unchanged   7308

# save the table
write.table(res_Naive_AllCAF_DEGs.df,
            "results/Naive_AllCAF_ENG_pos_vs_neg_DEGs.tsv",
            sep =  "\t", col.names = T, row.names = T, quote = F)


# highlight ECM markers
# In All DEGs
HighLight_Naive_AllCAF <- res_Naive_AllCAF_DEGs.df[c(rownames(ECM_genes_PosENG_UP_AllCAF_naive),
                                                     rownames(ECM_genes_NegENG_UP_AllCAF_naive)),]
HighLight_Naive_AllCAF$GeneSymbol <- rownames(HighLight_Naive_AllCAF)

# In significant DEGs 
HighLight_Naive_AllCAF_sig <- res_Naive_AllCAF_DEGs.df[c(rownames(ECM_genes_PosENG_UP_AllCAF_naive_list),
                                                     rownames(ECM_genes_NegENG_UP_AllCAF_naive_list)),]
HighLight_Naive_AllCAF_sig$GeneSymbol <- rownames(HighLight_Naive_AllCAF_sig)


# volcano plot
pdf('results/Naive_ENG_DEGs_volcanoPlot.pdf',
    width = 7.5, height = 7)
ggplot(res_Naive_AllCAF_DEGs.df, 
       aes(avg_log2FC, -log(p_val_adj,10))) + # -log10 conversion  
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p.adj")) + 
  ggtitle('Naive_AllCAF') +
  scale_color_manual(values = c("darkblue","red","gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype='dotted') + 
  geom_hline(yintercept = 2, linetype='dotted') + 
  theme_bw() +
  geom_text_repel(data = HighLight_Naive_AllCAF, 
                  aes(label = GeneSymbol),
                  size = 3, 
                  in.segment.length = 0, 
                  segment.size = 0.6, 
                  seed = 42, 
                  box.padding = 1,
                  max.overlaps = Inf) +
  xlim(c(-6,6)) + ylim(c(0,65))


ggplot(res_Naive_AllCAF_DEGs.df, 
       aes(avg_log2FC, -log(p_val_adj,10))) + # -log10 conversion  
  geom_point(aes(color = Expression), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p.adj")) + 
  ggtitle('Naive_AllCAF (Highlighted significant genes)') +
  scale_color_manual(values = c("darkblue","red","gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype='dotted') + 
  geom_hline(yintercept = 2, linetype='dotted') + 
  theme_bw() +
  geom_text_repel(data = HighLight_Naive_AllCAF_sig, 
                  aes(label = GeneSymbol),
                  size = 3, 
                  in.segment.length = 0, 
                  segment.size = 0.6, 
                  seed = 42, 
                  box.padding = 1,
                  max.overlaps = Inf) +
  xlim(c(-6,6)) + ylim(c(0,65))
dev.off()


