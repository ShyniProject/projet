library(dplyr)
library(tidyr)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(PPInfer)

dat <- read.table("GSE129081_small.csv", header = TRUE, sep = ",")
dat <- dat %>% filter(!is.na(padj))

# OPTIONAL: to update GO data
mapping_Ensembl_Entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)
write.table(mapping_Ensembl_Entrez, "EnsemblToEntrez_Drerio.txt")

#SEA analysis
DE <- dat %>% dplyr::select(id, padj) %>% filter(padj < 0.05)
DE_genes <- as.character(DE$id)
mapping_Ensembl_Entrez <- read.table("EnsemblToEntrez_Drerio.txt")
DE_genes <- (mapping_Ensembl_Entrez %>% filter(ensembl_gene_id %in% DE_genes))$entrezgene_id
DE_genes <- DE_genes[!is.na(DE_genes)]

ORA_dfs <- list()
cats <- c("BP", "CC", "MF")
for (cat in cats) {
  ORA_df <- enrichGO(DE_genes, 'org.Dr.eg.db', ont=cat, pvalueCutoff=0.05)@result
  ORA_df <- data.frame(GO_term = ORA_df$ID, desc = ORA_df$Description, ratio = ORA_df$GeneRatio, pval = ORA_df$qvalue)
  ORA_df <- ORA_df %>% separate(ratio, c("DE_n", "size"))
  ORA_df <- ORA_df %>% mutate(DE_n = as.numeric(DE_n), size = as.numeric(size))
  
  ORA_dfs[[cat]] <- ORA_df
}
# BP
ggplot(ORA_dfs[[1]][1:10,], aes(x = desc, y = DE_n/size, fill = pval)) + geom_col()
# CC
ggplot(ORA_dfs[[2]][1:10,], aes(x = desc, y = DE_n/size, fill = pval)) + geom_col()
# MF
ggplot(ORA_dfs[[3]][1:10,], aes(x = desc, y = DE_n/size, fill = pval)) + geom_col()

#GSEA
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
mapping_Ensembl_Entrez <- read.table("EnsemblToEntrez_Drerio.txt")
mapping_Ensembl_Entrez <- mapping_Ensembl_Entrez %>% filter(ensembl_gene_id %in% dat$id & 
                                                              !is.na(entrezgene_id) & 
                                                              !duplicated(ensembl_gene_id) &
                                                              !duplicated(entrezgene_id))
gsea_dat <- dat[dat$id %in% mapping_Ensembl_Entrez$ensembl_gene_id,] 
gsea_df <- data.frame(gsea_dat)
gsea_df <- gsea_df %>% mutate(FoldChange = 2^log2FoldChange, GeneID = mapping_Ensembl_Entrez$entrezgene_id) %>% 
  dplyr::select(FoldChange)
geneList <- gsea_df$FoldChange
names(geneList) <- mapping_Ensembl_Entrez$entrezgene_id
geneList <- sort(geneList, decreasing = TRUE)

ES_CC <- gseGO(geneList, OrgDb = "org.Dr.eg.db", pvalueCutoff = 1, ont = "CC")
ES_MF <- gseGO(geneList, OrgDb = "org.Dr.eg.db", pvalueCutoff = 1, ont = "MF")
ES_BP <- gseGO(geneList, OrgDb = "org.Dr.eg.db", pvalueCutoff = 1, ont = "BP")
GSEA.barplot(ES_CC@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
             sort = 'p.adjust', decreasing = TRUE, top = 20, 
             title = "Main differentially expressed Cellular Components (GO)")
GSEA.barplot(ES_MF@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
             sort = 'p.adjust', decreasing = TRUE, top = 20,
             title = "Main differentially expressed Molecular Functions (GO)")
GSEA.barplot(ES_BP@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
             sort = 'p.adjust', decreasing = TRUE, top = 20,
             title = "Main differentially expressed Biological Processes (GO)")
