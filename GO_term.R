library(dplyr)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(PPInfer)

ORA <- function(mapping, DE_genes) {
  ORA_df <- data.frame(term = unique(mapping$term))
  size <- c()
  DE_n <- c()
  pval <- c()
  for (term in ORA_df$term) {
    k <- sum(mapping$term == term)
    x <- sum(mapping$term == term & mapping$gene_name %in% DE_genes)
    size <- c(size, k)
    DE_n <- c(DE_n, x)
    pval <- c(pval, dhyper(x, length(DE_genes), length(unique(mapping$gene_name)), k))
  }
  ORA_df$size <- size
  ORA_df$DE_n <- DE_n
  ORA_df$pval <- pval
  ORA_df <- ORA_df %>% arrange(pval)
  return(ORA_df)
}

dat <- read.table("GSE129081_small.csv", header = TRUE, sep = ",")
dat <- dat %>% filter(!is.na(padj))

# OPTIONAL: to update GO data
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mapping_Hs_GO <- getBM(attributes = c("external_gene_name", "name_1006", "namespace_1003"), mart = mart, 
                       filters = c("with_go"), values = TRUE)
mapping_Hs_GO <- mapping_Hs_GO %>% filter(namespace_1003 != "" & namespace_1003 != "go")
write.table(mapping_Hs_GO, "HsGenesToGO.txt")
mapping_Ensembl_Entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)
write.table(mapping_Ensembl_Entrez, "EnsemblToEntrez_Drerio.txt")

#SEA analysis
DE <- dat %>% dplyr::select(Orthologous_human_gene, padj) %>% filter(padj < 0.05)
all <- dat %>% dplyr::select(Orthologous_human_gene, padj)
DE_genes <- as.character(DE$Orthologous_human_gene)
DE_genes <- DE_genes[!is.na(DE_genes)]
all_genes <- as.character(all$Orthologous_human_gene)
all_genes <- all_genes[!is.na(all_genes)]
mapping_Hs_GO <- read.table("HsGenesToGO.txt")
mapping_useful <- mapping_Hs_GO %>% filter(external_gene_name %in% all_genes)

ORA_dfs <- list()
cats <- c("biological_process", "cellular_component", "molecular_function")
for (cat in cats) {
  da <- mapping_useful %>% filter(namespace_1003 == cat)
  da <- data.frame(gene_name = da$external_gene_name, term = da$name_1006)
  
  ORA_dfs[[cat]] <- ORA(da, DE_genes)
}
# BP
ggplot(ORA_dfs[[1]][1:10,], aes(x = term, y = DE_n/size, fill = pval)) + geom_col()
# CC
ggplot(ORA_dfs[[2]][1:10,], aes(x = term, y = DE_n/size, fill = pval)) + geom_col()
# MF
ggplot(ORA_dfs[[3]][1:10,], aes(x = term, y = DE_n/size, fill = pval)) + geom_col()

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
