library(plyr)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(EnrichmentBrowser)
library(analytics)
library(clusterProfiler)
library(PPInfer)

dat <- read.table("GSE129081_small.csv", header = TRUE, sep = ",")
dat <- dat %>% filter(!is.na(padj))

#SEA
sea_dat <- dat %>% dplyr::select(id, padj) %>% filter(padj < 0.05)
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
mapping_Dr_Hs <- getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart = mart, 
                       filters = c("with_hsapiens_homolog"), values = TRUE)
gene_list <- as.character(sea_dat$id)
gene_list_Hs <- c()
for (gene in gene_list) {
  if(gene %in% mapping_Dr_Hs$ensembl_gene_id) {
    gene_list_Hs <- c(gene_list_Hs,
                      (mapping_Dr_Hs %>% filter(ensembl_gene_id == gene) %>% 
                         dplyr::select(hsapiens_homolog_associated_gene_name)))
  }
}
gene_list_Hs <- unname(unlist(gene_list_Hs))
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mapping_Hs_GO <- getBM(attributes = c("ensembl_gene_id", "name_1006", "namespace_1003"), mart = mart, 
                       filters = c("with_go", "external_gene_name"), values = list(TRUE, gene_list_Hs))
GO_CC <- as.factor((mapping_Hs_GO %>% filter(namespace_1003 == "cellular_component") %>% dplyr::select(name_1006))$name_1006)
GO_BP <- as.factor((mapping_Hs_GO %>% filter(namespace_1003 == "biological_process") %>% dplyr::select(name_1006))$name_1006)
GO_MF <- as.factor((mapping_Hs_GO %>% filter(namespace_1003 == "molecular_function") %>% dplyr::select(name_1006))$name_1006)

keep <- plyr::count(GO_CC)[order(plyr::count(GO_CC)$freq),][(length(levels(GO_CC))-11):length(levels(GO_CC)),]
ggplot(keep, aes(x="", y=freq, fill=x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set3")
keep <- plyr::count(GO_BP)[order(plyr::count(GO_BP)$freq),][(length(levels(GO_BP))-11):length(levels(GO_BP)),]
ggplot(keep, aes(x="", y=freq, fill=x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set3")
keep <- plyr::count(GO_MF)[order(plyr::count(GO_MF)$freq),][(length(levels(GO_MF))-11):length(levels(GO_MF)),]
ggplot(keep, aes(x="", y=freq, fill=x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set3")

#GSEA
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
mapping_Dr_Hs <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart, 
                       filters = c("ensembl_gene_id"), values = dat$id)
mapping_Dr_Hs <- mapping_Dr_Hs[!is.na(mapping_Dr_Hs$entrezgene_id),]
mapping_Dr_Hs <- mapping_Dr_Hs[!duplicated(mapping_Dr_Hs$ensembl_gene_id),]
mapping_Dr_Hs <- mapping_Dr_Hs[!duplicated(mapping_Dr_Hs$entrezgene_id),]
gsea_dat <- dat[dat$id %in% mapping_Dr_Hs$ensembl_gene_id,] 
gsea_df <- data.frame(gsea_dat)
gsea_df <- gsea_df %>% mutate(FoldChange = 2^log2FoldChange, GeneID = mapping_Dr_Hs$entrezgene_id) %>% dplyr::select(FoldChange)
geneList <- gsea_df$FoldChange
names(geneList) <- mapping_Dr_Hs$entrezgene_id
geneList <- sort(geneList, decreasing = TRUE)

ES_CC <- gseGO(geneList, OrgDb = "org.Dr.eg.db", pvalueCutoff = 0.5, ont = "CC")
ES_MF <- gseGO(geneList, OrgDb = "org.Dr.eg.db", pvalueCutoff = 0.5, ont = "MF")
ES_BP <- gseGO(geneList, OrgDb = "org.Dr.eg.db", pvalueCutoff = 0.5, ont = "BP")
GSEA.barplot(ES_CC@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
             sort = 'p.adjust', decreasing = TRUE, top = 20)
GSEA.barplot(ES_MF@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
             sort = 'p.adjust', decreasing = TRUE, top = 20)
GSEA.barplot(ES_BP@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
             sort = 'p.adjust', decreasing = TRUE, top = 20)