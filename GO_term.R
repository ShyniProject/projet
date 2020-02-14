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
                         select(hsapiens_homolog_associated_gene_name)))
  }
}
gene_list_Hs <- unname(unlist(gene_list_Hs))
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mapping_Hs_GO <- getBM(attributes = c("ensembl_gene_id", "name_1006", "namespace_1003"), mart = mart, 
                       filters = c("with_go", "external_gene_name"), values = list(TRUE, gene_list_Hs))
GO_CC <- as.factor((mapping_Hs_GO %>% filter(namespace_1003 == "cellular_component") %>% select(name_1006))$name_1006)
GO_BP <- as.factor((mapping_Hs_GO %>% filter(namespace_1003 == "biological_process") %>% select(name_1006))$name_1006)
GO_MF <- as.factor((mapping_Hs_GO %>% filter(namespace_1003 == "molecular_function") %>% select(name_1006))$name_1006)

keep <- count(GO_CC)[order(count(GO_CC)$freq),][(length(levels(GO_CC))-9):length(levels(GO_CC)),]
ggplot(keep, aes(x="", y=freq, fill=x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set3")
keep <- count(GO_BP)[order(count(GO_BP)$freq),][(length(levels(GO_BP))-9):length(levels(GO_BP)),]
ggplot(keep, aes(x="", y=freq, fill=x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set3")
keep <- count(GO_MF)[order(count(GO_MF)$freq),][(length(levels(GO_MF))-9):length(levels(GO_MF)),]
ggplot(keep, aes(x="", y=freq, fill=x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set3")

#GSEA
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl")
mapping_Dr_Hs <- getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart = mart, 
                       filters = c("with_hsapiens_homolog"), values = TRUE)
gene_list_Hs <- list()
for (gene in dat$id) {
  if(gene %in% mapping_Dr_Hs$ensembl_gene_id) {
    gene_list_Hs[[gene]] <- (mapping_Dr_Hs %>% filter(ensembl_gene_id == gene) %>% 
                               dplyr::select(hsapiens_homolog_associated_gene_name))$hsapiens_homolog_associated_gene_name
  }
}
gsea_dat <- dat[dat$id %in% names(gene_list_Hs),] 
gsea_dat$id <- as.character(gsea_dat$id)
gsea_df <- data.frame(Gene = '', baseMean = 0, log2FoldChange = 0, pvalue = 0, padj = 0, Biotype = '')
for (gene in names(gene_list_Hs)) {
  for (hs_gene in gene_list_Hs[[gene]]) {
    newrow <- unname(dat[dat$id == gene,])[2:7]
    names(newrow) <- names(gsea_df) 
    gsea_df <- rbind(gsea_df, hs_gene = newrow)
  }
}
gsea_df <- gsea_df[2:dim(gsea_df)[1],] %>% mutate(FoldChange = 2^log2FoldChange) %>% dplyr::select(FoldChange)
gsea_df <- as.matrix(gsea_df)
rownames(gsea_df) <- unlist(gene_list_Hs)
gsea_df <- rowmean(gsea_df)
names(geneList) <- rownames(gsea_df)
geneList <- sort(geneList, decreasing = TRUE)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mapping_Hs_GO <- getBM(attributes = c("external_gene_name", "name_1006", "namespace_1003"), mart = mart, 
                       filters = c("with_go", "external_gene_name"), values = list(TRUE, names(geneList)))
mapping_Hs_GO$external_gene_name <- toupper(mapping_Hs_GO$external_gene_name)
GO_CC <- mapping_Hs_GO %>% filter(namespace_1003 == "cellular_component") %>% dplyr::select(name_1006, external_gene_name)
GO_BP <- mapping_Hs_GO %>% filter(namespace_1003 == "biological_process") %>% dplyr::select(name_1006, external_gene_name)
GO_MF <- mapping_Hs_GO %>% filter(namespace_1003 == "molecular_function") %>% dplyr::select(name_1006, external_gene_name)

ES_CC <- GSEA(geneList, TERM2GENE = GO_CC)
ES_MF <- GSEA(geneList, TERM2GENE = GO_MF)
ES_BP <- GSEA(geneList, TERM2GENE = GO_BP)
GSEA.barplot(ES_CC@result, category = 'ID', score = 'pvalue', pvalue = 'NES', 
             sort = 'p.adjust', decreasing = TRUE, top = 20)
GSEA.barplot(ES_MF@result, category = 'ID', score = 'pvalue', pvalue = 'NES', 
             sort = 'p.adjust', decreasing = TRUE, top = 20)
GSEA.barplot(ES_BP@result, category = 'ID', score = 'pvalue', pvalue = 'NES', 
             sort = 'p.adjust', decreasing = TRUE, top = 20)
