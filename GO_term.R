# OPTIONAL: to update GO data
Update_Go_data <- function(){
  mapping_Ensembl_Entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)
  write.table(mapping_Ensembl_Entrez, "EnsemblToEntrez_Drerio.txt")
}

#SEA analysis
SEAanalysis <- function(pvals, pvalcut, organismDb){
  DE_pvals <- pvals[pvals < pvalcut]
  DE_genes <- names(DE_pvals)
  DE_genes <- DE_genes[!is.na(DE_genes)]
  
  ORA_dfs <- list()
  cats <- c("BP", "CC", "MF")
  for (cat in cats) {
    ORA_df <- enrichGO(DE_genes, organismDb, ont=cat, pvalueCutoff=pvalcut)@result
    ORA_df <- data.frame(GO_term = ORA_df$ID, desc = ORA_df$Description, ratio = ORA_df$GeneRatio, pval = ORA_df$qvalue)
    ORA_df <- ORA_df %>% tidyr::separate(ratio, c("DE_n", "size"))
    ORA_df <- ORA_df %>% mutate(DE_n = as.numeric(DE_n), size = as.numeric(size))
    
    ORA_dfs[[cat]] <- ORA_df
  }
  return(ORA_dfs)
}

#Plots
bp_SEA_Plot=function(ORA_dfs){
  ggplot(ORA_dfs[[1]][1:10,], aes(x = desc, y = DE_n/size, fill = pval)) + geom_col()
}

mf_SEA_Plot=function(ORA_dfs){
  ggplot(ORA_dfs[[3]][1:10,], aes(x = desc, y = DE_n/size, fill = pval)) + geom_col()
}

cc_SEA_Plot=function(ORA_dfs){
  ggplot(ORA_dfs[[2]][1:10,], aes(x = desc, y = DE_n/size, fill = pval)) + geom_col()
}




###########################################################################################################################
#GSEA

GSEAanalysis <- function(geneList, organismDb){
  ES_CC <- gseGO(geneList, OrgDb = organismDb, pvalueCutoff = 1, ont = "CC")
  ES_MF <- gseGO(geneList, OrgDb = organismDb, pvalueCutoff = 1, ont = "MF")
  ES_BP <- gseGO(geneList, OrgDb = organismDb, pvalueCutoff = 1, ont = "BP")
  return(list("CC" =ES_CC,"MF"= ES_MF,"BP" = ES_BP))
}

bp_GSEA_Plot <- function(data){
  bpdata = data$BP
  GSEA.barplot(bpdata@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
               sort = 'p.adjust', decreasing = TRUE, top = 20,
               title = "Main differentially expressed Biological Processes (GO)")
}

mf_GSEA_Plot=function(data){
  mfdata = data$MF
  GSEA.barplot(mfdata@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
               sort = 'p.adjust', decreasing = TRUE, top = 20,
               title = "Main differentially expressed Molecular Functions (GO)")
}

cc_GSEA_Plot=function(data){
  ccdata = data$CC
  GSEA.barplot(ccdata@result, category = 'Description', score = 'NES', pvalue = 'pvalue', 
               sort = 'p.adjust', decreasing = TRUE, top = 20, 
               title = "Main differentially expressed Cellular Components (GO)")
}




