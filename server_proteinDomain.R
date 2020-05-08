proteinDomain <- function(input, output, session, org, organismsDbKegg, pvalues, logF, entrezID)
{
  ########################################################
  ##########  Protein Domains  ##########
  ########################################################
  ## Protein domains DATA
  library(biomaRt)
  # pDomains = read.csv("mart_export_INTERPRO.txt")
  
  # dataset = searchDatasets(mart = ensembl, pattern = "dre")
  ensembl = useEnsembl(biomart="ensembl")
  dataset = searchDatasets(mart = ensembl, pattern = organismsDbKegg[org()])
  ensembl = useEnsembl(biomart="ensembl", dataset = dataset$dataset)
  iprData <- getBM(attributes=c('ensembl_gene_id', 'interpro', 'interpro_description'), mart = ensembl)
  
  pD.TERM2GENE = dplyr::select(iprData, interpro, ensembl_gene_id)
  # pD.TERM2GENE$Gene.stable.ID <- mapIds(orgDb, as.vector(pD.TERM2GENE$ensembl_gene_id), 'ENTREZID', 'ENSEMBL')
  pD.TERM2NAME = dplyr::select(iprData, interpro, interpro_description)
  
  gene <- pvalues
  names(gene) <- names(entrezID)
  gene = gene[gene<0.05]
  
  ## ANALYSIS ##
  pDomains.SEA = enricher(names(gene), TERM2GENE = pD.TERM2GENE, TERM2NAME = pD.TERM2NAME)
  
  ## TABLE ## 
  proteinDomains.SEA_description = dplyr::select(pDomains.SEA@result, ID, Description)
  proteinDomains.SEA_description$ID = paste0("<a href=https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001064/", 
                                             proteinDomains.SEA_description$ID," target='_blank'>",proteinDomains.SEA_description$ID,"</a>")
  proteinDomains.SEA_value = dplyr::select(pDomains.SEA@result, everything(), -ID, -Description , -GeneRatio, -geneID,-BgRatio) %>% round(4)
  proteinDomains.SEA = cbind(proteinDomains.SEA_description, proteinDomains.SEA_value)
  proteinDomains.SEAtoSAVE = dplyr::select(proteinDomains.SEA, everything(), -ID)
  
  output$proteinDomains.SEA.Table = DT::renderDataTable({
    proteinDomains.SEA
  }, escape = F) # escape FALSE to make url
  
  output$dl.pDomains <- downloadHandler(
    filename = "proteinDomainsResults_SEA.csv",
    content = function(filename) {
      write.csv(proteinDomains.SEAtoSAVE, filename, row.names = T)
    }
  )
  ## PLOTS ##
  output$dotPlot.pDomains <- renderPlot({ clusterProfiler::dotplot(pDomains.SEA, showCategory=input$categNb_DP.D, ) + ggtitle("dotplot for SEA") })
  
  ########################################################
  ##########  Motifs  ##########
  ########################################################
  organismsMsigdbr = c("org.Hs.eg.db"="Homo sapiens","org.Mm.eg.db"="Mus musculus","org.Rn.eg.db"="Rattus norvegicus",
                       "org.Sc.sgd.db"="Saccharomyces cerevisiae","org.Dm.eg.db"="Drosophila melanogaster",
                       "org.Dr.eg.db"="Danio rerio","org.Bt.eg.db"="Bos taurus","org.Ce.eg.db"="Caenorhabditis elegans",
                       "org.Gg.eg.db"="Caenorhabditis elegans","org.Cf.eg.db"="Canis lupus familiaris","org.Ss.eg.db"="Sus scrofa")
  
  ## MOTIFS (GENE SET) DATA
  motif.geneSet = msigdbr(species = organismsMsigdbr[org()], category = "C3")
  head(motif.geneSet)
  motif.TERM2GENE = motif.geneSet %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
  
  ## GSEA
  ## ANALYSIS ##
  withProgress(message = 'Motif GSEA ... ', value = 0, {
    motifs.GSEA <- GSEA(logF, TERM2GENE=motif.TERM2GENE, pvalueCutoff = 0.15, nPerm = 10000, pAdjustMethod = "BH")
    head(motifs.GSEA@result$ID)
    
    ## TABLE ##
    incProgress(1/2, detail = "Motif dynamic results table...")    
    Motif.GSEA_description = dplyr::select(motifs.GSEA@result, ID, Description)
    Motif.GSEA_description$ID = paste0("<a href=https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=", Motif.GSEA_description$ID," target='_blank'>",Motif.GSEA_description$ID,"</a>")
    Motif.GSEA_value = dplyr::select(motifs.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
    Motif.GSEA = cbind(Motif.GSEA_description, Motif.GSEA_value)
    
    output$Motif.GSEA.Table = DT::renderDataTable({
      Motif.GSEA
    }, escape = F) # escape FALSE to make url
    
    ## PLOTS ##
    incProgress(2/2, detail = "Results plots...")    
    output$dotPlot.Motif <- renderPlot({ dotplot(motifs.GSEA, showCategory=input$categNb_DP.M, ) + ggtitle("dotplot for GSEA") })
    output$ridgePlot.Motif <- renderPlot({ ridgeplot(motifs.GSEA, showCategory=input$categNb_RP.M) })
  })
}