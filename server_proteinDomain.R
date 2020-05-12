ORA <- function(mapping, DE_genes) 
{
  ORA_df <- data.frame(interpro = unique(mapping$interpro))
  size <- c()
  DE_n <- c()
  pval <- c()
  for (interpro in ORA_df$interpro) {
    k <- sum(mapping$interpro == interpro)
    x <- sum(mapping$interpro == interpro & mapping$entrezgene_id %in% DE_genes)
    size <- c(size, k)
    DE_n <- c(DE_n, x)
    pval <- c(pval, dhyper(x, length(DE_genes), length(unique(mapping$entrezgene_id)), k))
  }
  ORA_df$size <- size
  ORA_df$DE_n <- DE_n
  ORA_df$pval <- pval
  ORA_df <- ORA_df %>% arrange(pval)
  return(ORA_df)
}

proteinDomain <- function(input, output, session, org, organismsDbKegg, pvalues, logF, entrezID, minGS, maxGS, nPerm, pvAdjustMethod)
{
  ########################################################
  ##########  Protein Domains  ##########
  ########################################################
  ## Protein domains DATA
  library(biomaRt)
  # dataset = searchDatasets(mart = ensembl, pattern = "dre")
  ensembl = useEnsembl(biomart="ensembl")
  dataset = searchDatasets(mart = ensembl, pattern = organismsDbKegg[org])
  ensembl = useEnsembl(biomart="ensembl", dataset = dataset$dataset)

  iprData <- getBM(attributes=c('entrezgene_id', 'interpro', 'interpro_description'), mart = ensembl,
                   filter="entrezgene_id",
                   values = entrezID,
                   uniqueRows = TRUE)
  
  pD.TERM2GENE = dplyr::select(iprData, interpro, entrezgene_id)
  # pD.TERM2GENE$Gene.stable.ID <- mapIds(orgDb, as.vector(pD.TERM2GENE$ensembl_gene_id), 'ENTREZID', 'ENSEMBL')
  pD.TERM2NAME = dplyr::select(iprData, interpro, interpro_description)
  
  gene <- pvalues[!is.na(names(pvalues))]
  gene <- gene[gene<input$pv.DomainsSEA]
  
  # names(gene) <- entrezID
  observe({gene_vec()})
  gene_vec <- eventReactive(input$pv.DomainsSEA,{
    gene <- gene[gene<input$pv.DomainsSEA]
    gene
  })
    
    
    ## ANALYSIS ##
    # pD.TERM2GENE <- pD.TERM2GENE[pD.TERM2GENE$interpro != "",]
    pd = (pD.TERM2GENE[pD.TERM2GENE$interpro != "",])
    print("D2BUT SEA")
    pDomains.SEA <- ORA(mapping = pd, DE_genes = names(gene))
  
  
  # pDomains.SEA = enricher(names(gene), 
  #                         TERM2GENE = pD.TERM2GENE,
  #                         TERM2NAME = pD.TERM2NAME,
  #                         minGSSize = minGS,
  #                         maxGSSize = maxGS,
  #                         pAdjustMethod = pvAdjustMethod)
  
  ## TABLE ## 
  # proteinDomains.SEA_description = dplyr::select(pDomains.SEA@result, ID, Description)
  # proteinDomains.SEA_description$ID = paste0("<a href=https://www.ebi.ac.uk/interpro/entry/InterPro/", 
  #                                            proteinDomains.SEA_description$ID," target='_blank'>",proteinDomains.SEA_description$ID,"</a>")
  # proteinDomains.SEA_value = dplyr::select(pDomains.SEA@result, everything(), -ID, -Description , -GeneRatio, -geneID,-BgRatio) %>% round(5)
  # proteinDomains.SEA = cbind(proteinDomains.SEA_description, proteinDomains.SEA_value)
  # proteinDomains.SEAtoSAVE = dplyr::select(proteinDomains.SEA, everything(), -ID)
  pDomains.SEAtoPlot <- pDomains.SEA
  pDomains.SEAtoPlot$interpro = paste0("<a href=https://www.ebi.ac.uk/interpro/entry/InterPro/",
                                       pDomains.SEAtoPlot$interpro," target='_blank'>", pDomains.SEAtoPlot$interpro,"</a>")
  # print(head(pDomains.SEA, n = 5))
  output$proteinDomains.SEA.Table = DT::renderDataTable({
    # proteinDomains.SEA
    pDomains.SEAtoPlot
  }, escape = F) # escape FALSE to make url
  
  output$dl.pDomains <- downloadHandler(
    filename = "proteinDomainsResults_SEA.csv",
    content = function(filename) {
      write.csv(pDomains.SEA, filename, row.names = T)
    }
  )
  ## PLOTS ##
  # DP.pDomains <- clusterProfiler::dotplot(pDomains.SEA, showCategory=input$categNb_DP.D, ) + ggtitle("dotplot for SEA")
  pDomains.SEA
  DP.pDomains <- ggplot(dplyr::arrange(pDomains.SEA, pval)[1:12,], aes(x = interpro, y = DE_n/size, fill = pval)) + geom_col() + ggtitle("Protein domains dotplot for SEA")
  # DP.pDomains <- ggplot(dplyr::arrange(pDomains.SEA, pval)[1:input$categNb_DP.D,], aes(x = interpro, y = DE_n/size, fill = pval)) + geom_col() + ggtitle("Protein domains dotplot for SEA")
  output$dotPlot.pDomains <- renderPlot({ DP.pDomains })
  
  ## DOWNLOAD PLOT
  output$dl.SEAproteinDomainsDotPlot <- downloadHandler(
    filename = "proteinDomainsDotPlot_SEA.pdf",
    content = function(file) 
    {
      ggsave(file, plot=DP.pDomains)
    }
  )
  
  
  ########################################################
  ##########  Motifs  ##########
  ########################################################
  organismsMsigdbr = c("org.Hs.eg.db"="Homo sapiens","org.Mm.eg.db"="Mus musculus","org.Rn.eg.db"="Rattus norvegicus",
                       "org.Sc.sgd.db"="Saccharomyces cerevisiae","org.Dm.eg.db"="Drosophila melanogaster",
                       "org.Dr.eg.db"="Danio rerio","org.Bt.eg.db"="Bos taurus","org.Ce.eg.db"="Caenorhabditis elegans",
                       "org.Gg.eg.db"="Caenorhabditis elegans","org.Cf.eg.db"="Canis lupus familiaris","org.Ss.eg.db"="Sus scrofa")
  
  ## MOTIFS (GENE SET) DATA
  motif.geneSet = msigdbr(species = organismsMsigdbr[org], category = "C3")
  head(motif.geneSet)
  motif.TERM2GENE = motif.geneSet %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
  
  ## GSEA
  ## ANALYSIS ##
  withProgress(message = 'Motif GSEA ... ', value = 0, {
    motifs.GSEA <- GSEA(logF, 
                        TERM2GENE=motif.TERM2GENE, 
                        pvalueCutoff = 0.15, 
                        nPerm = nPerm, 
                        pAdjustMethod = pvAdjustMethod, 
                        minGSSize = minGS,
                        maxGSSize = maxGS)
    head(motifs.GSEA@result$ID)
    
    ## TABLE ##
    incProgress(1/2, detail = "Motif dynamic results table...")    
    Motif.GSEA_description = dplyr::select(motifs.GSEA@result, ID, Description)
    Motif.GSEA_description$ID = paste0("<a href=https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=", Motif.GSEA_description$ID," target='_blank'>",Motif.GSEA_description$ID,"</a>")
    Motif.GSEA_value = dplyr::select(motifs.GSEA@result, setSize:p.adjust, -NES) %>% round(5)
    Motif.GSEA = cbind(Motif.GSEA_description, Motif.GSEA_value)
    Motif.GSEAtoSAVE = dplyr::select(Motif.GSEA, everything(), -Description, -setSize)
    
    output$Motif.GSEA.Table = DT::renderDataTable({
      Motif.GSEAtoSAVE
    }, escape = F) # escape FALSE to make url
    
    ## PLOTS ##
    incProgress(2/2, detail = "Results plots...")   
    Motif.DP <- dotplot(motifs.GSEA, showCategory=input$categNb_DP.M ) + ggtitle("Motif dotplot for GSEA") 
    output$dotPlot.Motif <- renderPlot({ dotplot(motifs.GSEA, showCategory=input$categNb_DP.M ) + ggtitle("Motif dotplot for GSEA")  })
    Motif.RP <- ridgeplot(motifs.GSEA, showCategory=input$categNb_RP.M) + ggtitle("Motif Ridge Plot for GSEA") 
    output$ridgePlot.Motif <- renderPlot({ ridgeplot(motifs.GSEA, showCategory=input$categNb_RP.M) + ggtitle("Motif Ridge Plot for GSEA")  })
    
    ## DOWNLOAD DOTPLOT
    output$dl.GSEAmotifDotPlot <- downloadHandler(
      filename = "motifDotPlot_GSEA.pdf",
      content = function(file) 
      {
        ggsave(file, plot=Motif.DP)
      }
    )
    
    output$dl.SEAmotifRidgePlot <- downloadHandler(
      filename = "motifRidgePlot_GSEA.pdf",
      content = function(file) 
      {
        ggsave(file, plot=Motif.RP)
      }
    )
  })
}