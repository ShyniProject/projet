function(input, output, session) {
  
  observe({
    analysis()
  })
  
  analysis <- eventReactive(input$data,{
      
      ########################################################
      ##########  data preprocesses  ##########
      ########################################################
      ## RNAseq data ##
      d <- read.csv(input$data$datapath)
      entrezID <- mapIds(org.Dr.eg.db, as.vector(d$id), 'ENTREZID', 'ENSEMBL')
      
      ## for SEA : p-values vector and gene ids
      pvalues <- d$padj 
      names(pvalues) <- entrezID
      
      ## for GSEA : LogF vector and gene ids
      logF = d$log2FoldChange 
      names(logF) <- entrezID
      logF = logF[which(!is.na(names(logF)))] ; logF = na.exclude(logF)#17 809 samples
      ## delete duplicated genes -> median value computed
      if(length(which(duplicated(names(logF))==T))>0)
      {
        df = data.frame("id" = names(logF), "value" = logF)
        df %>%
          group_by(id) %>%
          summarise_all(median) %>%
          data.frame() -> newdf
        logF = newdf$value
        names(logF) = newdf$id # 17 777
      }
      
      logF = sort(logF, decreasing = TRUE) # decreasing vector to gsea
    
      
      ########################################################
      ##########  KEGG  ##########
      ########################################################
      
      ### SEA ###
        ## ANALYSIS ##
      withProgress(message = 'SEA ... ', value = 0, {

        gene = pvalues[pvalues<0.05]
        Kegg.SEA = enrichKEGG(gene = names(gene),
                   organism     = 'dre',
                   pvalueCutoff = 0.05)
        ## TABLE ##
        incProgress(2/5, detail = "dynamic results table...")    
        K.SEA_description = dplyr::select(Kegg.SEA@result, ID, Description)
        K.SEA_description$ID = paste0("<a href=https://www.kegg.jp/kegg-bin/show_pathway?", K.SEA_description$ID," target='_blank'>",K.SEA_description$ID,"</a>")
        K.SEA_value = dplyr::select(Kegg.SEA@result, everything(), -ID, -Description , -GeneRatio, -geneID,-BgRatio) %>% round(4)
        K.SEA = cbind(K.SEA_description, K.SEA_value)
        K.SEAtoSAVE = dplyr::select(K.SEA, everything(), -ID)
        output$K.SEA.Table = DT::renderDataTable({
          K.SEA
        }, escape = F) # escape FALSE to make url
        output$dl.KEGG_SEA <- downloadHandler(
          filename = "KEGGResults_SEA.csv",
          content = function(filename) {
            write.csv(K.SEAtoSAVE, filename, row.names = T)
          }
        )
        ## PLOTS ##
        incProgress(4/5, detail = "Enrichment dot plot...")
        output$dotPlot.KEGG_SEA <- renderPlot({ clusterProfiler::dotplot(Kegg.SEA, showCategory=input$categNb_DP_SEAK) + ggtitle("dotplot for SEA") })
        
        ## PATHVIEW - SEA##
        incProgress(5/5, detail = "Pathway map afterSEA visualization...")
        #selectBox's update with pathway ids
        input$pathwayChoice_SEA
        updateSelectInput(session, "pathwayChoice_SEA", choices = row.names(K.SEA_description), selected = NULL)
        
        observe({
          if (input$pathwayChoice_SEA != "")
          {
            output$pathwayViewer_SEA <- renderImage({
              pathPNG <- pathview(gene.data  = gene,
                                  pathway.id = input$pathwayChoice_SEA,
                                  species = "dre")
              
              list(src = paste0(input$pathwayChoice_SEA, ".pathview.png"),
                   contentType = 'image/png',
                   width = session$clientData$output_pathwayViewer_width*session$clientData$pixelratio*0.7,
                   height = session$clientData$output_pathwayViewer_height*session$clientData$pixelratio
              )
              
            })}})
      })
      
      ### GSEA ###
        ## ANALYSIS
        withProgress(message = 'GSEA ... ', value = 0, {
              Kegg.GSEA <-  gseKEGG(geneList = logF,
              organism     = 'dre',
              nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              # pAdjustMethod = "BH",
               verbose      = FALSE) 

        ## TABLE ##
        incProgress(2/5, detail = "dynamic results table...")    
        K.GSEA_description = dplyr::select(Kegg.GSEA@result, ID, Description)
        K.GSEA_description$ID = paste0("<a href=https://www.kegg.jp/kegg-bin/show_pathway?", K.GSEA_description$ID," target='_blank'>",K.GSEA_description$ID,"</a>")
        K.GSEA_value = dplyr::select(Kegg.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
        K.GSEA = cbind(K.GSEA_description, K.GSEA_value)
        K.GSEAtoSAVE = dplyr::select(K.GSEA, everything(), -ID)
        output$K.GSEA.Table = DT::renderDataTable({
          K.GSEA
        }, escape = F) # escape FALSE to make url
        output$dl.KEGG <- downloadHandler(
          filename = "KEGGResults_GSEA.csv",
          content = function(filename) {
            write.csv(K.GSEAtoSAVE, filename, row.names = T)
          }
        )
        ## PLOTS ##
        incProgress(3/5, detail = "Enrichment dot plot...")
        output$dotPlot.KEGG <- renderPlot({ clusterProfiler::dotplot(Kegg.GSEA, showCategory=input$categNb_DP) + ggtitle("dotplot for GSEA") })
        incProgress(4/5, detail = "Enrichment ridge plot...")
        output$ridgePlot.kegg <- renderPlot({ clusterProfiler::ridgeplot(Kegg.GSEA, showCategory=input$categNb_RP) })
        
        ## Pathway viewer ##

        #selectBox's update with pathway ids
        input$pathwayChoice
        updateSelectInput(session, "pathwayChoice", choices = row.names(K.GSEA_description), selected = NULL)
        
        incProgress(5/5, detail = "Pathway map visualization...")
        observe({
        if (input$pathwayChoice != "")
        {
          # data("demo.paths")
          output$pathwayViewer <- renderImage({
            pathPNG <- pathview(gene.data  = logF,
                                pathway.id = input$pathwayChoice,
                                species = "dre")

            list(src = paste0(input$pathwayChoice, ".pathview.png"),
                 contentType = 'image/png',
                 width = session$clientData$output_pathwayViewer_width*session$clientData$pixelratio*0.7,
                 height = session$clientData$output_pathwayViewer_height*session$clientData$pixelratio
            )

        })}})
      }) # GSEA progress bar end
        
        ########################################################
        ##########  Protein Domains  ##########
        ########################################################
        ## Protein domains DATA
        
        pDomains = read.csv("mart_export.txt") 
        pD.TERM2GENE = dplyr::select(pDomains, Interpro.ID, Gene.stable.ID)
        pD.TERM2GENE$Gene.stable.ID <- mapIds(org.Dr.eg.db, as.vector(pD.TERM2GENE$Gene.stable.ID), 'ENTREZID', 'ENSEMBL')
        # pD.TERM2NAME = unique(pDomains$Interpro.Description) ; names(pD.TERM2NAME) = unique(pDomains$Interpro.ID)
        pD.TERM2NAME = dplyr::select(pDomains, Interpro.ID, Interpro.Description)
        
        gene = pvalues[pvalues<0.05]
        
        ## ANALYSIS ##
        pDomains.SEA = enricher(names(gene), TERM2GENE = pD.TERM2GENE, TERM2NAME = pD.TERM2NAME)
        
        ## TABLE ## 
        proteinDomains.SEA_description = dplyr::select(pDomains.SEA@result, ID, Description)
        proteinDomains.SEA_description$ID = paste0("<a href=https://www.ebi.ac.uk/interpro/entry/InterPro/IPR001064/", proteinDomains.SEA_description$ID," target='_blank'>",proteinDomains.SEA_description$ID,"</a>")
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
        output$dotPlot.pDomains <- renderPlot({ clusterProfiler::dotplot(pDomains.SEA, showCategory=input$categNb_DP.D, ) + ggtitle("dotplot for GSEA") })
        
        ########################################################
        ##########  Motifs  ##########
        ########################################################
        ## MOTIFS (GENE SET) DATA
        motif.geneSet = msigdbr(species = "Danio rerio", category = "C3")
        motif.TERM2GENE = motif.geneSet %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
        
        ### SEA ###
        ## ANALYSIS ##
        withProgress(message = 'SEA ... ', value = 0, {})
        ## PLOTS ##
        
        ## GSEA
        ## ANALYSIS ##
        withProgress(message = 'GSEA ... ', value = 0, {
        motifs.GSEA <- GSEA(logF, TERM2GENE=motif.TERM2GENE, pvalueCutoff = 0.15, nPerm = 10000, pAdjustMethod = "BH")
        head(motifs.GSEA@result$ID) # 0
        
        ## TABLE ##
        Motif.GSEA_description = dplyr::select(motifs.GSEA@result, ID, Description)
        Motif.GSEA_description$ID = paste0("<a href=https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=", Motif.GSEA_description$ID," target='_blank'>",Motif.GSEA_description$ID,"</a>")
        Motif.GSEA_value = dplyr::select(motifs.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
        Motif.GSEA = cbind(Motif.GSEA_description, Motif.GSEA_value)
        Motif.GSEAtoSAVE = dplyr::select(K.GSEA, everything(), -ID)
        
        output$Motif.GSEA.Table = DT::renderDataTable({
          Motif.GSEA
        }, escape = F) # escape FALSE to make url
        
        output$dl.KEGG <- downloadHandler(
          filename = "MotifResults_GSEA.csv",
          content = function(filename) 
          {
            write.csv(Motif.GSEAtoSAVE, filename, row.names = T)
          }
        )
        ## PLOTS ##
        output$dotPlot.Motif <- renderPlot({ dotplot(motifs.GSEA, showCategory=input$categNb_DP.M, ) + ggtitle("dotplot for GSEA") })
        output$ridgePlot.Motif <- renderPlot({ ridgeplot(motifs.GSEA, showCategory=input$categNb_RP.M) })
  }) })
}
