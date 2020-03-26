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
        
        ################### A FAIRE #################
        
      })
        ## PLOTS ##
      
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
        K.GSEA_description = select(Kegg.GSEA@result, ID, Description)
        K.GSEA_description$ID = paste0("<a href=https://www.kegg.jp/kegg-bin/show_pathway?", K.GSEA_description$ID," target='_blank'>",K.GSEA_description$ID,"</a>")
        K.GSEA_value = select(Kegg.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
        K.GSEA = cbind(K.GSEA_description, K.GSEA_value)
        
        output$K.GSEA.Table = DT::renderDataTable({
          K.GSEA
        }, escape = F) # escape FALSE to make url
        
        ## PLOTS ##
        incProgress(3/5, detail = "Enrichment dot plot...")
        output$dotPlot.KEGG <- renderPlotly({ dotplot(Kegg.GSEA, showCategory=input$categNb_DP) + ggtitle("dotplot for GSEA") })
        incProgress(4/5, detail = "Enrichment ridge plot...")
        output$ridgePlot.kegg <- renderPlot({ ridgeplot(Kegg.GSEA, showCategory=input$categNb_RP) })
        
        ## Pathway viewer ##

        #select box update with pathway ids
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
                                species    = "dre")

            list(src = paste0(input$pathwayChoice, ".pathview.png"),
                 contentType = 'image/png',
                 width = session$clientData$output_pathwayViewer_width*session$clientData$pixelratio*0.7,
                 height = session$clientData$output_pathwayViewer_height*session$clientData$pixelratio
            )

        })}})
      }) # GSEA progress bar finished
        
        ########################################################
        ##########  Motifs  ##########
        ########################################################
        ## MOTIFS (GENE SET) DATA
        motif.geneSet = msigdbr(species = "Danio rerio", category = "C3")
        head(motif.geneSet)
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
        Motif.GSEA_description = select(motifs.GSEA@result, ID, Description)
        Motif.GSEA_description$ID = paste0("<a href=https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=", Motif.GSEA_description$ID," target='_blank'>",Motif.GSEA_description$ID,"</a>")
        Motif.GSEA_value = select(motifs.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
        Motif.GSEA = cbind(Motif.GSEA_description, Motif.GSEA_value)
        
        output$Motif.GSEA.Table = DT::renderDataTable({
          Motif.GSEA
        }, escape = F) # escape FALSE to make url
        output$dotPlot.Motif <- renderPlot({ dotplot(motifs.GSEA, showCategory=input$categNb_DP.M, ) + ggtitle("dotplot for GSEA") })
        output$ridgePlot.Motif <- renderPlot({ ridgeplot(motifs.GSEA, showCategory=input$categNb_RP.M) })
  }) })
}
