function(input, output, session) {
  
  observe({
    analysis()
  })
  
  analysis <- eventReactive(input$data,{
    
    ########################################################
    ##########  data type check  ##########
    ########################################################
    formats = c(
      'text/csv',
      'text/comma-separated-values',
      'text/tab-separated-values',
      'text/plain',
      'csv',
      'tsv',
      "application/vnd.ms-excel"
    )
    
    if(!input$data$type %in% formats)
    {
      shinyalert("Oops!", paste0("It seems you have uploaded a ",input$data$type," file \n Please refresh page."), type = "warning", animation = T)
    } 
    
    
    
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
      observe({SEA()})
      SEA <- eventReactive(input$pv.KEGG,{
        gene = pvalues[pvalues<input$pv.KEGG]
        Kegg.SEA = enrichKEGG(gene = names(gene),
                              organism     = 'dre',
                              pvalueCutoff = input$pv.KEGG)
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
        
        ## DOWNLOAD
        output$dl.KEGG_SEA <- downloadHandler(
          filename = "KEGGResults_SEA.csv",
          content = function(filename) {
            write.csv(K.SEAtoSAVE, filename, row.names = T)
          }
        )
        ## PLOTS ##
        incProgress(4/5, detail = "Enrichment dot plot...")
        output$dotPlot.KEGG_SEA <- renderPlot({ clusterProfiler::dotplot(Kegg.SEA, showCategory=input$categNb_DP_SEAK) + ggtitle("dotplot for SEA") })
        
        ## PATHVIEW - SEA##
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
      })})
    
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
      K.GSEA_description$ID = paste0("<a href=https://www.kegg.jp/kegg-bin/show_pathway?",
                                     K.GSEA_description$ID," target='_blank'>",K.GSEA_description$ID,"</a>")
      K.GSEA_value = dplyr::select(Kegg.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
      K.GSEA = cbind(K.GSEA_description, K.GSEA_value)
      
      output$K.GSEA.Table = DT::renderDataTable({
        K.GSEA
      }, escape = F) # escape FALSE to make url
      
      ## PLOTS ##
      incProgress(3/5, detail = "Enrichment dot plot...")
      output$dotPlot.KEGG <- renderPlot({ dotplot(Kegg.GSEA, showCategory=input$categNb_DP) + ggtitle("dotplot for GSEA") })
      incProgress(4/5, detail = "Enrichment ridge plot...")
      output$ridgePlot.kegg <- renderPlot({ ridgeplot(Kegg.GSEA, showCategory=input$categNb_RP) })
      
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
                                species    = "dre")
            
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
    
    pDomains = read.csv("mart_export_INTERPRO.txt") 
    pD.TERM2GENE = dplyr::select(pDomains, Interpro.ID, Gene.stable.ID)
    pD.TERM2GENE$Gene.stable.ID <- mapIds(org.Dr.eg.db, as.vector(pD.TERM2GENE$Gene.stable.ID), 'ENTREZID', 'ENSEMBL')
    # pD.TERM2NAME = unique(pDomains$Interpro.Description) ; names(pD.TERM2NAME) = unique(pDomains$Interpro.ID)
    pD.TERM2NAME = dplyr::select(pDomains, Interpro.ID, Interpro.Description)
    
    gene = pvalues[pvalues<0.05]
    
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
      head(motifs.GSEA@result$ID) 
      
      ## TABLE ##
      Motif.GSEA_description = dplyr::select(motifs.GSEA@result, ID, Description)
      Motif.GSEA_description$ID = paste0("<a href=https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=", Motif.GSEA_description$ID," target='_blank'>",Motif.GSEA_description$ID,"</a>")
      Motif.GSEA_value = dplyr::select(motifs.GSEA@result, setSize:p.adjust, -NES) %>% round(3)
      Motif.GSEA = cbind(Motif.GSEA_description, Motif.GSEA_value)
      
      output$Motif.GSEA.Table = DT::renderDataTable({
        Motif.GSEA
      }, escape = F) # escape FALSE to make url
      
      ## PLOTS ##
      output$dotPlot.Motif <- renderPlot({ dotplot(motifs.GSEA, showCategory=input$categNb_DP.M, ) + ggtitle("dotplot for GSEA") })
      output$ridgePlot.Motif <- renderPlot({ ridgeplot(motifs.GSEA, showCategory=input$categNb_RP.M) })
      
      
    })
    ########################################################
    ##########        Volcano & MA plots          ##########
    ########################################################
    source(file = "MA_Volcano_plot.R")
    
    withProgress(message = 'Volcano & MA Plots ... ', value = 0, {
      observe({pv.GO()})
      pv.GO <- eventReactive(c(input$logFcCut,input$padjcut),{
        dataFiltered <- data.filtering(input =  d,
                                     logFcCut = input$logFcCut,
                                     padjCut = input$padjCut) # Filtre et travail sur les donnees 

      incProgress(2/3, detail = "Volcano Plot...")
      
      output$volcanoPlot <- renderPlot({ 
        VolcanoPlot(dataFiltered)
      })
      
      observe({clickVolc()})
      clickVolc <- eventReactive(input$plot_click.Volcano,{
        
        
        clicked.Volcano <-  nearPoints(dataFiltered, input$plot_click.Volcano, xvar = "log2FoldChange", yvar = "negLogpadj")
        #Prend en compte un nombre de point proche du pointeur 
        
        output$clickedPoints.Volcano <- renderTable({
          clicked.Volcano
        }, rownames = T)})
      
      incProgress(3/3, detail = "MA Plot...")
      output$maPlot <- renderPlot({
        MaPlot(dataFiltered)
      })
      
      
      observe({clickMA()})
      clickMA <- eventReactive(input$plot_click.MA,{
        clicked.MA <- nearPoints(dataFiltered, input$plot_click.MA, xvar = "logBaseMean", yvar = "log2FoldChange")
        
        #output those points into a table
        output$clickedPoints.MA <- renderTable({clicked.MA}, rownames = T)
        
      })
      })
    })
    
    ########################################################
    ##########                GO                  ##########
    ########################################################
    source(file = "GO_term.R")
    
    withProgress(message = 'GO ... ', value = 0, {
      dataFilteredGO <-  DataFilterGO(d)
      
      ########################## GO SEA #######################################
      incProgress(2/5, detail = "SEA analysis...")
      observe({SEA_GO()})
      SEA_GO <- eventReactive(input$pv.GO,{
      SEA_result <- SEAanalysis(dataFilteredGO,input$pv.GO)
      #### Plots ####
      incProgress(3/5, detail = "SEA plots...")
      
      output$SEA_bp <- renderPlot({bp_SEA_Plot(SEA_result)})
      
      
      output$SEA_mf <- renderPlot({
        mf_SEA_Plot(SEA_result)
      })
      
      output$SEA_cc <- renderPlot({
        cc_SEA_Plot(SEA_result)
      })
      
      # cats <- c("BP", "CC", "MF")
      ## BP TABLE ##
      GO_BP.SEA_description = SEA_result$BP
      GO_BP.SEA_description$GO_term= paste0("<a href=http://amigo.geneontology.org/amigo/term/", GO_BP.SEA_description$GO_term," target='_blank'>",GO_BP.SEA_description$GO_term, "</a>")
      GO_BP.SEA_description$pval = GO_BP.SEA_description$pval %>% round(3)
      output$GO_BP.SEA.Table = DT::renderDataTable({
        GO_BP.SEA_description
      }, escape = F) # escape FALSE t
      ## DOWNLOAD
      output$dl.SEAGO_BP <- downloadHandler(
        filename = "SEA_GO_BP.csv",
        content = function(filename) {
          write.csv(GO_BP.SEA_description, filename, row.names = T)
        }
      )
      ## CC TABLE ##
      GO_CC.SEA_description = SEA_result$CC
      GO_CC.SEA_description$GO_term= paste0("<a href=http://amigo.geneontology.org/amigo/term/", GO_CC.SEA_description$GO_term," target='_blank'>",GO_CC.SEA_description$GO_term, "</a>")
      GO_CC.SEA_description$pval = GO_CC.SEA_description$pval %>% round(3)
      output$GO_CC.SEA.Table = DT::renderDataTable({
        GO_CC.SEA_description
      }, escape = F) # escape FALSE t   
      ## DOWNLOAD
      output$dl.SEAGO_CC <- downloadHandler(
        filename = "SEA_GO_CC.csv",
        content = function(filename) {
          write.csv(GO_CC.SEA_description, filename, row.names = T)
        }
      )
      ## MF TABLE ##
      GO_MF.SEA_description = SEA_result[[1]]
      GO_MF.SEA_description$GO_term= paste0("<a href=http://amigo.geneontology.org/amigo/term/", GO_MF.SEA_description$GO_term," target='_blank'>",GO_MF.SEA_description$GO_term, "</a>")
      GO_MF.SEA_description$pval = GO_MF.SEA_description$pval %>% round(3)
      
      output$GO_MF.SEA.Table = DT::renderDataTable({
        GO_MF.SEA_description
      }, escape = F) # escape FALSE to make url
      
      ## DOWNLOAD
      output$dl.SEAGO_MF <- downloadHandler(
        filename = "SEA_GO_MF.csv",
        content = function(filename) {
          write.csv(GO_MF.SEA_description, filename, row.names = T)
        }
      )
  }) # end GO reactive
      

      
      ########################## GO GSEA #######################################
      
      
      incProgress(4/5, detail = "GSEA analysis...")
      GSEA_result <- GSEAanalysis(dataFilteredGO)
      
      
      #### Plots ####
      incProgress(5/5, detail = "GSEA Plots...")
      
      output$GSEA_bp <- renderPlot({
        bp_GSEA_Plot(GSEA_result)
      })
      
      output$GSEA_mf <- renderPlot({
        mf_GSEA_Plot(GSEA_result)
      })
      
      output$GSEA_cc <- renderPlot({
        cc_GSEA_Plot(GSEA_result)
      })
      
      ## TABLE - Biological Processes ##
      GO.GSEA_description = GSEA_result[["BP"]]@result
      GO.GSEA_description = dplyr::select(GSEA_result[["BP"]]@result, ID, Description)
      GO.GSEA_description$ID = paste0("<a href=http://amigo.geneontology.org/amigo/term/", GO.GSEA_description$ID," target='_blank'>",GO.GSEA_description$ID,"</a>")
      GO.GSEA_value = dplyr::select(GSEA_result[["BP"]]@result, setSize:p.adjust, -NES) %>% round(3)
      GO.GSEA = cbind(GO.GSEA_description, GO.GSEA_value)

      output$GO_BP.GSEA.Table = DT::renderDataTable({
        GO.GSEA
      }, escape = F) # escape FALSE to make url
      
      ## DOWNLOAD
      output$dl.GSEAGO_BP <- downloadHandler(
        filename = "GSEA_GO_BP.csv",
        content = function(filename) {
          write.csv(GSEA_result[["BP"]]@resul, filename, row.names = T)
        }
    )
      
      ## TABLE - Molecular Functions ##
      GO.GSEA_description = GSEA_result[["MF"]]@result
      GO.GSEA_description = dplyr::select(GSEA_result[["MF"]]@result, ID, Description)
      GO.GSEA_description$ID = paste0("<a href=http://amigo.geneontology.org/amigo/term/", GO.GSEA_description$ID," target='_blank'>",GO.GSEA_description$ID,"</a>")
      GO.GSEA_value = dplyr::select(GSEA_result[["MF"]]@result, setSize:p.adjust, -NES) %>% round(3)
      GO.GSEA = cbind(GO.GSEA_description, GO.GSEA_value)

      output$GO_MF.GSEA.Table = DT::renderDataTable({
        GO.GSEA
      }, escape = F) # escape FALSE to make url
      
      ## DOWNLOAD
      output$dl.GSEAGO_MF <- downloadHandler(
        filename = "GSEA_GO_MF.csv",
        content = function(filename) {
          write.csv(GSEA_result[["MF"]]@resul, filename, row.names = T)
        }
    )
      
      ## TABLE - Cellular Components ##
      GO.GSEA_description = GSEA_result[["CC"]]@result
      GO.GSEA_description = dplyr::select(GSEA_result[["CC"]]@result, ID, Description)
      GO.GSEA_description$ID = paste0("<a href=http://amigo.geneontology.org/amigo/term/", GO.GSEA_description$ID," target='_blank'>",GO.GSEA_description$ID,"</a>")
      GO.GSEA_value = dplyr::select(GSEA_result[["CC"]]@result, setSize:p.adjust, -NES) %>% round(3)
      GO.GSEA = cbind(GO.GSEA_description, GO.GSEA_value)

      output$GO_CC.GSEA.Table = DT::renderDataTable({
        GO.GSEA
      }, escape = F) # escape FALSE to make url
      
      ## DOWNLOAD
      output$dl.GSEAGO_CC <- downloadHandler(
        filename = "GSEA_GO_CC.csv",
        content = function(filename) {
          write.csv(GSEA_result[["CC"]]@resul, filename, row.names = T)
        }
    )
      
  })
  } 
)}
