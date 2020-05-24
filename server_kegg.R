########################################################
##########  KEGG  ##########
########################################################
kegg <- function(input, output, session, org, organismsDbKegg, pvalues, logF, minGS, maxGS, nPerm, pvAdjustMethod)
{
  ### SEA ###
  ## ANALYSIS ##
  gene = pvalues[pvalues<0.05]
  Kegg.SEA = enrichKEGG(gene = names(gene),
                        organism = organismsDbKegg[org],
                        keyType = "ncbi-geneid",
                        pvalueCutoff = 0.05,
                        minGSSize = minGS,
                        maxGSSize = maxGS,
                        pAdjustMethod = pvAdjustMethod)
  observe({SEA()})
  SEA <- eventReactive(input$pv.KEGG,{
    withProgress(message = 'SEA ... ', value = 0, {
      gene = pvalues[pvalues<input$pv.KEGG]
      Kegg.SEA = enrichKEGG(gene = names(gene),
                            organism     = organismsDbKegg[org],
                            keyType = "ncbi-geneid",
                            pvalueCutoff = input$pv.KEGG, 
                            minGSSize = minGS,
                            maxGSSize = maxGS,
                            pAdjustMethod = pvAdjustMethod)
      ## TABLE ##
      incProgress(1/3, detail = "dynamic results table...")    
      K.SEA_description = dplyr::select(Kegg.SEA@result, ID, Description)
      K.SEA_description$ID = paste0("<a href=https://www.kegg.jp/kegg-bin/show_pathway?", K.SEA_description$ID," target='_blank'>",K.SEA_description$ID,"</a>")
      K.SEA_value = dplyr::select(Kegg.SEA@result, everything(), -ID, -Description , -GeneRatio, -geneID,-BgRatio) %>% round(5)
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
      incProgress(2/3, detail = "Enrichment dot plot...")
      keggDP <- clusterProfiler::dotplot(Kegg.SEA, showCategory=input$categNb_DP_SEAK) + ggtitle("KEGG enrichment dotplot for SEA") 
      output$dotPlot.KEGG_SEA <- renderPlot({ clusterProfiler::dotplot(Kegg.SEA, showCategory=input$categNb_DP_SEAK) + ggtitle("KEGG enrichment dotplot for SEA")  })
      
      ## DOWNLOAD PLOT
      output$dl.KEGG_dotPlot <- downloadHandler(
        filename = "KEGGResults_SEAdotplot.pdf",
        content = function(file) 
        {
          ggsave(file, plot=keggDP)
        }
      )
      
      ## PATHVIEW - SEA##
      incProgress(3/3, detail = "Pathway map afterSEA visualization...")
      #selectBox's update with pathway ids
      input$pathwayChoice_SEA
      updateSelectInput(session, "pathwayChoice_SEA", choices = row.names(K.SEA_description), selected = NULL)
      
      observe({
        if (input$pathwayChoice_SEA != "")
        {
          output$pathwayViewer_SEA <- renderImage({
            pathPNG <- pathview(gene.data  = gene,
                                pathway.id = input$pathwayChoice_SEA,
                                species = organismsDbKegg[org])
            
            list(src = paste0(input$pathwayChoice_SEA, ".pathview.png"),
                 contentType = 'image/png',
                 width = session$clientData$output_pathwayViewer_width*session$clientData$pixelratio*0.65,
                 height = session$clientData$output_pathwayViewer_height*session$clientData$pixelratio
            )
          })
          ## DOWNLOAD PATHWAY VIEW
          filename <- paste0(input$pathwayChoice_SEA, ".png")
          output$dl.KEGG_pathway <- downloadHandler(
            filename = filename,
            content = function(file) 
            {
              file.copy(paste0(getwd(),'/',input$pathwayChoice_SEA, ".png"), file)
            }
          )
        }})
    })})
  
  ### GSEA ###
  ## ANALYSIS
  withProgress(message = 'GSEA ... ', value = 0, {
    Kegg.GSEA <-  gseKEGG(geneList = logF,
                          organism     = organismsDbKegg[org],
                          keyType = "ncbi-geneid",
                          nPerm        = nPerm,
                          minGSSize    = minGS,
                          maxGSSize    = maxGS,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = pvAdjustMethod,
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
    DPKeggGSEA <- dotplot(Kegg.GSEA, showCategory=input$categNb_DP) + ggtitle("KEGG enrichment dotplot for GSEA")
    output$dotPlot.KEGG <- renderPlot({ dotplot(Kegg.GSEA, showCategory=input$categNb_DP) + ggtitle("KEGG enrichment dotplot for GSEA") })
    
    incProgress(4/5, detail = "Enrichment ridge plot...")
    RPKeggGSEA <- ridgeplot(Kegg.GSEA, showCategory=input$categNb_RP) + ggtitle("KEGG enrichment ridge plot for GSEA")
    output$ridgePlot.kegg <- renderPlot({ ridgeplot(Kegg.GSEA, showCategory=input$categNb_RP) + ggtitle("KEGG enrichment ridge plot for GSEA") })
    
    ## DOWNLOAD GSEA DOTPLOT
    output$dl.KEGG_dotPlotGSEA <- downloadHandler(
      filename = "KEGG_dotPlotGSEA.pdf",
      content = function(file) 
      {
        ggsave(file, plot=DPKeggGSEA)
      }
    )
    
    ## DOWNLOAD GSEA RIDGE PLOT
    output$dl.KEGG_ridgeGSEA <- downloadHandler(
      filename = "KEGG_ridgeGSEA.pdf",
      content = function(file) 
      {
        ggsave(file, plot=RPKeggGSEA)
      }
    )
    
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
                              species = organismsDbKegg[org])
          
          list(src = paste0(input$pathwayChoice, ".pathview.png"),
               contentType = 'image/png',
               width = session$clientData$output_pathwayViewer_width*session$clientData$pixelratio*0.7,
               height = session$clientData$output_pathwayViewer_height*session$clientData$pixelratio
          )
        })
        ## DOWNLOAD PATHWAY VIEW GSEA
        filename <- paste0(input$pathwayChoice, ".png")
        output$dl.KEGG_pathwayGSEA <- downloadHandler(
          filename = filename,
          content = function(file) 
          {
            file.copy(paste0(getwd(),'/',input$pathwayChoice, ".png"), file)
          }
        )      
    }})
  }) # GSEA progress bar end
}