summary <- function(input, output, session, RNAseqDF)
{
  ########################################################
  ##########        Volcano & MA plots          ##########
  ########################################################
  source(file = "MA_Volcano_plot.R")
  
  d <- RNAseqDF
  observe({pv.GO()})
  pv.GO <- eventReactive(c(input$logFcCut,input$padjcut),{
    withProgress(message = 'Volcano & MA Plots ... ', value = 0, {
      dataFiltered <- data.filtering(input =  d,
                                     logFcCut = input$logFcCut,
                                     padjCut = input$padjCut) # Filtre et travaille sur les donnees 
      
      incProgress(1/2, detail = "Volcano Plot...")
      
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
      
      incProgress(2/2, detail = "MA Plot...")
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
}