#######################################################
## Project: Fonction shiny project Volcano et MA plot #
## Script purpose: Tester les fonction de génération  #
## du volcano et du MA plot                           #
## Date: 23/03/2020                                   #
## Author: Arnoux Jérôme                              #
#######################################################


ui <- fluidPage(
  h1("Clickable Volcano Plots!"),
  sliderInput('cpmCut', label="log(CPM) cutoff",0,10,2, width="200px"), # Sélection du seuil pour le foldchange
  sliderInput('padjCut', label="padj cutoff",0,1,0.05, width="200px"),  # Sélection du seuil de la pvalue ajusté
  plotOutput('volcanoPlot',click='plot_click.Volcano'), #VolcanoPlot
  tableOutput('clickedPoints.Volcano'), # Tableau correspondant à la zone cliqué du VolcanoPlot
  plotOutput('maPlot',click='plot_click.MA'), # MA plot
  tableOutput('clickedPoints.MA') # Tableau correspondant à la zone cliqué du MA plot
  )

server <- function(input, output) {
  dataFrame<- reactive({
    read_csv(file="https://raw.githubusercontent.com/ShyniProject/Rendu/master/GSE129081_small.csv", col_names = TRUE)
  })
  dataFilter <- reactive({
    data.filtering(dataFrame()) # Filtre et travail sur les données 
  })
  
  output$volcanoPlot <- renderPlot({ 
    VolcanoPlot(dataFilter())
  })
  clicked.Volcano <- reactive({
    nearPoints(dataFilter(), input$plot_click.Volcano, xvar = "log2FoldChange", yvar = "negLogpadj")
  }) #Prend en compte un nombre de point proche du pointeur 
  
  output$clickedPoints.Volcano <- renderTable({
    clicked.Volcano()
  }, rownames = T)
  
  output$maPlot <- renderPlot({
    MaPlot(dataFilter())
  })
  clicked.MA <- reactive({
    nearPoints(dataFilter(), input$plot_click.MA, xvar = "logBaseMean", yvar = "log2FoldChange")
  })

  #output those points into a table
  output$clickedPoints.MA <- renderTable({
    clicked.MA()
  }, rownames = T)
}

shinyApp(ui, server)