# Read data
summary = tabItem(tabName = "summary",
                  
                  h1("Put your csv file"),
                  
                  ## FILE INPUT
                  fileInput("data", label = NULL,
                            buttonLabel = "Browse...",
                            placeholder = "No file selected",
                            accept = c('text/csv',
                                       'text/comma-separated-values',
                                       'text/tab-separated-values',
                                       'text/plain',
                                       '.csv',
                                       '.tsv')),
                  useShinyalert(), #if data type not supported only
                  
                  h1("Summary"),
                  fluidRow(
                    column(6,sliderInput('logFcCut', label="log(CPM) cutoff",0,10,2, width="200px")), # Sélection du seuil pour le foldchange
                    column(6,sliderInput('padjCut', label="padj cutoff",0,1,0.05, width="200px")),), 
                  fluidRow(
                    column(width = 12,
                           box(collapsible = TRUE,
                               title = "Volcano plot", status = "info", width = NULL,
                               column(width = 11,
                                      plotOutput("volcanoPlot",click='plot_click.Volcano') %>% withSpinner(color = "#b68f40", type = 6)
                               )
                           ), 
                           box(collapsible = TRUE,
                               title = "Selected point in Volcano Plot", status = "info", width = NULL,
                               column(width = 11,
                                      tableOutput("clickedPoints.Volcano") 
                               )
                           ),
                           
                    ),
                  ),
                  fluidRow(
                    column(width = 12,
                           box(collapsible = TRUE,
                               title = "MA plot", status = "info", width = NULL,
                               column(width = 11,
                                      plotOutput("maPlot",click='plot_click.MA') %>% withSpinner(color = "#b68f40", type = 6)
                               )
                           ), 
                           box(collapsible = TRUE,
                               title = "Selected point in MA Plot", status = "info", width = NULL,
                               column(width = 1,
                                      tableOutput("clickedPoints.MA") 
                               )
                           ),
                           
                    ),
                  )
)