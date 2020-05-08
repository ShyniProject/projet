# Read data
summary = tabItem(tabName = "summary",
                  
                fluidRow(
                  widgetUserBox(
                    title = strong("Enrichment Analysis in Shiny Environment"),
                    subtitle = strong("EASE"),
                    type = NULL,
                    src = "https://www.flaticon.com/premium-icon/icons/svg/1882/1882123.svg",
                    color = "yellow",
                    # background = TRUE,
                    # backgroundUrl = "https://image.flaticon.com/icons/svg/307/307757.svg",
                    closable = F,
                    fluidRow(column(width = 12, align = "center", 
                                    h2("Welcome to our enrichment analysis shiny application for RNA-seq data !"),
                                    h1("Put your .csv file\n"),
                                    em("csv or text/csv files only with \"id\" for ENSEMBL genes ID, \"padj\" for p-values adjusted and \"log2FoldChange\" columns.\n\n"),
                                    ## FILE INPUT
                                    column(width = 7, align = "center", offset = 2,
                                           fileInput("data", label = NULL,
                                                     buttonLabel = "Browse...",
                                                     placeholder = "No file selected",
                                                     accept = c('text/csv',
                                                                'text/comma-separated-values',
                                                                'text/tab-separated-values',
                                                                'text/plain',
                                                                '.csv',
                                                                '.tsv'))),
                                    h1("Choose organism\n"),
                    )),
                    useShinyalert(), #if data type not supported only
                    footer = NULL,
                    footer_padding = F,
                    width = 12
                  )),
                
                selectInput("organismDb","Organism:", choices = NULL, selected = NULL),
                
                h1("Whole data inspection"),
                fluidRow(
                  column(6,sliderInput('logFcCut', label="log(CPM) cutoff",0,10,2, width="200px")), # S?lection du seuil pour le foldchange
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