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
                                    em("For input format, please refer to github's README available at\n"),
                                    appButton(
                                      url = "https://github.com/ShinyProject/projet",
                                      label = "Our Github",
                                      icon = "fa fa-github",
                                      enable_badge = TRUE,
                                      badgeColor = "red",
                                      badgeLabel = "new"
                                    ),
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
                                    column(width = 7, align = "center", offset = 2,
                                      h1("Choose organism\n"),
                                      selectInput("organismDb",icon("paw"), choices = NULL, selected = NULL)
                                    ),
                                    column(width = 7, align = "center", offset = 2,
                                      h1("Choose ID database origin\n"),
                                      selectInput("geneidDb",icon("database"), choices = NULL, selected = NULL)
                                    ),
                                    column(width = 7, align = "center", offset = 2,
                                      h1("Enrichment parameters\n"),
                                        box(collapsible = TRUE, collapsed = TRUE,
                                               title = span(icon("gear")), status = "primary", width = NULL,
                                               column(width = 6,
                                                      numericInput("minGSSize", "Minimal number of genes in pathways:", value = 5, min = 1, max = 5000, step = 1)
                                               ),
                                               column(width = 6,
                                                      numericInput("maxGSSize", "Maximal number of genes in pathways:", value = 800, min = 1, max = 5000, step = 1)
                                               ),
                                               column(width = 6,
                                                      numericInput("nPerm", "Permutation number for GSEA's p-values:", value = 10000, min = 1, max = 5000, step = 1)
                                               ),
                                               column(width = 6,
                                                      selectInput("pvAdjust","p-value correction method:", choices = NULL, selected = NULL)
                                               ),
                                        ),
                                        actionBttn(
                                          inputId = "start",
                                          label = "Go!",
                                          color = "success",
                                          style = "material-flat",
                                          icon = icon("rocket", lib = "font-awesome"),
                                      ),
                                    )
                    )),
                    useShinyalert(), #if data type not supported only
                    footer = NULL,
                    footer_padding = F,
                    width = 12
                  )),
                
                
                h1("Whole data inspection"),
                fluidRow(
                  column(6,sliderInput('logFcCut', label="log(CPM) cutoff",0,10,2, width="200px")), # S?lection du seuil pour le foldchange
                  column(6,sliderInput('padjCut', label="padj cutoff",0,1,0.05, width="200px")),), 
                fluidRow(
                  column(width = 12,
                         box(collapsible = TRUE,
                             title = "Volcano plot", status = "info", width = NULL,
                             column(width = 11,
                                    plotOutput("volcanoPlot",click='plot_click.Volcano') %>% withSpinner(color = "#b68f40", type = 6),
                                    downloadButton("dl.volcanoPlot", "Download")
                             )
                         ), 
                         box(collapsible = TRUE,
                             title = "Selected point in Volcano Plot", status = "info", width = NULL,
                             column(width = 11,
                                    tableOutput("clickedPoints.Volcano"),
                             )
                         ),
                         
                  ),
                ),
                fluidRow(
                  column(width = 12,
                         box(collapsible = TRUE,
                             title = "MA plot", status = "info", width = NULL,
                             column(width = 11,
                                    plotOutput("maPlot",click='plot_click.MA') %>% withSpinner(color = "#b68f40", type = 6),
                                    downloadButton("dl.maPlot", "Download")
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