pathways = tabItem(tabName = "pathways",
  h1("Pathways"),
  h2("SEA"),
  ## pvalue
  sliderInput("pv.KEGG", label = h3("p-value threshold : "), min = 0, 
              max = 0.1, value = 0.05),
  fluidRow(
    column(width = 12,
           box(collapsible = TRUE,
               title = span(icon("th-list"),"Enrichment results - SEA"), status = "primary", width = NULL,
               column(width = 10,
                      DT::dataTableOutput("K.SEA.Table") %>% withSpinner()
               ),
               column(width = 2,
                      downloadButton("dl.KEGG_SEA", "Download")
               )
           ),
           box(collapsible = TRUE,
               title = span(icon("sort-amount-up"), "Dot plot - SEA"), status = "primary", width = NULL,
               column(width = 3,
                      numericInput("categNb_DP_SEAK", "Category number : ", 10, min = 1, max = 40, step = 5)
               ),
               column(width = 9,
                      plotOutput("dotPlot.KEGG_SEA") %>% withSpinner(type = 6)
               )
           ),     
    ),
    column(width = 12,
           box(collapsible = TRUE,
               title = span(icon("bezier-curve"), "Pathway Viewer - SEA"), status = "primary", width = NULL,
               column(width = 12,
                      selectInput("pathwayChoice_SEA", "Choose pathway(s) to see", "", selected = NULL, multiple = F,
                                  selectize = TRUE, size = NULL)
               ),
               column(width = 12,
                      imageOutput("pathwayViewer_SEA") %>% withSpinner()
               )
           )
    )
  ),
  
  h2("GSEA"),
  fluidRow(
    column(width = 12,
      box(collapsible = TRUE,
        title = span(icon("th-list"),"Enrichment results - GSEA"), status = "primary", width = NULL,
          DT::dataTableOutput("K.GSEA.Table") %>% withSpinner()
        ),
      box(collapsible = TRUE,
        title = span(icon("sort-amount-up"), "Dot plot - GSEA"), status = "primary", width = NULL,
        column(width = 3,
          numericInput("categNb_DP", "Category number : ", 10, min = 1, max = 40, step = 5)
        ),
        column(width = 9,
          plotOutput("dotPlot.KEGG") %>% withSpinner(type = 6)
          )
        ),     
    ),
    column(width = 12,
      box(collapsible = TRUE,
        title = span(icon("sort-amount-up"), "Ridge plot - GSEA"), status = "primary",width = NULL,
        column(width = 2,
          numericInput("categNb_RP", "Category number : ", 10, min = 1, max = 40, step = 5)
        ),
        column(width = 10,
          plotOutput("ridgePlot.kegg") %>% withSpinner(type = 6)
          )
        ),
    ),
    column(width = 12,
      box(collapsible = TRUE,
        title = span(icon("bezier-curve"), "Pathway Viewer - GSEA"), status = "primary", width = NULL,
        column(width = 12,
          selectInput("pathwayChoice", "Choose pathway(s) to see", "", selected = NULL, multiple = F,
            selectize = TRUE, size = NULL)
          ),
        column(width = 12,
          imageOutput("pathwayViewer") %>% withSpinner()
          )
        )
    )
  )
)
