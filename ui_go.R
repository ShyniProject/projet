GO <-tabItem(tabName = "GO",
              h1("Gene Ontology"),
             
             ########## Slider to SEA p-value threshold numeric input 
              h2("SEA"),
              sliderInput("pv.GO", label = h3("p-value threshold : "), min = 0, 
                          max = 0.1, value = 0.05),
             
             ########## Biologicial processes - GSEA
             ## TABLE
             fluidRow(
               box(collapsible = TRUE,
                   title = span(icon("th-list"), "SEA table - Biological Process"), status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO_BP.SEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.SEAGO_BP", "Download")   
                   )
               ),
             ),
             
             ## PLOT
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = span(icon("sort-amount-up"),"Biological Process Plot"), status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("SEA_bp") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             
             ########## Molecular functions - SEA
             ## TABLE
             fluidRow(
               box(collapsible = TRUE,
                   title = span(icon("th-list"), "SEA table : Molecular functions"), status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO_MF.SEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.SEAGO_MF", "Download")   
                   )
               ),
             ),
             
             ## PLOT
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = span(icon("sort-amount-up"),"Molecular Function Plot"), status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("SEA_mf") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             
             ########## Cellular Components - SEA
             ## TABLE
             fluidRow(
               box(collapsible = TRUE,
                   title = span(icon("th-list"), "SEA table : Cellular Components"), status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO_CC.SEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.SEAGO_CC", "Download")   
                   )
               ),
             ),
             
             ## PLOT
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = span(icon("sort-amount-up"),"Cellular Component Plot"), status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("SEA_cc") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 

               )),
              h2("GSEA"),
              fluidRow(
              ########## Biologicial processes - SEA
              ## TABLE
               box(collapsible = TRUE,
                   title = span(icon("th-list"),"GSEA table : Biological Functions"), status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO_BP.GSEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.GSEAGO_BP", "Download")   
                   )
               ),
             ),
             ## PLOT
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = span(icon("sort-amount-up"),"Biological Process Plot"), status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("GSEA_bp") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             
             ########## Molecular Functions - GSEA
             ## TABLE
             fluidRow(
              box(collapsible = TRUE,
                 title = span(icon("th-list"),"GSEA table : Molecular Functions"), status = "info", width = NULL,
                 column(width = 10,
                        DT::dataTableOutput("GO_MF.GSEA.Table") %>% withSpinner(color = "#b68f40")
                 ),
                 column(width = 2,
                        downloadButton("dl.GSEAGO_MF", "Download")   
                        )
                ),
              ),
              ## PLOT
              fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = span(icon("sort-amount-up"),"Molecular Functions Plot"), status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("GSEA_mf") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),

             ########## Cellular Components - GSEA
             ## TABLE
             fluidRow(
               box(collapsible = TRUE,
                   title = span(icon("th-list"), "GSEA table : Cellular Components"), status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO_CC.GSEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.GSEAGO_CC", "Download")   
                   )
               ),
             ),
             ## PLOT
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = span(icon("sort-amount-up"),"Cellular Components Plot"), status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("GSEA_cc") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
)