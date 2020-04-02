GO <-tabItem(tabName = "GO",
              h1("Gene Ontology"),
              h2("SEA"),
              sliderInput("pv.GO", label = h3("p-value threshold : "), min = 0, 
                          max = 0.1, value = 0.05),
             fluidRow(
               box(collapsible = TRUE,
                   title = "SEA table", status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO.SEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.SEAGO", "Download")   
                   )
               ),
             ),
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = "Biological Process Plot", status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("SEA_bp") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = "Molecular Function Plot", status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("SEA_mf") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = "Cellular Component Plot", status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("SEA_cc") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 

               )),
              h2("GSEA"),
             fluidRow(
               box(collapsible = TRUE,
                   title = "GSEA table", status = "info", width = NULL,
                   column(width = 10,
                          DT::dataTableOutput("GO.GSEA.Table") %>% withSpinner(color = "#b68f40")
                   ),
                   column(width = 2,
                          downloadButton("dl.GSEAGO", "Download")   
                   )
               ),
             ),
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = "Biological Process Plot", status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("GSEA_bp") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = "Molecular Function Plot", status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("GSEA_mf") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
             fluidRow(
               column(width = 12,
                      box(collapsible = TRUE,
                          title = "Cellular Component Plot", status = "info", width = NULL,
                          column(width = 9,
                                 plotOutput("GSEA_cc") %>% withSpinner(color = "#b68f40", type = 6)
                          )
                      ), 
               )),
)