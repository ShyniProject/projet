domains <- tabItem(tabName = "domains",
                   h1("Domains"),
                   h2("SEA"),
                   sliderInput("pv.DomainsSEA", label = h3("p-value threshold : "), min = 0, 
                               max = 0.1, value = 0.05),
                   fluidRow(
                     column(width = 12,
                            box(collapsible = TRUE,
                                title =  span(icon("th-list"),"Enrichment results - SEA"), status = "info", width = NULL,
                                column(width = 10,
                                       DT::dataTableOutput("proteinDomains.SEA.Table") %>% withSpinner(color = "#b68f40")
                                ),
                                column(width = 2,
                                       downloadButton("dl.pDomains", "Download")   
                                )
                            ),
                            box(collapsible = TRUE,
                                title = span(icon("sort-amount-up"), "Dot plot - SEA"), status = "info", width = NULL,
                                column(width = 3,
                                       numericInput("categNb_DP.D", "Category number : ", 10, min = 1, max = 40, step = 5),
                                       downloadButton("dl.SEAproteinDomainsDotPlot", "Download")  
                                ),
                                column(width = 9,
                                       plotOutput("dotPlot.pDomains") %>% withSpinner(color = "#b68f40", type = 6)
                                )
                            ),     
                     ),
                   ),
                   
                   
                   
                   h1("Motifs"),
                   h2("GSEA"),
                   fluidRow(
                     column(width = 12,
                            box(collapsible = TRUE,
                                title =  span(icon("th-list"),"Enrichment results - GSEA"), status = "info", width = NULL,
                                column(width = 10,
                                       DT::dataTableOutput("Motif.GSEA.Table") %>% withSpinner(color = "#b6cca1")
                                ),
                                column(width = 2,
                                       downloadButton("dl.Motif", "Download")   
                                )
                            ),
                            box(collapsible = TRUE,
                                title = span(icon("sort-amount-up"), "Dot plot - GSEA"), status = "info", width = NULL,
                                column(width = 3,
                                       numericInput("categNb_DP.M", "Category number : ", 10, min = 1, max = 40, step = 5),
                                       downloadButton("dl.GSEAmotifDotPlot", "Download")  
                                ),
                                column(width = 9,
                                       plotOutput("dotPlot.Motif") %>% withSpinner(color = "#b6cca1", type = 6)
                                )
                            ),     
                     ),
                     column(width = 12,
                            box(collapsible = TRUE,
                                title = span(icon("sort-amount-up"), "Ridge plot - GSEA"), status = "info",width = NULL,
                                column(width = 2,
                                       numericInput("categNb_RP.M", "Category number : ", 10, min = 1, max = 40, step = 5),
                                       downloadButton("dl.SEAmotifRidgePlot", "Download")  
                                ),
                                column(width = 10,
                                       plotOutput("ridgePlot.Motif") %>% withSpinner(color = "#b6cca1", type = 6)
                                )
                            )
                     )
                   )
)