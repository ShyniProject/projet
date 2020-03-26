domains <- tabItem(tabName = "domains",
                   h2("Domains"),
                   h3("SEA"),
                   h3("GSEA"),
                   
                   h2("Motifs"),
                   h3("SEA"),
                   sliderInput("pv.Domain", label = h3("p-value threshold : "), min = 0, 
                               max = 0.1, value = 0.05),
      
                   h3("GSEA"),
                   fluidRow(
                     column(width = 12,
                            box(collapsible = TRUE,
                                title = "Enrichment results - GSEA", status = "info", width = NULL,
                                DT::dataTableOutput("Motif.GSEA.Table") %>% withSpinner(color = "#14D6E2")
                            ),
                            #   ),),
                            # fluidRow(
                            # column(width = 12,
                            box(collapsible = TRUE,
                                title = "Dot plot - GSEA", status = "info", width = NULL,
                                column(width = 3,
                                       numericInput("categNb_DP.M", "Category number : ", 10, min = 1, max = 40, step = 5)
                                ),
                                column(width = 9,
                                       plotOutput("dotPlot.Motif") %>% withSpinner(color = "#14D6E2", type = 6)
                                )
                            ),     
                     ),
                     column(width = 12,
                            box(collapsible = TRUE,
                                title = "Ridge plot - GSEA", status = "info",width = NULL,
                                column(width = 2,
                                       numericInput("categNb_RP.M", "Category number : ", 10, min = 1, max = 40, step = 5)
                                ),
                                column(width = 10,
                                       plotOutput("ridgePlot.Motif") %>% withSpinner(color = "#14D6E2", type = 6)
                                )
                            )
  )
