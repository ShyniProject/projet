## app.R ##
library(shiny)
library(shinydashboard)
library(shinyWidgets)

sidebar <- dashboardSidebar(
    sidebarMenu(
        h3("Paramètres", align = "center", color = "#000000"),
        fileInput("file", "Votre fichier",
                  buttonLabel = "Importer...",
                  placeholder = "Aucun fichier",
                  accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),   
        sliderInput("slider", "Nombre de lignes", 1, 100, 50),
        radioGroupButtons(
            inputId = "method", label = "Analysis :", 
            choices = c("SEA", "GSEA"), 
            justified = TRUE, status = "primary",
            checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
        ),
        hr(),
        fluidRow(column(width = 6,
                            column(6, offset = 5,
                                   # put all your filters in here
                                  # actionButton(inputId ="foo", label = icon("sign-in-alt", lib="font-awesome"), color = "success",class = "btn-secondary", width = NULL)
                                  actionBttn(
                                      inputId = "bttn1",
                                      label = "Go!",
                                      color = "royal",
                                      style = "material-flat",
                                      icon = icon("rocket", lib = "font-awesome"),
                                    ),
                                  )
                        )
                )
    )
)

Accueil <- tabPanel("Accueil",
         h1("Analyses rapides"),
         h2("Aperçu"),
)
GO <- tabPanel("GO",
         h1("Gene Ontology"),
         h2("Aperçu"),
         h2("Effectifs par domaine")
)
Pathways <- tabPanel("Pathways",
         h1("Header 1")
)
Domain <- tabPanel("Domains",
         h1("Header 1")
)

body <- dashboardBody(
    tabsetPanel(
        Accueil,
        GO,
        Pathways,
        Domain
    ),
    tags$head(tags$style(HTML('
    /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #f1f1f1;
                              color:black;
        }
        
    /* main sidebar */
        .sidebar-menu {
                              color:black;
        }
    /* body */
        .content-wrapper, .right-side {
        background-color: #ffffff;
        }                              
    ')))
)

dashboardPage(
    dashboardHeader(title = "Test Brillant"),
    sidebar,
    body
)

