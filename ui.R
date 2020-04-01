library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyalert) # for file format checking
library(ggplot2)
library(plotly) # to make dynamic ggplots
library(DT) ## datatable results
library(clusterProfiler) # gsea analysis
library(enrichplot) #plots after gsea
library(org.Dr.eg.db) # mapIds
library(KEGG.db)
library(KEGGREST)
library(msigdbr) # motifs
library(dplyr)
library(pathview) # kegg pathways png
library(DOSE)


source("summary.r")
source("pathways.r")
source("domains.r")



##############################Go Panel #############################


GO <- tabItem(tabName = "GO",
               h1("Gene Ontology"),
               h2("SEA"),
               h3("Biological Process"),
               plotOutput('SEA_bp'),
               h3("Molecular Function"),
               plotOutput('SEA_mf'),
               h3("Cellular Component"),
               plotOutput('SEA_cc'),
               h2("GSEA"),
               h3("Biological Process"),
               plotOutput('GSEA_bp'),
               h3("Molecular Function"),
               plotOutput('GSEA_mf'),
               h3("Cellular Component"),
               plotOutput('GSEA_cc'),
)

dashboardPage(
    dashboardHeader(title = "Brillant"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Summary", tabName = "summary"),
            menuItem("Pathways", tabName = "pathways", icon = icon("poll")),
            menuItem("Domains", tabName = "domains"),
            menuItem("GO", tabName = "GO"))
        
    ),
    dashboardBody(
        tabItems(
            summary,
            pathways,
            domains,
            GO
        )
    )
)
