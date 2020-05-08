library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library(shinyalert)# for file format checking
library(ggplot2)
# library(plotly) # to make dynamic ggplots
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
library(PPInfer)

source("ui_summary.R")
source("ui_kegg.R")
source("ui_go.R")
source("ui_proteinDomain.R")

# source("summary.r")
# source("pathways.r")
# source("domains.r")
# source("GO_item.R")

dashboardPage(
    dashboardHeader(title = "EASE"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Whole data inspection", tabName = "summary", icon = icon("poll")),
            menuItem("Pathways", tabName = "pathways", icon = icon("project-diagram")),
            menuItem("GO terms", tabName = "GO", icon = icon("sitemap")),
            menuItem("Protein Domains", tabName = "domains", icon = icon("bullseye")))
        
        
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