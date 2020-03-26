library(shiny)
library(shinydashboard)
library(shinycssloaders)
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

dashboardPage(
    dashboardHeader(title = "Brillant"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Summary", tabName = "summary"),
            menuItem("Pathways", tabName = "pathways", icon = icon("poll")),
            menuItem("Domains", tabName = "domains"))
        
    ),
    dashboardBody(
        tabItems(
            summary,
            pathways,
            domains
        )
    )
)