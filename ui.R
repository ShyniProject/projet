library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyalert)# for file format checking
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
source("GO_item.R")


dashboardPage(
    dashboardHeader(title = "EASE"),
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