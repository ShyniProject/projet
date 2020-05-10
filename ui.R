# ===============================================
# Uncomment lines below to install all packages 
# ===============================================
# install.packages("anyLib")
# anyLib.packages=c("shiny", "shinydashboard", "shinydashboardPlus", "shinycssloaders", "shinyalert",
#                   "ggplot2", "DT", "clusterProfiler", "enrichplot", "org.Dr.eg.db", "KEGG.db", "KEGGREST",
#                   "msigdbr", "dplyr", "pathview", "DOSE", "PPInfer")
# anyLib::anyLib(anyLib.packages, autoUpdate = T)
# ===============================================

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library(shinyalert)# for file format checking
library(ggplot2)
library(DT) ## datatable results
library(clusterProfiler) # gsea analysis
library(enrichplot) #plots after gsea
library(KEGG.db)
library(KEGGREST)
library(msigdbr) # motifs
library(dplyr)
library(pathview) # kegg pathways png
library(DOSE)
library(PPInfer)
## annotation files for mapID()
library(org.Dr.eg.db) 
library(org.Dm.eg.db) 


source("ui_summary.R")
source("ui_kegg.R")
source("ui_go.R")
source("ui_proteinDomain.R")

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
            GO,
            domains
        )
    )
)