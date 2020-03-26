#######################################################
## Project: Fonction shiny project Volcano et MA plot #
## Script purpose: Créer des fonctions pour filtrer   #
## et modifier les données d'entrér pour générer un   #
## Date: 23/03/2020                                   #
## Author: Arnoux Jérôme                              #
#######################################################

##Fonction 
data.filtering=function(input,logFcCut,padjCut){
  mutate.input=input %>% 
    mutate(SigFC=ifelse(abs(log2FoldChange)>logFcCut, 1, 0)) %>%
    mutate(SigPadj=ifelse(padj<padjCut,1,0)) %>%
    mutate(Legendary=ifelse(SigFC==1,
                      ifelse(SigPadj==1, #SigFC vaut 1
                             paste0("logFC>",logFcCut," and padj<",padjCut),
                             paste0("logFC>",logFcCut)),
                      ifelse(SigPadj==1, #SigFC=0
                             paste0("padj<",padjCut),
                             paste0("logFC<",logFcCut," and padj>",padjCut)))) %>%
    mutate(negLogpadj=-log10(padj)) %>%
    mutate(logBaseMean=log(baseMean,2)) %>%
    select(id:Orthologous_human_gene,Sig:logBaseMean)
  return(mutate.input)
}
View(data.filtering(data,2,0.05))
#VolcanoPlot
VolcanoPlot=function(data){
  ggplot(data,aes(x=log2FoldChange, y=negLogpadj, color=Sig)) +
    geom_point() +
    coord_cartesian() +
    ylab("-log10 padj") +
    xlab("log2 fold change")
}
#MAPlot
MaPlot=function(data){
  ggplot(data,aes(x=logBaseMean, y=log2FoldChange, color=Sig)) +
    geom_point() +
    coord_cartesian() +
    ylab("log2 FC") +
    xlab("log2 CPM")
}


## Execution de test

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")
# BiocManager::install("SummarizedExperiment")

# library(edgeR)
# library(SummarizedExperiment)
# library(readr)
# library(dplyr)
# library(ggplot2)
#get the data
#data <- read_csv(file="https://raw.githubusercontent.com/ShyniProject/Rendu/master/GSE129081_small.csv",col_names = TRUE)
# head(data)
# # A tibble: 6 x 8
# id     Gene   baseMean log2FoldChange    pvalue      padj Biotype Orthologous_hum…
# <chr>  <chr>     <dbl>          <dbl>     <dbl>     <dbl> <chr>   <chr>           
#   1 ENSDA… si:dk…    1278.          -2.29 7.74e-105 1.58e-100 protei… COL28A1         
# 2 ENSDA… pdxka      353.          -2.44 4.81e- 66 4.90e- 62 protei… PDXK            
# 3 ENSDA… 5S_rR…     284.          -2.51 9.64e- 47 6.55e- 43 rRNA    NA              
# 4 ENSDA… si:dk…    1310.          -1.33 6.15e- 33 3.13e- 29 protei… NA              
# 5 ENSDA… rabl3      274.          -1.36 4.10e- 29 1.67e- 25 protei… RABL3           
# 6 ENSDA… elovl6    2588.          -1.59 7.64e- 26 2.59e- 22 protei… ELOVL6    

# dataFiltered=data.filtering(data)
# VolcanoPlot(dataFiltered)
# MaPlot(dataFiltered)