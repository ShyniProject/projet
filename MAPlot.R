###### MA Plot

#library
library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library("limma")

BiocManager::install("affy")
library(affy)
DOC
browseVignettes("limma")
limmaUsersGuide()

#Load Data
data <- fread("https://raw.githubusercontent.com/ShyniProject/Rendu/master/GSE129081_small.csv", header = T)
#head(data)
# 
# 
ma.plot(MA$A, MA$M,
               xlab = "log2(count)",
               ylab = "log2(foldchange)",
               main = "MA-Plot",
               status = status,
               values=values,
               hl.col=col, cex =1,

               )
# #Def Variable
MA = new("MAList")
MA$A= log(data$baseMean,base = 2)
MA$M= data$log2FoldChange
status = rep("Gene",nrow(data))
threshold.padj=0.01
status[which(!is.na(data$padj))]="NA"
status[which(data$padj>threshold.padj)]="A"
values <- c("A","NA", "Gene")
col <- c("blue","black", "red")
#Plot
plotMA(object = MA,
       xlab = "log2(count)",
       ylab = "log2(foldchange)",
       main = "MA-Plot",
       status = status,
       values=values,
       hl.col=col
       #       array = MA$A,
       #       coef = MA$M,
       #       status = ,
       #       zero.weights = ,
       )

library(ggplot2)
ggplot(data,aes(x=log(baseMean,base = 2), y=log2FoldChange)) +
  geom_point() +
  coord_cartesian() +
  ylab("log2 FC") +
  xlab("log2 CPM")
