###### Volcano Plot

#library
library(data.table)
library(calibrate)

#Load Data
data <- fread("https://raw.githubusercontent.com/ShyniProject/Rendu/master/GSE129081_small.csv", header = T)
#head(data)

#Def Variable
threshold.padj=0.05
threshold.log2=1
x.limit=c(-2.5,2)
x.legend.pos=-1.25
y.legend.pos=100
col1="red"
col2="blue"
col3="green"

# Make a basic volcano plot
with(data, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=x.limit)) #pch donne la forme des points sur le graph. On pourrait proposer ça en option

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(data, padj<threshold.padj ), points(log2FoldChange, -log10(pvalue), pch=20, col=col1))
with(subset(data, abs(log2FoldChange)>threshold.log2), points(log2FoldChange, -log10(pvalue), pch=20, col=col2))
with(subset(data, padj<threshold.padj & abs(log2FoldChange)>threshold.log2), points(log2FoldChange, -log10(pvalue), pch=20, col=col3))

# Label points with the textxy function from the calibrate plot
with(subset(data, padj<threshold.padj & abs(log2FoldChange)>threshold.log2),textxy(log2FoldChange, -log10(pvalue), labs=Orthologous_human_gene, cex=.8))
legend(x.legend.pos, y.legend.pos,
       legend=c(paste("padj<",threshold.padj,sep = ""), 
                paste("|log2FoldChange|>",threshold.log2,sep = ""),
                paste("padj<",threshold.padj," &"," |log2FoldChange|>",threshold.log2, sep = "")),
       col=c(col1, col2, col3), lty=1:2, cex=0.8, pch=20)

