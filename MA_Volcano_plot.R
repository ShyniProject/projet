#######################################################
## Project: Fonction shiny project Volcano et MA plot #
## Script purpose: Tester les fonction de génération  #
## du volcano et du MA plot                           #
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
    dplyr::select(id:Orthologous_human_gene,Legendary:logBaseMean)
    return(mutate.input)
}

#VolcanoPlot
VolcanoPlot=function(data){
    ggplot(data,aes(x=log2FoldChange, y=negLogpadj, color=Legendary)) +
    geom_point() +
    coord_cartesian() +
    ylab("-log10 padj") +
    xlab("log2 fold change")
}
#MAPlot
MaPlot=function(data){
  ggplot(data,aes(x=logBaseMean, y=log2FoldChange, color=Legendary)) +
    geom_point() +
    coord_cartesian() +
    ylab("log2 FC") +
    xlab("log2 CPM")
}

