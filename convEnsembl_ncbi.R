gene_ensembl = read.csv("GSE129081_small.csv", header = TRUE)
geneEnsembl <- as.data.frame(gene_ensembl)
ensembl.genes <- as.vector(geneEnsembl['id'])

if(interactive()){
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
}


mart <- useDataset("drerio_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=ensembl.genes,
  mart=mart)

keggId <- keggConv("dre", "ncbi-geneid") #danio rerio name "dre" for KEGGREST)




                        