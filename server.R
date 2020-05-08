function(input, output, session) 
{
  source("server_orgChoice.R")
  source("server_summary.R")
  source("server_kegg.R")
  source("server_go.R")
  source("server_proteinDomain.R")
  
  ##### ATTENTION pour l'instant l'upload du fichier lance les analyses, avant le choix de l'espèce donc.
  org <- organismChoice(input, output, session)
  
  observe({
    analysis()
  })
  
  analysis <- eventReactive(input$data,{
    orgDb <- eval(parse(text = org(), keep.source=FALSE))
    
    ########################################################
    ##########  data type check  ##########
    ########################################################
    formats = c(
      'text/csv',
      'text/comma-separated-values',
      'text/tab-separated-values',
      'text/plain',
      'csv',
      'tsv',
      "application/vnd.ms-excel"
    )
    
    if(!input$data$type %in% formats)
    {
      shinyalert("Oops!", paste0("It seems you have uploaded a ",input$data$type," file \n Please refresh page."), type = "warning", animation = T)
    } 
    
    ########################################################
    ##########  data preprocesses  ##########
    ########################################################
    ## RNAseq data for test ##
    # d <- read.csv("bims/GSE129081_small.csv")
    # entrezID <- mapIds(org.Dr.eg.db, as.vector(d$id), 'ENTREZID', 'ENSEMBL')
    
    d <- read.csv(input$data$datapath)
    entrezID <- mapIds(orgDb, as.vector(d$id), 'ENTREZID', 'ENSEMBL')
    
    ## for SEA : p-values vector and gene ids
    pvalues <- d$padj 
    names(pvalues) <- entrezID
    
    ## for GSEA : LogF vector and gene ids
    logF = d$log2FoldChange 
    names(logF) <- entrezID
    logF = logF[which(!is.na(names(logF)))] ; logF = na.exclude(logF)#17 809 samples
    ## delete duplicated genes -> median value computed
    if(length(which(duplicated(names(logF))==T))>0)
    {
      df = data.frame("id" = names(logF), "value" = logF)
      df %>%
        group_by(id) %>%
        summarise_all(median) %>%
        data.frame() -> newdf
      logF = newdf$value
      names(logF) = newdf$id # 17 777
    }
    
    logF = sort(logF, decreasing = TRUE) # decreasing vector to gsea
    
    ########################################################
    ##########  SUMMARY  ##########
    ########################################################
    summary(input, output, session, d)
    
    ########################################################
    ##########  KEGG  ##########
    ########################################################
    organismsDbKegg = c("org.Hs.eg.db"="hsa","org.Mm.eg.db"="mmu","org.Rn.eg.db"="rno",
                        "org.Sc.sgd.db"="sce","org.Dm.eg.db"="dme","org.At.tair.db"="ath",
                        "org.Dr.eg.db"="dre","org.Bt.eg.db"="bta","org.Ce.eg.db"="cel",
                        "org.Gg.eg.db"="gga","org.Cf.eg.db"="cfa","org.Ss.eg.db"="ssc",
                        "org.Mmu.eg.db"="mcc","org.EcK12.eg.db"="eck","org.Xl.eg.db"="xla",
                        "org.Pt.eg.db"="ptr","org.Ag.eg.db"="aga","org.Pf.plasmo.db"="pfa",
                        "org.EcSakai.eg.db"="ecs")
    
    kegg(input, output, session, org(), organismsDbKegg, pvalues, logF)

    ########################################################
    ##########  Protein Domains  ##########
    ########################################################
    # proteinDomain(input, output, session, org, organismsDbKegg, pvalues, logF, entrezID)
      
    ########################################################
    ##########                GO                  ##########
    ########################################################      
    # go(input, output, session, org, orgDb, pvalues, logF)
        
  } 
)}
