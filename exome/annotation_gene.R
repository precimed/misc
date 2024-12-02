options(connectionObserver = NULL)
library("AnnotationDbi")
library("org.Hs.eg.db")
hs <- org.Hs.eg.db
#columns(org.Hs.eg.db) # returns list of available keytypes
meanp1=fread('depression_gwas_v1.csv_meanP.genes.out.txt',data.table = F)
meanp1Entrez2Symbol <- select(hs, keys = as.character(meanp1$GENE),columns = c("ENTREZID", "SYMBOL","ENSEMBL",'GENENAME','ENZYME'),
                              keytype = "ENTREZID")
meanp1$GENE <- as.character(meanp1$GENE)
meanp1$SYMBOL=meanp1Entrez2Symbol$SYMBOL[match(meanp1$GENE,meanp1Entrez2Symbol$ENTREZID)]
meanp1$ENSEMBL=meanp1Entrez2Symbol$ENSEMBL[match(meanp1$GENE,meanp1Entrez2Symbol$ENTREZID)]
meanp1$GENENAME=meanp1Entrez2Symbol$GENENAME[match(meanp1$GENE,meanp1Entrez2Symbol$ENTREZID)]
#meanp2 <- meanp1[meanp1$P<.05/15073 & !duplicated(meanp1$SYMBOL) & !is.na(meanp1$SYMBOL),]
meanp2=meanp1[str_detect(meanp1$SYMBOL, 'KCNQ'),]
meanp3=meanp2[complete.cases(meanp2),]
fwrite(meanp3,file = 'kcn_related_genes_meanp.csv',sep = '\t',row.names = F,quote = F,na = 'NA')

map_function=function(genes){
  meanp1=fread(genes,data.table = F)
  hs <- org.Hs.eg.db
  meanp1Entrez2Symbol <- AnnotationDbi::select(hs, keys = as.character(meanp1$GENE),columns = c("ENTREZID", "SYMBOL","ENSEMBL",'GENENAME'),
                                keytype = "ENTREZID")
  meanp1$GENE <- as.character(meanp1$GENE)
  meanp1$SYMBOL=meanp1Entrez2Symbol$SYMBOL[match(meanp1$GENE,meanp1Entrez2Symbol$ENTREZID)]
  meanp1$ENSEMBL=meanp1Entrez2Symbol$ENSEMBL[match(meanp1$GENE,meanp1Entrez2Symbol$ENTREZID)]
  meanp1$GENENAME=meanp1Entrez2Symbol$GENENAME[match(meanp1$GENE,meanp1Entrez2Symbol$ENTREZID)]
  meanp1$P_adjusted=p.adjust(meanp1$P,method = 'fdr',n=nrow(meanp1))
  return(meanp1)
}

mgm1=map_function('depression_gwas_v2.csv_meanP.genes.out.txt')
mgm2=mgm1[str_detect(mgm1$SYMBOL, 'KCNQ'),]
mgm3=mgm2[complete.cases(mgm2),]

#Dementia inhouse
#Study 1: 1,2,3,4 as case and 17 as control
#Study 2: 1,2,3,4,5,6,10,11,12,13,14 as case and 17 as control

