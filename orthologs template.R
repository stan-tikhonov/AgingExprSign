#Map yeast genes to mouse orthologs
#A) Find orthologs
library(biomaRt)
ensembl <- useMart("ensembl")
yeast_dataset = useDataset("scerevisiae_gene_ensembl", mart=ensembl)
mouse_dataset = useDataset("mmusculus_gene_ensembl", mart=ensembl)
yeast = useMart("ensembl","scerevisiae_gene_ensembl")
mouse = useMart("ensembl","mmusculus_gene_ensembl")
listAttributes(yeast_dataset)
Yeast_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                   values=union(as.character(signatures$Yeast$Deletion$Max$Entrez),as.character(signatures$Yeast$Species$Median$Entrez)),
                                   mart=yeast_dataset,attributesL=c("entrezgene_id","mgi_symbol"), martL=mouse_dataset)
colnames(Yeast_to_mouse_orthologs) <- c("Yeast_Entrez","Mouse_Entrez","Mouse_gene")
#write.csv(Yeast_to_mouse_orthologs,"Yeast_to_mouse_orthology.csv")
Mouse_to_yeast_orthologs <- getLDS(attributes=c("entrezgene_id","mgi_symbol"), filters="entrezgene_id",
                                   values=as.character(unique(Yeast_to_mouse_orthologs$Mouse_Entrez)),
                                   mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=yeast_dataset)
colnames(Mouse_to_yeast_orthologs) <- c("Mouse_Entrez","Mouse_gene","Yeast_Entrez")
#write.csv(Mouse_to_yeast_orthologs,"Mouse_to_yeast_orthology.csv")

#B) Leave only genes that uniquely map to mouse orthologs and which are uniquely mapped to by mouse orthologs
Yeast_mouse_unique <- Yeast_to_mouse_orthologs
dupl_Yeast_Entrez <- unique(Yeast_to_mouse_orthologs$Yeast_Entrez[which(duplicated(Yeast_to_mouse_orthologs$Yeast_Entrez))])
unique_Yeast_Entrez <- unique(Yeast_to_mouse_orthologs$Yeast_Entrez[!(Yeast_to_mouse_orthologs$Yeast_Entrez %in% dupl_Yeast_Entrez)])

dupl_Mouse_Entrez <- unique(Mouse_to_yeast_orthologs$Mouse_Entrez[which(duplicated(Mouse_to_yeast_orthologs$Mouse_Entrez))])
unique_Mouse_Entrez <- unique(Mouse_to_yeast_orthologs$Mouse_Entrez[!(Mouse_to_yeast_orthologs$Mouse_Entrez %in% dupl_Mouse_Entrez)])
length(unique_Mouse_Entrez)

Yeast_mouse_unique <- Yeast_to_mouse_orthologs[Yeast_to_mouse_orthologs$Yeast_Entrez %in% unique_Yeast_Entrez & Yeast_to_mouse_orthologs$Mouse_Entrez %in% unique_Mouse_Entrez,]
dim(Yeast_mouse_unique)
Mouse_Yeast_unique <- Mouse_to_yeast_orthologs[Mouse_to_yeast_orthologs$Mouse_Entrez %in% unique_Mouse_Entrez & Mouse_to_yeast_orthologs$Yeast_Entrez %in% unique_Yeast_Entrez,]
dim(Mouse_Yeast_unique)
rownames(Yeast_mouse_unique) <- Yeast_mouse_unique$Yeast_Entrez
#write.csv(Yeast_mouse_unique,"Yeast_mouse_orthologs_unique.csv")

#C) Add corresponding mouse ids
signatures$Yeast$Deletion$Max
for (a in names(signatures$Yeast)){
  for (b in names(signatures$Yeast[[a]])){
    print(a)
    print(b)
    print(nrow(signatures$Yeast[[a]][[b]]))
    signatures$Yeast[[a]][[b]] <- signatures$Yeast[[a]][[b]][!is.na(signatures$Yeast[[a]][[b]]$Entrez),]
    signatures$Yeast[[a]][[b]] <- signatures$Yeast[[a]][[b]][!duplicated(signatures$Yeast[[a]][[b]]$Entrez),]
    print(nrow(signatures$Yeast[[a]][[b]]))
    rownames(signatures$Yeast[[a]][[b]]) <- as.character(signatures$Yeast[[a]][[b]]$Entrez)
    signatures$Yeast[[a]][[b]]$Mouse_Entrez <- Yeast_mouse_unique[rownames(signatures$Yeast[[a]][[b]]),]$Mouse_Entrez
    signatures$Yeast[[a]][[b]]$Mouse_Gene <- Yeast_mouse_unique[rownames(signatures$Yeast[[a]][[b]]),]$Mouse_gene
    signatures$Yeast[[a]][[b]]$Zscore <- -log10(signatures$Yeast[[a]][[b]]$P.Value)*sign(signatures$Yeast[[a]][[b]]$logFC)
  }
}
head(signatures$Yeast$Deletion$Max)
setwd("D:/Documents/Сколково/PhD/analysis/Comparison of signatures/results/")
#dput(signatures,"Signatures_mouse.R")
