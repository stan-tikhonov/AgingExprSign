# subgroup analysis

# get signatures:
load("agingsignatures_v3.RData")

library(tidyverse)
library(reshape2)
library(VennDiagram)
library(annotate)
library(org.Hs.eg.db)

##### subgroup analysis

# correlation heatmap

# species:
totalgenes = rownames(agingsignatures_v3[["Human"]])
for (name in c("Rat", "Mouse")){
  totalgenes = union(totalgenes, rownames(agingsignatures_v3[[name]]))
}

signaturematrix = matrix(nrow = length(totalgenes), ncol = 3)

colnames(signaturematrix) = c("Human", "Rat", "Mouse")
rownames(signaturematrix) = totalgenes

for (name in colnames(signaturematrix)){
  signaturematrix[rownames(agingsignatures_v3[[name]]), name] = agingsignatures_v3[[name]]$logFC
}
cormatrixspecies <- round(cor(signaturematrix, method = "spearman", use = "pairwise.complete.obs"),2)

# tissues:
totalgenes = rownames(agingsignatures_v3[["Liver"]])
for (name in c("Muscle", "Brain")){
  totalgenes = union(totalgenes, rownames(agingsignatures_v3[[name]]))
}

signaturematrix = matrix(nrow = length(totalgenes), ncol = 3)

colnames(signaturematrix) = c("Liver", "Muscle", "Brain")
rownames(signaturematrix) = totalgenes

for (name in colnames(signaturematrix)){
  signaturematrix[rownames(agingsignatures_v3[[name]]), name] = agingsignatures_v3[[name]]$logFC
}
cormatrixtissues <- round(cor(signaturematrix, method = "spearman", use = "pairwise.complete.obs"),2)

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Reorder the correlation matrix
cormatrix <- reorder_cormat(cormatrixtissues) # here, type in cormatrixtissues or cormatrixspecies
cormatrix = apply(cormatrix, 2, rev)
upper_tri <- get_upper_tri(cormatrix)
# Melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat <- melt(cormatrix, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.6,0.6), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)

# Venn diagrams

vennpvals = as.data.frame(t(c(1, 1, 1, 1)))
for (i in 1:(length(names(agingsignatures_v3))-1)){
  for (j in (i+1):length(names(agingsignatures_v3))){
    chisqtable = matrix(nrow = 2, ncol = 2)
    chisqtable[1,1] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval > 0.05, logFC < 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval > 0.05, logFC < 0))))
    chisqtable[1,2] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval > 0.05, logFC < 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval < 0.05, logFC < 0))))
    chisqtable[2,1] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval < 0.05, logFC < 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval > 0.05, logFC < 0))))
    chisqtable[2,2] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval < 0.05, logFC < 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval < 0.05, logFC < 0))))
    vennpvals = rbind(vennpvals, c(names(agingsignatures_v3)[i], names(agingsignatures_v3)[j], fisher.test(chisqtable)$p.value, "down"))
    
    chisqtable = matrix(nrow = 2, ncol = 2)
    chisqtable[1,1] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval > 0.05, logFC > 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval > 0.05, logFC > 0))))
    chisqtable[1,2] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval > 0.05, logFC > 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval < 0.05, logFC > 0))))
    chisqtable[2,1] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval < 0.05, logFC > 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval > 0.05, logFC > 0))))
    chisqtable[2,2] = length(intersect(rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[i]]], adj_pval < 0.05, logFC > 0)), rownames(subset(agingsignatures_v3[[names(agingsignatures_v3)[j]]], adj_pval < 0.05, logFC > 0))))
    vennpvals = rbind(vennpvals, c(names(agingsignatures_v3)[i], names(agingsignatures_v3)[j], fisher.test(chisqtable)$p.value, "up"))
  }
}

venn.diagram(
  x = list(rownames(subset(agingsignatures_v3[["Human"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Rat"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Mouse"]], logFC < 0 & adj_pval < 0.05))),
  category.names = c("Human" , "Rat", "Mouse"),
  filename = 'plots/signatureplots/venn_diagramm_species_down.png', imagetype = "png",  
  output=TRUE
)

venn.diagram(
  x = list(rownames(subset(agingsignatures_v3[["Human"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Rat"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Mouse"]], logFC > 0 & adj_pval < 0.05))),
  category.names = c("Human" , "Rat", "Mouse"),
  filename = 'plots/signatureplots/venn_diagramm_species_up.png', imagetype = "png",  
  output=TRUE
)
venn.diagram(
  x = list(rownames(subset(agingsignatures_v3[["Liver"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Muscle"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Brain"]], logFC < 0 & adj_pval < 0.05))),
  category.names = c("Liver" , "Muscle", "Brain"),
  filename = 'plots/signatureplots/venn_diagramm_tissues_down.png', imagetype = "png",  
  output=TRUE
)
venn.diagram(
  x = list(rownames(subset(agingsignatures_v3[["Liver"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Muscle"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures_v3[["Brain"]], logFC > 0 & adj_pval < 0.05))),
  category.names = c("Liver" , "Muscle", "Brain"),
  filename = 'plots/signatureplots/venn_diagramm_tissues_up.png', imagetype = "png",  
  output=TRUE
)

# preparing data for GSEA:

source("FUN.Entrez_converter.R")
agingsignatures_v3_for_gsea = list()
for (name in names(agingsignatures_v3)){
  agingsignatures_v3_for_gsea[[name]] = entrez_converter(agingsignatures_v3[[name]], from = "Mouse", to = "Human")
  symbols = getSYMBOL(rownames(agingsignatures_v3_for_gsea[[name]]), data = "org.Hs.eg.db")
  agingsignatures_v3_for_gsea[[name]]$genesymbol = symbols[rownames(agingsignatures_v3_for_gsea[[name]])]
}


for (name in names(agingsignatures_v3_for_gsea)){
  positive = agingsignatures_v3_for_gsea[[name]] %>% filter(logFC > 0) %>% arrange(adj_pval)
  negative = agingsignatures_v3_for_gsea[[name]] %>% filter(logFC < 0) %>% arrange(desc(adj_pval))
  positive$adj_pval = abs(log(positive$adj_pval))
  negative$adj_pval = abs(log(negative$adj_pval)) * -1
  tableforgsea = rbind(positive, negative)
  tableforgsea = tableforgsea[c("genesymbol","adj_pval")]
  rownames(tableforgsea) = NULL
  colnames(tableforgsea) = NULL
  write.table(tableforgsea, file = paste0("GSEA_", name, ".rnk"), row.names = F, quote = F, sep = "\t")
}

# gene clusterisation based on logFCmatrix:

totalgenes = c()
for (name in names(agingsignatures_v3)){
  totalgenes = union(totalgenes, rownames(agingsignatures_v3[[name]]))
}

logFCforclus = matrix(ncol = length(agingsignatures_v3), nrow = length(totalgenes))
rownames(logFCforclus) = totalgenes
colnames(logFCforclus) = names(agingsignatures_v3)
for (name in names(agingsignatures_v3)){
  logFCforclus[rownames(agingsignatures_v3[[name]]), name] = agingsignatures_v3[[name]]$logFC
}

logFCforclus[is.na(logFCforclus)] = 0
logFCforclus = as.data.frame(logFCforclus)
dd = dist(logFCforclus, method = "manhattan")
hc <- hclust(dd,method = "average")
clusteredshit <-logFCforclus[hc$order,]
clusteredshit = clusteredshit %>% rownames_to_column(var = "id")
meltedshit = gather(clusteredshit, dataset, logFC, -id, factor_key = T)
meltedshit$id = factor(meltedshit$id, levels = as.character(meltedshit$id[1:length(rownames(logFCforclus))]))
ggheatmap <- ggplot(meltedshit, aes(dataset, id, fill = logFC))+
  geom_tile()+
  scale_fill_gradient2(trans = scales::pseudo_log_trans(base = exp(0.1), sigma = 0.1), low = "blue4", high = "red4", mid = "white", 
                       midpoint = 0, space = "Lab", 
                       name="LogFC")
ggheatmap


