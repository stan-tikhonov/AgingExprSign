# subgroup analysis

# get signatures:
load("agingsignatures.RData")

library(tidyverse)
library(reshape2)
library(VennDiagram)
library(annotate)
library(org.Hs.eg.db)

##### subgroup analysis

# correlation heatmap

# species:
totalgenes = rownames(agingsignatures[["Human"]])
for (name in c("Rat", "Mouse")){
  totalgenes = union(totalgenes, rownames(agingsignatures[[name]]))
}

signaturematrix = matrix(nrow = length(totalgenes), ncol = 3)

colnames(signaturematrix) = c("Human", "Rat", "Mouse")
rownames(signaturematrix) = totalgenes

for (name in colnames(signaturematrix)){
  signaturematrix[rownames(agingsignatures[[name]]), name] = agingsignatures[[name]]$logFC
}
cormatrixspecies <- round(cor(signaturematrix, method = "spearman", use = "pairwise.complete.obs"),2)

# tissues:
totalgenes = rownames(agingsignatures[["Liver"]])
for (name in c("Muscle", "Brain")){
  totalgenes = union(totalgenes, rownames(agingsignatures[[name]]))
}

signaturematrix = matrix(nrow = length(totalgenes), ncol = 3)

colnames(signaturematrix) = c("Liver", "Muscle", "Brain")
rownames(signaturematrix) = totalgenes

for (name in colnames(signaturematrix)){
  signaturematrix[rownames(agingsignatures[[name]]), name] = agingsignatures[[name]]$logFC
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
cormatrix <- reorder_cormat(cormatrixspecies) # here, type in cormatrixtissues or cormatrixspecies
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
for (i in 1:(length(names(agingsignatures))-1)){
  for (j in (i+1):length(names(agingsignatures))){
    chisqtable = matrix(nrow = 2, ncol = 2)
    chisqtable[1,1] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval > 0.05, logFC < 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval > 0.05, logFC < 0))))
    chisqtable[1,2] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval > 0.05, logFC < 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval < 0.05, logFC < 0))))
    chisqtable[2,1] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval < 0.05, logFC < 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval > 0.05, logFC < 0))))
    chisqtable[2,2] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval < 0.05, logFC < 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval < 0.05, logFC < 0))))
    vennpvals = rbind(pvals, c(names(agingsignatures)[i], names(agingsignatures)[j], fisher.test(chisqtable)$p.value, "down"))
    
    chisqtable = matrix(nrow = 2, ncol = 2)
    chisqtable[1,1] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval > 0.05, logFC > 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval > 0.05, logFC > 0))))
    chisqtable[1,2] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval > 0.05, logFC > 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval < 0.05, logFC > 0))))
    chisqtable[2,1] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval < 0.05, logFC > 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval > 0.05, logFC > 0))))
    chisqtable[2,2] = length(intersect(rownames(subset(agingsignatures[[names(agingsignatures)[i]]], adj_pval < 0.05, logFC > 0)), rownames(subset(agingsignatures[[names(agingsignatures)[j]]], adj_pval < 0.05, logFC > 0))))
    vennpvals = rbind(pvals, c(names(agingsignatures)[i], names(agingsignatures)[j], fisher.test(chisqtable)$p.value, "up"))
  }
}

venn.diagram(
  x = list(rownames(subset(agingsignatures[["Human"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Rat"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Mouse"]], logFC < 0 & adj_pval < 0.05))),
  category.names = c("Human" , "Rat", "Mouse"),
  filename = 'plots/signatureplots/venn_diagramm_species_down.png', imagetype = "png",  
  output=TRUE
)

venn.diagram(
  x = list(rownames(subset(agingsignatures[["Human"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Rat"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Mouse"]], logFC > 0 & adj_pval < 0.05))),
  category.names = c("Human" , "Rat", "Mouse"),
  filename = 'plots/signatureplots/venn_diagramm_species_up.png', imagetype = "png",  
  output=TRUE
)
venn.diagram(
  x = list(rownames(subset(agingsignatures[["Liver"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Muscle"]], logFC < 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Brain"]], logFC < 0 & adj_pval < 0.05))),
  category.names = c("Liver" , "Muscle", "Brain"),
  filename = 'plots/signatureplots/venn_diagramm_tissues_down.png', imagetype = "png",  
  output=TRUE
)
venn.diagram(
  x = list(rownames(subset(agingsignatures[["Liver"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Muscle"]], logFC > 0 & adj_pval < 0.05)), rownames(subset(agingsignatures[["Brain"]], logFC > 0 & adj_pval < 0.05))),
  category.names = c("Liver" , "Muscle", "Brain"),
  filename = 'plots/signatureplots/venn_diagramm_tissues_up.png', imagetype = "png",  
  output=TRUE
)

# preparing data for GSEA:

source("FUN.Entrez_converter.R")
agingsignatures_for_gsea = list()
for (name in names(agingsignatures)){
  agingsignatures_for_gsea[[name]] = entrez_converter(agingsignatures[[name]], from = "Mouse", to = "Human")
  symbols = getSYMBOL(rownames(agingsignatures_for_gsea[[name]]), data = "org.Hs.eg.db")
  agingsignatures_for_gsea[[name]]$genesymbol = symbols[rownames(agingsignatures_for_gsea[[name]])]
}


for (name in names(agingsignatures_for_gsea)){
  positive = agingsignatures_for_gsea[[name]] %>% filter(logFC > 0) %>% arrange(adj_pval)
  negative = agingsignatures_for_gsea[[name]] %>% filter(logFC < 0) %>% arrange(desc(adj_pval))
  positive$adj_pval = abs(log(positive$adj_pval))
  negative$adj_pval = abs(log(negative$adj_pval)) * -1
  tableforgsea = rbind(positive, negative)
  tableforgsea = tableforgsea[c("genesymbol","adj_pval")]
  rownames(tableforgsea) = NULL
  colnames(tableforgsea) = NULL
  write.table(tableforgsea, file = paste0("GSEA_", name, ".gct"))
}





