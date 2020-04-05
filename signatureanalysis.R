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

# for all signatures:

totalgenes = rownames(agingsignatures_v3[["Human"]])
for (name in names(agingsignatures_v3)){
  totalgenes = union(totalgenes, rownames(agingsignatures_v3[[name]]))
}

signaturematrix = matrix(nrow = length(totalgenes), ncol = length(names(agingsignatures_v3)))

colnames(signaturematrix) = names(agingsignatures_v3)
rownames(signaturematrix) = totalgenes

for (name in colnames(signaturematrix)){
  signaturematrix[rownames(agingsignatures_v3[[name]]), name] = agingsignatures_v3[[name]]$logFC
}
cormatrixall <- round(cor(signaturematrix, method = "spearman", use = "pairwise.complete.obs"),2)


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
cormatrix <- reorder_cormat(cormatrixall) # here, type in cormatrixtissues or cormatrixspecies
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
library(ggdendro)
library(grid)

totalgenes = c()
for (name in names(agingsignatures_v3)){
  totalgenes = union(totalgenes, rownames(agingsignatures_v3[[name]]))
}

logFCforclus = matrix(ncol = length(agingsignatures_v3), nrow = length(totalgenes))
rownames(logFCforclus) = totalgenes
colnames(logFCforclus) = names(agingsignatures_v3)
pvalsforclus = matrix(ncol = length(agingsignatures_v3), nrow = length(totalgenes))
rownames(pvalsforclus) = totalgenes
colnames(pvalsforclus) = names(agingsignatures_v3)
for (name in names(agingsignatures_v3)){
  logFCforclus[rownames(agingsignatures_v3[[name]]), name] = agingsignatures_v3[[name]]$logFC
  pvalsforclus[rownames(agingsignatures_v3[[name]]), name] = agingsignatures_v3[[name]]$adj_pval
}

load("deminglistforsignatures.RData")
# get the deming coefs:
minimums = c()
for (i in 1:10){
  #deminglistforsignatures[[i]] = deming_minimizer(logFCmatrixchosen)
  minimums = c(minimums, deminglistforsignatures[[i]]$minimum)
}
kres = deminglistforsignatures[[which.min(minimums)]]$coefs

for (i in 1:length(colnames(logFCforclus))){
  logFCforclus[,i] = logFCforclus[,i] / kres[i]
}

logFCforclus = na.omit(logFCforclus)
logFCforclus = as.data.frame(logFCforclus)
pvalsforclus = na.omit(pvalsforclus)
pvalsforclus = as.data.frame(pvalsforclus)
# filter genes by significance:
# geometric mean of p-values:
pvalsforclus$geommin = (pvalsforclus$Human * pvalsforclus$Rat * pvalsforclus$Mouse * pvalsforclus$Brain * pvalsforclus$Muscle * pvalsforclus$Liver * pvalsforclus$All) ^ (1/7)
logFCforclusfiltered = logFCforclus[rownames(pvalsforclus %>% rownames_to_column("Row.names") %>% filter(geommin < 0.2) %>% column_to_rownames("Row.names")),]

logFCforclusfiltered = as.matrix(logFCforclusfiltered)
genedendro = as.dendrogram(hclust(d = dist(x = logFCforclusfiltered, method = "euclidean"), method = "complete"))
dendro.plot <- ggdendrogram(data = genedendro, rotate = TRUE)
print(dendro.plot)

dendroorder <- order.dendrogram(genedendro)

tempshit = as.data.frame(logFCforclusfiltered) %>% rownames_to_column(var = "id")
meltedshit = gather(tempshit, dataset, logFC, -id, factor_key = T)
meltedshit$id = factor(meltedshit$id, levels = tempshit$id[dendroorder], ordered = T)
heatmap.plot <- ggplot(meltedshit, aes(dataset, id, fill = logFC))+
  geom_tile()+
  scale_fill_gradient2(trans = scales::pseudo_log_trans(base = exp(0.1), sigma = 0.1), low = "blue4", high = "red4", mid = "white", 
                       midpoint = 0, space = "Lab", 
                       name="LogFC")
print(heatmap.plot)

grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.5, width = 0.2, height = 1.087))

#dd = dist(logFCforclus, method = "manhattan")
#hc <- hclust(dd,method = "average")
#clusteredshit <-logFCforclus[hc$order,]
#clusteredshit = clusteredshit %>% rownames_to_column(var = "id")
#meltedshit = gather(clusteredshit, dataset, logFC, -id, factor_key = T)
#meltedshit$id = factor(meltedshit$id, levels = as.character(meltedshit$id[1:length(rownames(logFCforclus))]))
#heatmap.plot <- ggplot(meltedshit, aes(dataset, id, fill = logFC))+
#  geom_tile()+
#  scale_fill_gradient2(trans = scales::pseudo_log_trans(base = exp(0.1), sigma = 0.1), low = "blue4", high = "red4", mid = "white", 
#                       midpoint = 0, space = "Lab", 
#                       name="LogFC")
#ggheatmap

# gsea heatmap:
totalposrownames = c()
totalnegrownames = c()
for (name in names(agingsignatures_v3)){
  pos = read.table(paste0("./GSEA_result/", name, "/", name, "_positive.xls"), sep = "\t", header = T)
  neg = read.table(paste0("./GSEA_result/", name, "/", name, "_negative.xls"), sep = "\t", header = T)
  pos = pos %>% filter(FDR.q.val < 0.05)
  totalposrownames = c(totalposrownames, as.character(pos$NAME))
  neg = neg %>% filter(FDR.q.val < 0.05)
  totalnegrownames = c(totalnegrownames, as.character(neg$NAME))
}
postable = matrix(nrow = length(totalposrownames), ncol = 7)
rownames(postable) = totalposrownames
colnames(postable) = names(agingsignatures_v3)
negtable = matrix(nrow = length(totalnegrownames), ncol = 7)
rownames(negtable) = totalnegrownames
colnames(negtable) = names(agingsignatures_v3)
for (name in names(agingsignatures_v3)){
  pos = read.table(paste0("./GSEA_result/", name, "/", name, "_positive.xls"), sep = "\t", header = T)
  neg = read.table(paste0("./GSEA_result/", name, "/", name, "_negative.xls"), sep = "\t", header = T)
  pos = pos %>% filter(FDR.q.val < 0.05)
  neg = neg %>% filter(FDR.q.val < 0.05)
  postable[as.character(pos$NAME), name] = pos$NES
  negtable[as.character(neg$NAME), name] = neg$NES
}
totaltable = rbind(na.omit(postable), na.omit(negtable))
totaltable = as.data.frame(totaltable)
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
totaltable = totaltable %>% rownames_to_column("Row.names")
melted_cormat <- gather(totaltable, key = "key", value = "value", -Row.names)
melted_cormat$Row.names = factor(melted_cormat$Row.names, levels = melted_cormat$Row.names[1:41])
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(key, Row.names, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, space = "Lab") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)


tissuetable = rbind(na.omit(postable[,c("Brain", "Muscle", "Liver")]), na.omit(negtable[,c("Brain", "Muscle", "Liver")]))
tissuetable = as.data.frame(tissuetable)
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
tissuetable = tissuetable %>% rownames_to_column("Row.names")
melted_cormat <- gather(tissuetable, key = "key", value = "value", -Row.names)
melted_cormat$Row.names = factor(melted_cormat$Row.names, levels = melted_cormat$Row.names[1:70])
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(key, Row.names, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, space = "Lab") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)

speciestable = rbind(na.omit(postable[,c("Human", "Mouse", "Rat")]), na.omit(negtable[,c("Human", "Mouse", "Rat")]))
speciestable = as.data.frame(speciestable)
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
speciestable = speciestable %>% rownames_to_column("Row.names")
melted_cormat <- gather(speciestable, key = "key", value = "value", -Row.names)
melted_cormat$Row.names = factor(melted_cormat$Row.names, levels = melted_cormat$Row.names[1:73])
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(key, Row.names, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, space = "Lab") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)


##### deming regression for signatures

# make cortestsign matrix:

adjmethod = "BH"
cormethod = "pearson"
thres = "750"
corthres = 0.1
corpvalsign = data.frame()
cormatrixsign = data.frame()
for (i in 1:length(agingsignatures_v3)){
  for (j in i:length(agingsignatures_v3)){
    topA = agingsignatures_v3[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj_pval)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = agingsignatures_v3[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj_pval)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    cormatrixsign[names(agingsignatures_v3)[i], names(agingsignatures_v3)[j]] = cor(agingsignatures_v3[[i]][totalrownames,]$logFC, agingsignatures_v3[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    cormatrixsign[names(agingsignatures_v3)[j], names(agingsignatures_v3)[i]] = cor(agingsignatures_v3[[i]][totalrownames,]$logFC, agingsignatures_v3[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    corpvalsign[names(agingsignatures_v3)[i], names(agingsignatures_v3)[j]] = as.numeric(cor.test(agingsignatures_v3[[i]][totalrownames,]$logFC, agingsignatures_v3[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    corpvalsign[names(agingsignatures_v3)[j], names(agingsignatures_v3)[i]] = as.numeric(cor.test(agingsignatures_v3[[i]][totalrownames,]$logFC, agingsignatures_v3[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    #    mergedmatrix = agingsignatures_v3[[i]]["logFC"]
    #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
    #    mergedmatrix = merge(mergedmatrix, agingsignatures_v3[[j]]["logFC"], by=0, all=TRUE)
    #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
    #    cormatrixdenoised[names(agingsignatures_v3)[i], names(agingsignatures_v3)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
  }
}
coradjpvalsign = data.frame()
vec = as.vector(corpvalsign[upper.tri(corpvalsign, diag = F)])
vec = p.adjust(vec, method = adjmethod)
tempmatrix = matrix(0, length(agingsignatures_v3), length(agingsignatures_v3))
tempmatrix[upper.tri(tempmatrix, diag = F)] = vec
tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
#tempmatrix[lower.tri(tempmatrix, diag = F)] = vec
#tempmatrix = t(tempmatrix)
#tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
coradjpvalsign = tempmatrix
colnames(coradjpvalsign) = colnames(corpvalsign)
rownames(coradjpvalsign) = rownames(corpvalsign)
cortestsign = data.frame()
for (colname in colnames(cormatrixsign)){
  for (rowname in rownames(cormatrixsign)){
    if (coradjpvalsign[rowname, colname] < 0.05){
      if (cormatrixsign[rowname, colname] > corthres){
        cortestsign[rowname, colname] = 1
      } else if (cormatrixsign[rowname, colname] < -corthres){
        cortestsign[rowname, colname] = -1
      } else {
        cortestsign[rowname, colname] = 0
      }
    } else{
      cortestsign[rowname, colname] = 0
    }
  }
}

# make totalrownamematrix:
totalrownamematrix = list()
for (i in 1:(length(agingsignatures_v3)-1)){
  for (j in (i + 1):length(agingsignatures_v3)){
    topA = agingsignatures_v3[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj_pval)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = agingsignatures_v3[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj_pval)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    #tempdata = matrix(nrow = length(totalrownames), ncol = 2)
    #rownames(tempdata) = totalrownames
    #tempdata[, 1] = agingsignatures_v3[[i]][totalrownames,]$logFC
    #tempdata[, 2] = agingsignatures_v3[[j]][totalrownames,]$logFC
    #tempdata = na.omit(tempdata)
    #totalrownames = rownames(tempdata)
    rownamesA = rownames(subset(agingsignatures_v3[[i]], rownames(agingsignatures_v3[[i]]) %in% totalrownames))
    rownamesB = rownames(subset(agingsignatures_v3[[j]], rownames(agingsignatures_v3[[j]]) %in% totalrownames))
    totalrownamematrix[[names(agingsignatures_v3)[i]]][[names(agingsignatures_v3)[j]]] = intersect(rownamesA, rownamesB)
    totalrownamematrix[[names(agingsignatures_v3)[j]]][[names(agingsignatures_v3)[i]]] = intersect(rownamesA, rownamesB)
  }
}



deminglistforsignatures = list()
source("FUN.Deming_minimizer.R")
minimums = c()
for (i in 1:10){
  deminglistforsignatures[[i]] = deming_minimizer(logFCforclus, totalrownamematrix)
  minimums = c(minimums, deminglistforsignatures[[i]]$minimum)
  print(paste0("Done with ", as.character(i), "th run"))
}
  
save(deminglistforsignatures, file = "deminglistforsignatures.RData")

# visualisation of coef distributions:

# here a is coefs and b is minimums
a = as.data.frame(deminglistforsignatures[[1]]$coefs)
for (i in 2:length(deminglistforsignatures)){
  a = cbind(a, deminglistforsignatures[[i]]$coefs)
}
colnames(a) = 1:10
a$group <- row.names(a)
a.m <- melt(a, id.vars = "group")
b = data.frame()
for (i in 1:10){
  b = rbind(b, deminglistforsignatures[[i]]$minimum)
}
b = cbind(b, 1:10)
colnames(b) = c("minimum", "variable")
temp = cortestsign
normcoef = sum(temp[upper.tri(temp, diag = F)] == 1)
b$minimum = b$minimum / normcoef
b$variable = as.factor(b$variable)
a.m$minimum = left_join(a.m, b, by = "variable")
print(ggplot(a.m$minimum, aes(group, value)) + geom_boxplot() + geom_jitter(aes(color = as.factor(round(minimum, 4)))))


##### PCA of datasets

# first obtain logFCmatrixregr via signature_builder.R

logFCmatrixregr = logFCmatrixregr[,-which(colnames(logFCmatrixregr) %in% c("Mouse_GSE11291_Heart_Male", "Mouse_GSE11291_Brain_Male", "Mouse_GSE11291_Muscle_Male"))]

# normalize:
for (i in 1:length(colnames(logFCmatrixregr))){
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}
load("deminglist.RData")

# get the deming coefs:
minimums = c()
for (i in 1:10){
  #deminglistforsignatures[[i]] = deming_minimizer(logFCmatrixchosen)
  minimums = c(minimums, deminglist[["All"]][[i]]$minimum)
}
kres = deminglist[["All"]][[which.min(minimums)]]$coefs

for (i in 1:length(colnames(logFCmatrixregr))){
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / kres[i]
}

logFCmatrixregr[is.na(logFCmatrixregr)] = 0

exprforpca = data.frame(t(scale(t(logFCmatrixregr))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
species = sub("([^_]+)_.*", "\\1", colnames(logFCmatrixregr))
tissues = sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\3", colnames(logFCmatrixregr))
cluster_plot + geom_point(aes(color = tissues, shape = species))
# outliers:
cluster_plot + geom_point(aes(color = tissues, shape = species)) + geom_text(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# MDS of datasets

source("FUN.Cormatricesmaker.R")
cormethod = "spearman"
thres = "750"
res = cormatricesmaker(logFCunlisted, cormethod, signifgenesthres = thres)
cortestsign = res$cortestsign
cormatrixsign = res$cormatrixsign

mdsfordatasets <- cmdscale((1 - cor(logFCmatrixregr, use = "complete.obs", method = "spearman"))/2)
mdsfordatasets = cbind(mdsfordatasets, species, tissues)
mdsfordatasets = as.data.frame(mdsfordatasets)
mdsfordatasets$V1 = as.numeric(as.character(mdsfordatasets$V1))
mdsfordatasets$V2 = as.numeric(as.character(mdsfordatasets$V2))
plot(mdsfordatasets[,1],mdsfordatasets[,2], pch=16, cex=1.5, xlab="Coordinate 1", ylab="Coordinate 2")
ggplot(mdsfordatasets, aes(x = V1, y = V2, color = species)) + geom_point()

# tSNE of datasets
library(M3C)
tsne(logFCmatrixregr, labels = as.factor(tissues)) # or as.factor(species) for coloring by species








