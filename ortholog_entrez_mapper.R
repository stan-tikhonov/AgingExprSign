library(biomaRt)
library(tidyverse)

# create entrez mappings between mouse and rat, and mouse and human

# default:
ensembl <- useMart("ensembl")
# in case www.ensembl.org is under maintenance:
ensembl = useEnsembl("ensembl", host = "uswest.ensembl.org")

datasets <- listDatasets(ensembl)
rat_dataset = useDataset("rnorvegicus_gene_ensembl", mart=ensembl)
human_dataset = useDataset("hsapiens_gene_ensembl", mart=ensembl)
mouse_dataset = useDataset("mmusculus_gene_ensembl", mart=ensembl)

# default:
rat = useMart("ensembl","rnorvegicus_gene_ensembl")
human = useMart("ensembl", "hsapiens_gene_ensembl")
mouse = useMart("ensembl","mmusculus_gene_ensembl")

# in case www.ensembl.org is under maintenance:
rat = useEnsembl("ensembl","rnorvegicus_gene_ensembl", host = "uswest.ensembl.org")
human = useEnsembl("ensembl", "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
mouse = useEnsembl("ensembl","mmusculus_gene_ensembl", host = "uswest.ensembl.org")

Rat_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                   mart=rat_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Rat_to_mouse_orthologs) <- c("Rat_Entrez","Mouse_Entrez")
Mouse_to_rat_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                   values=as.character(unique(Rat_to_mouse_orthologs$Mouse_Entrez)),
                                   mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=rat_dataset)
colnames(Mouse_to_rat_orthologs) <- c("Mouse_Entrez","Rat_Entrez")
Rat_to_mouse_orthologs1 = Rat_to_mouse_orthologs %>% group_by(Rat_Entrez) %>% filter(n() == 1)
Rat_to_mouse_orthologs1 = na.omit(Rat_to_mouse_orthologs1)
Mouse_to_rat_orthologs1 = subset(Mouse_to_rat_orthologs, Mouse_Entrez %in% Rat_to_mouse_orthologs1$Mouse_Entrez)
Mouse_to_rat_orthologs1 = Mouse_to_rat_orthologs1 %>% group_by(Mouse_Entrez) %>% filter(n() == 1)
mouse_rat_entrez_map = na.omit(Mouse_to_rat_orthologs1)
mouse_rat_entrez_map = as.data.frame(mouse_rat_entrez_map)
mouse_rat_entrez_map = mouse_rat_entrez_map %>% mutate_all(as.character)

Human_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                 mart=human_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Human_to_mouse_orthologs) <- c("Human_Entrez","Mouse_Entrez")
Mouse_to_human_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                 values=as.character(unique(Human_to_mouse_orthologs$Mouse_Entrez)),
                                 mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=human_dataset)
colnames(Mouse_to_human_orthologs) <- c("Mouse_Entrez","Human_Entrez")
Human_to_mouse_orthologs1 = Human_to_mouse_orthologs %>% group_by(Human_Entrez) %>% filter(n() == 1)
Human_to_mouse_orthologs1 = na.omit(Human_to_mouse_orthologs1)
Mouse_to_human_orthologs1 = subset(Mouse_to_human_orthologs, Mouse_Entrez %in% Human_to_mouse_orthologs1$Mouse_Entrez)
Mouse_to_human_orthologs1 = Mouse_to_human_orthologs1 %>% group_by(Mouse_Entrez) %>% filter(n() == 1)
mouse_human_entrez_map = na.omit(Mouse_to_human_orthologs1)
mouse_human_entrez_map = as.data.frame(mouse_human_entrez_map)
mouse_human_entrez_map = mouse_human_entrez_map %>% mutate_all(as.character)

#correlation matrix and heatmap
logFCmatrix = logFClist$Mouse$GSE6591$Lung$Male$DBA2J["logFC"]
logFCmatrix = logFCmatrix %>% rename(Mouse_GSE6581_Lung_Male_DBA2J = logFC)
logFCmatrix = merge(logFCmatrix, logFClist$Mouse$GSE6591$Lung$Male$C57BL6J["logFC"], by=0)
logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
logFCmatrix = logFCmatrix %>% rename(Mouse_GSE6581_Lung_Male_C57BL6J = logFC)
logFClist1 = logFClist
logFClist1$Mouse$GSE6591 = NULL
totalrownames = rownames(logFCmatrix)

#correct version :)
#j = 1
for (species in names(logFClist1)){
  for (dataset in names(logFClist1[[species]])){
    for (tissue in names(logFClist1[[species]][[dataset]])){
      for (sex in names(logFClist1[[species]][[dataset]][[tissue]])){
        if (species == "Rat"){
          exd <- logFClist1[[species]][[dataset]][[tissue]][[sex]]
          rownames(exd) = as.character(rownames(exd))
          exd$Rat_Entrez = rownames(exd)
          exd <- left_join(exd, mouse_rat_entrez_map)
          exd = na.omit(exd, cols=Mouse_Entrez)
          logFClist1[[species]][[dataset]][[tissue]][[sex]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
        }
        if (species == "Human"){
          exd <- logFClist1[[species]][[dataset]][[tissue]][[sex]]
          rownames(exd) = as.character(rownames(exd))
          exd$Human_Entrez = rownames(exd)
          exd <- left_join(exd, mouse_human_entrez_map)
          exd = na.omit(exd, cols=Mouse_Entrez)
          logFClist1[[species]][[dataset]][[tissue]][[sex]] <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
        }
        totalrownames = union(totalrownames, rownames(logFClist1[[species]][[dataset]][[tissue]][[sex]]))
#        logFCmatrix = merge(logFCmatrix, logFClist1[[species]][[dataset]][[tissue]][[sex]]["logFC"], by=0, all=TRUE)
#        logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
#        colnames(logFCmatrix)[j + 2] = paste0(species, "_", dataset, "_", tissue, "_", sex)
#        j = j + 1
      }
    }
  }
}

# replenish logFClist1:
logFClist1$Mouse$GSE6591 = logFClist$Mouse$GSE6591

# create denoised correlation matrix:
logFCunlisted = list()
logFCunlisted[["Mouse_GSE6591_Lung_Male_DBA2J"]] = logFClist1$Mouse$GSE6591$Lung$Male$DBA2J
logFCunlisted[["Mouse_GSE6591_Lung_Male_C57BL6J"]] = logFClist1$Mouse$GSE6591$Lung$Male$C57BL6J
logFClist1$Mouse$GSE6591 = NULL
for (species in names(logFClist1)){
  for (dataset in names(logFClist1[[species]])){
    for (tissue in names(logFClist1[[species]][[dataset]])){
      for (sex in names(logFClist1[[species]][[dataset]][[tissue]])){
        logFCunlisted[[paste(species, dataset, tissue, sex, sep = "_")]] = logFClist1[[species]][[dataset]][[tissue]][[sex]]
      }
    }
  }
}

# create logFCmatrix:
logFCmatrix = matrix(nrow = length(totalrownames), ncol = length(logFCunlisted))
rownames(logFCmatrix) = totalrownames
colnames(logFCmatrix) = names(logFCunlisted)
logFCmatrix[rownames(logFClist$Mouse$GSE6591$Lung$Male$DBA2J["logFC"]), "Mouse_GSE6591_Lung_Male_DBA2J"] = logFClist$Mouse$GSE6591$Lung$Male$DBA2J$logFC
logFCmatrix[rownames(logFClist$Mouse$GSE6591$Lung$Male$C57BL6J["logFC"]), "Mouse_GSE6591_Lung_Male_C57BL6J"] = logFClist$Mouse$GSE6591$Lung$Male$C57BL6J$logFC
logFClist1$Mouse$GSE6591 = NULL
for (species in names(logFClist1)){
  for (dataset in names(logFClist1[[species]])){
    for (tissue in names(logFClist1[[species]][[dataset]])){
      for (sex in names(logFClist1[[species]][[dataset]][[tissue]])){
        logFCmatrix[rownames(logFClist1[[species]][[dataset]][[tissue]][[sex]]), paste(species, dataset, tissue, sex, sep = "_")] = logFClist1[[species]][[dataset]][[tissue]][[sex]]$logFC
      }
    }
  }
}
logFCmatrix = as.data.frame(logFCmatrix)
logFCmatrixunfiltered = logFCmatrix
# filter genes with lots of NAs (half of the number of columns or more)
logFCmatrix$NACount = rowSums(is.na(logFCmatrix))
ggplot(logFCmatrix, aes(x = NACount)) + geom_density()
# choose the threshold according to the plot's bimodal distribution:
logFCmatrix = logFCmatrix %>% rownames_to_column(var = "row.names")
logFCmatrix = logFCmatrix %>% filter(NACount < 35) # the treshold is set according to the plot, 15k genes that are left is OK
logFCmatrix = logFCmatrix %>% column_to_rownames(var = "row.names")
logFCmatrix$NACount = NULL

# create normal cormatrix
cormatrixnafiltered <- round(cor(logFCmatrix, method = "spearman", use = "complete.obs"),2)


# this is for choosing one species:
j = 1
for (species in names(logFClist1)){
  for (dataset in names(logFClist1[[species]])){
    for (tissue in names(logFClist1[[species]][[dataset]])){
      for (sex in names(logFClist1[[species]][[dataset]][[tissue]])){
        logFCmatrix = merge(logFCmatrix, logFClist1[[species]][[dataset]][[tissue]][[sex]]["logFC"], by=0)
        logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
        colnames(logFCmatrix)[j + 1] = paste0(species, "_", dataset, "_", tissue, "_", sex)
        j = j + 1
      }
    }
  }
}

# replace NAs by zeros:
#logFCmatrix = logFCmatrix %>% mutate_all(~replace(., is.na(.), 0))

cormatrixsign = list()
cormatrixsign[["100"]] = data.frame()
cormatrixsign[["200"]] = data.frame()
cormatrixsign[["300"]] = data.frame()
cormatrixsign[["400"]] = data.frame()
cormatrixsign[["500"]] = data.frame()
cormatrixsign[["750"]] = data.frame()
cormatrixsign[["1000"]] = data.frame()
cormatrixsign[["1500"]] = data.frame()
cormatrixsign[["2000"]] = data.frame()
for (thres in names(cormatrixsign)){
  for (i in 1:length(logFCunlisted)){
    for (j in i:length(logFCunlisted)){
      topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
      topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topA = topA %>% column_to_rownames(var = "row.names")
      topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
      topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topB = topB %>% column_to_rownames(var = "row.names")
      totalrownames = union(rownames(topA), rownames(topB))
      cormatrixsign[[thres]][names(logFCunlisted)[i], names(logFCunlisted)[j]] = cor(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = "spearman", use = "complete.obs")
      cormatrixsign[[thres]][names(logFCunlisted)[j], names(logFCunlisted)[i]] = cor(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = "spearman", use = "complete.obs")
  #    mergedmatrix = logFCunlisted[[i]]["logFC"]
  #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
  #    mergedmatrix = merge(mergedmatrix, logFCunlisted[[j]]["logFC"], by=0, all=TRUE)
  #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
  #    cormatrixdenoised[names(logFCunlisted)[i], names(logFCunlisted)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
    }
  }
}
cormatrixsign[["all"]] = as.data.frame(cor(logFCmatrixunfiltered, method = "spearman", use = "complete.obs"))

# performing cor.test to find statisticaly significant correlations
corpvalsign = list()
corpvalsign[["100"]] = data.frame()
corpvalsign[["200"]] = data.frame()
corpvalsign[["300"]] = data.frame()
corpvalsign[["400"]] = data.frame()
corpvalsign[["500"]] = data.frame()
corpvalsign[["750"]] = data.frame()
corpvalsign[["1000"]] = data.frame()
corpvalsign[["1500"]] = data.frame()
corpvalsign[["2000"]] = data.frame()
for (thres in names(corpvalsign)){
  for (i in 1:length(logFCunlisted)){
    for (j in i:length(logFCunlisted)){
      topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
      topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topA = topA %>% column_to_rownames(var = "row.names")
      topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
      topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topB = topB %>% column_to_rownames(var = "row.names")
      totalrownames = union(rownames(topA), rownames(topB))
      corpvalsign[[thres]][names(logFCunlisted)[i], names(logFCunlisted)[j]] = as.numeric(cor.test(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = "spearman")$p.value)
      corpvalsign[[thres]][names(logFCunlisted)[j], names(logFCunlisted)[i]] = as.numeric(cor.test(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = "spearman")$p.value)
      #    mergedmatrix = logFCunlisted[[i]]["logFC"]
      #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
      #    mergedmatrix = merge(mergedmatrix, logFCunlisted[[j]]["logFC"], by=0, all=TRUE)
      #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
      #    cormatrixdenoised[names(logFCunlisted)[i], names(logFCunlisted)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
    }
  }
}

corpvalsign[["all"]] = data.frame()
for (colnum1 in 1:length(colnames(logFCmatrixunfiltered))){
  for (colnum2 in colnum1:length(colnames(logFCmatrixunfiltered))){
    corpvalsign[["all"]][colnames(logFCmatrixunfiltered)[colnum1], colnames(logFCmatrixunfiltered)[colnum2]] = as.numeric(cor.test(logFCmatrixunfiltered[, colnum1], logFCmatrixunfiltered[, colnum2], method = "spearman")$p.value)
    corpvalsign[["all"]][colnames(logFCmatrixunfiltered)[colnum2], colnames(logFCmatrixunfiltered)[colnum1]] = as.numeric(cor.test(logFCmatrixunfiltered[, colnum1], logFCmatrixunfiltered[, colnum2], method = "spearman")$p.value)
  }
}

coradjpvalsign = list()
coradjpvalsign[["100"]] = data.frame()
coradjpvalsign[["200"]] = data.frame()
coradjpvalsign[["300"]] = data.frame()
coradjpvalsign[["400"]] = data.frame()
coradjpvalsign[["500"]] = data.frame()
coradjpvalsign[["750"]] = data.frame()
coradjpvalsign[["1000"]] = data.frame()
coradjpvalsign[["1500"]] = data.frame()
coradjpvalsign[["2000"]] = data.frame()
coradjpvalsign[["all"]] = data.frame()
for (thres in names(corpvalsign)){
  vec = as.vector(corpvalsign[[thres]][upper.tri(corpvalsign[[thres]], diag = F)])
  vec = p.adjust(vec, method = "BH")
  tempmatrix = matrix(0, 63, 63)
  tempmatrix[lower.tri(tempmatrix, diag = F)] = vec
  tempmatrix = t(tempmatrix)
  tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
  coradjpvalsign[[thres]] = tempmatrix
  colnames(coradjpvalsign[[thres]]) = colnames(corpvalsign[[thres]])
  rownames(coradjpvalsign[[thres]]) = rownames(corpvalsign[[thres]])
}

library(reshape2)
reorder_cormat <- function(cormat,method="complete"){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd,method = method)
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
cormatrix_2 <- reorder_cormat(cormatrixnafiltered,method="average")
cormatrix_2 = apply(cormatrix_2, 2, rev)
upper_tri <- get_upper_tri(cormatrix_2)
# Melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat <- melt(cormatrix_2, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                       midpoint = 0, limit = c(-0.4,0.4), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmap

# Reorder the correlation matrix
cormatrix_2 <- reorder_cormat(cormatrixsign[["200"]],method="average")
cormatrix_2 = apply(cormatrix_2, 2, rev)
upper_tri <- get_upper_tri(cormatrix_2)
# Melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat <- melt(cormatrix_2, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                       midpoint = 0, limit = c(-0.8,0.8), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmap

# analyze human:
logFCmatrixhuman = dplyr::select(logFCmartixunfiltered, contains("Human"))
# filter genes with lots of NAs (half of the number of columns or more)
logFCmatrixhuman$NACount = rowSums(is.na(logFCmatrixhuman))
ggplot(logFCmatrixhuman, aes(x = NACount)) + geom_density()
# choose the threshold according to the plot's bimodal distribution:
logFCmatrixhuman = logFCmatrixhuman %>% rownames_to_column(var = "row.names")
logFCmatrixhuman = logFCmatrixhuman %>% filter(NACount < 2) # the treshold is set according to the plot, 15k genes that are left is OK
logFCmatrixhuman = logFCmatrixhuman %>% column_to_rownames(var = "row.names")
logFCmatrixhuman$NACount = NULL

cormatrixhuman <- round(cor(logFCmatrixhuman, method = "spearman", use = "complete.obs"),2)

# Reorder the correlation matrix
cormatrix_2 <- reorder_cormat(cormatrixhuman,method="average")
cormatrix_2 = apply(cormatrix_2, 2, rev)
upper_tri <- get_upper_tri(cormatrix_2)
# Melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat <- melt(cormatrix_2, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                       midpoint = 0, limit = c(-0.8,0.8), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmap