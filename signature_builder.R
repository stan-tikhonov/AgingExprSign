# signature builder v. 2.0

library(biomaRt)
library(tidyverse)
library(deming)
library(reshape2)
library(metafor)

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


#remove bad datasets:
logFClist1 = logFClist
logFClist1$Mouse$GSE53959 = NULL
logFClist1$Human$GSE40645 = NULL
logFClist1$Human$GSE5086 = NULL

# convert to orthologs and make totalrownames:
totalrownames = union(rownames(logFClist$Mouse$GSE6591$Lung$Male$DBA2J), rownames(logFClist$Mouse$GSE6591$Lung$Male$C57BL6J))
logFClist1$Mouse$GSE6591 = NULL
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

# create logFCunlisted:
logFCunlisted = list()
logFCunlisted[["Mouse_GSE6591_Lung_Male_DBA2J"]] = logFClist$Mouse$GSE6591$Lung$Male$DBA2J
logFCunlisted[["Mouse_GSE6591_Lung_Male_C57BL6J"]] = logFClist$Mouse$GSE6591$Lung$Male$C57BL6J
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

SEmatrixregr = matrix(nrow = length(totalrownames), ncol = length(logFCunlisted))
rownames(SEmatrixregr) = totalrownames
colnames(SEmatrixregr) = names(logFCunlisted)

for (name in names(logFCunlisted)){
  logFCmatrix[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$logFC
  SEmatrixregr[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$SE
}

logFCmatrix = as.data.frame(logFCmatrix)
SEmatrixregr = as.data.frame(SEmatrixregr)
logFCmatrixregr = logFCmatrix

# get bad boys:
logFCmatrix$NACount = rowSums(is.na(logFCmatrix))
ggplot(logFCmatrix, aes(x = NACount)) + geom_density()

badboys = subset(rownames(logFCmatrix), logFCmatrix$NACount >=35)

# make cortestsign matrix:

cormethod = "pearson"
corpvalsign = data.frame()
cormatrixsign = data.frame()
thres = "750"
for (i in 1:length(logFCunlisted)){
  for (j in i:length(logFCunlisted)){
    topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    cormatrixsign[names(logFCunlisted)[i], names(logFCunlisted)[j]] = cor(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    cormatrixsign[names(logFCunlisted)[j], names(logFCunlisted)[i]] = cor(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
    corpvalsign[names(logFCunlisted)[i], names(logFCunlisted)[j]] = as.numeric(cor.test(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    corpvalsign[names(logFCunlisted)[j], names(logFCunlisted)[i]] = as.numeric(cor.test(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
    #    mergedmatrix = logFCunlisted[[i]]["logFC"]
    #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
    #    mergedmatrix = merge(mergedmatrix, logFCunlisted[[j]]["logFC"], by=0, all=TRUE)
    #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
    #    cormatrixdenoised[names(logFCunlisted)[i], names(logFCunlisted)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
  }
}

coradjpvalsign = data.frame()
vec = as.vector(corpvalsign[upper.tri(corpvalsign, diag = F)])
vec = p.adjust(vec, method = "BH")
tempmatrix = matrix(0, 58, 58)
tempmatrix[lower.tri(tempmatrix, diag = F)] = vec
tempmatrix = t(tempmatrix)
tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
coradjpvalsign = tempmatrix
colnames(coradjpvalsign) = colnames(corpvalsign)
rownames(coradjpvalsign) = rownames(corpvalsign)

cortestsign = data.frame()
for (colname in colnames(cormatrixsign)){
  for (rowname in rownames(cormatrixsign)){
    if (coradjpvalsign[rowname, colname] < 0.05){
      if (cormatrixsign[rowname, colname] > 0.1){
        cortestsign[rowname, colname] = 1
      } else if (cormatrixsign[rowname, colname] < -0.1){
        cortestsign[rowname, colname] = -1
      } else {
        cortestsign[rowname, colname] = 0
      }
    } else{
      cortestsign[rowname, colname] = 0
    }
  }
}

# normalize by sd:

for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}

# make totalrownamematrix:
totalrownamematrix = list()
for (i in 1:(length(logFCunlisted)-1)){
  for (j in (i + 1):length(logFCunlisted)){
    topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    #tempdata = matrix(nrow = length(totalrownames), ncol = 2)
    #rownames(tempdata) = totalrownames
    #tempdata[, 1] = logFCunlisted[[i]][totalrownames,]$logFC
    #tempdata[, 2] = logFCunlisted[[j]][totalrownames,]$logFC
    #tempdata = na.omit(tempdata)
    #totalrownames = rownames(tempdata)
    rownamesA = rownames(subset(logFCunlisted[[i]], rownames(logFCunlisted[[i]]) %in% totalrownames))
    rownamesB = rownames(subset(logFCunlisted[[j]], rownames(logFCunlisted[[j]]) %in% totalrownames))
    totalrownamematrix[[names(logFCunlisted)[i]]][[names(logFCunlisted)[j]]] = intersect(rownamesA, rownamesB)
  }
}

# get source table:
sourcedata = as.data.frame(colnames(logFCmatrixregr))
rownames(sourcedata) = colnames(logFCmatrixregr)
colnames(sourcedata) = "kekkekkek"
sourcedata = sourcedata %>% separate(kekkekkek, c(NA, "dataset", NA), sep = "_")


# run deming minimization:
source("FUN.Deming_minimizer.R")
bigres = list()
minimums = c()
for (i in 1:10){
  bigres[[i]] = deming_minimizer(logFCmatrixchosen)
  minimums = c(minimums, bigres[[i]]$minimum)
}
kres = bigres[[which.min(minimums)]]$coefs

# normalize by deming coefficients:

for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / kres[i]
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / kres[i]
}

# discard bad boys:
'%notin%' = Negate('%in%')
logFCmatrixregr = subset(logFCmatrixregr, rownames(logFCmatrixregr) %notin% badboys)
SEmatrixregr = subset(SEmatrixregr, rownames(SEmatrixregr) %notin% badboys)

# run mixed-effect model:

source("FUN.Signature_builder.R")
signature = signature_builder(logFCmatrixchosen)

# correlation heatmap

alexsignatures = dget("Signatures_mouse_genes.R")

signatureforall = subset(signature, adj_pval < 0.1)

signaturegenes = rownames(signatureforall)
for (el in alexsignatures$Species){
  signaturegenes = union(signaturegenes, rownames(subset(el, FDR < 0.1)))
}
for (el in alexsignatures$Interventions){
  signaturegenes = union(signaturegenes, rownames(subset(el, FDR < 0.1)))
}

signaturematrix = matrix(nrow = length(signaturegenes), ncol = 12)

colnames(signaturematrix) = c("agingall", names(alexsignatures$Species), names(alexsignatures$Interventions))
rownames(signaturematrix) = signaturegenes

signaturematrix[rownames(signatureforall), "agingall"] = signatureforall$logFC
for (name in names(alexsignatures$Species)){
  signaturematrix[rownames(subset(alexsignatures$Species[[name]], FDR < 0.1)), name] = subset(alexsignatures$Species[[name]], FDR < 0.1)$logFC
}
for (name in names(alexsignatures$Interventions)){
  signaturematrix[rownames(subset(alexsignatures$Interventions[[name]], FDR < 0.1)), name] = subset(alexsignatures$Interventions[[name]], FDR < 0.1)$logFC
}

cormatrix <- round(cor(signaturematrix, method = "spearman", use = "pairwise.complete.obs"),2)
cormatrix[is.na(cormatrix)] = 0

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
cormatrix <- reorder_cormat(cormatrix)
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
ggheatmap

signaturegenes = rownames(signature)
for (el in alexsignatures$Species){
  signaturegenes = union(signaturegenes, rownames(el))
}
for (el in alexsignatures$Interventions){
  signaturegenes = union(signaturegenes, rownames(el))
}

signaturematrix = matrix(nrow = length(signaturegenes), ncol = 12)

colnames(signaturematrix) = c("agingall", names(alexsignatures$Species), names(alexsignatures$Interventions))
rownames(signaturematrix) = signaturegenes

signaturematrix[rownames(signatureforall), "agingall"] = signatureforall$logFC
for (name in names(alexsignatures$Species)){
  signaturematrix[rownames(alexsignatures$Species[[name]]), name] = alexsignatures$Species[[name]]$logFC
}
for (name in names(alexsignatures$Interventions)){
  signaturematrix[rownames(alexsignatures$Interventions[[name]]), name] = alexsignatures$Interventions[[name]]$logFC
}

cormatrix <- round(cor(signaturematrix, method = "spearman", use = "pairwise.complete.obs"),2)

# Reorder the correlation matrix
cormatrix <- reorder_cormat(cormatrix)
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
ggheatmap
