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
#logFClist1$Mouse$`E-MTAB-3374` = NULL
#logFClist1$Human$GSE9103 = NULL

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
tempmatrix = matrix(0, length(logFCunlisted), length(logFCunlisted))
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

logFCmatrixchosen = logFCmatrixregr
SEmatrixchosen = SEmatrixregr

##### THIS IS FOR THE COMPLETE SIGNATURE (10 MINIMIZATION RUNS)

# run deming minimization:
source("FUN.Deming_minimizer.R")

# running it 10 times:
bigres = list()
minimums = c()
for (i in 1:10){
  bigres[[i]] = deming_minimizer(logFCmatrixchosen)
  minimums = c(minimums, bigres[[i]]$minimum)
  print(paste0("Opa! ", i, "th minimization done."))
}
kres = bigres[[which.min(minimums)]]$coefs

# normalize by deming coefficients:

for (i in 1:length(colnames(logFCmatrixchosen))){
  SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
  logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
}

# discard bad boys:
'%notin%' = Negate('%in%')
logFCmatrixshosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %notin% badboys)
SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %notin% badboys)

# run mixed-effect model:

source("FUN.Signature_builder.R")
signature = signature_builder(logFCmatrixchosen)



##### constructing category table:
categorytable = matrix(nrow = 6, ncol = 4)
rownames(categorytable) = c("Brain", "Liver", "Heart", "Kidney", "Muscle", "Lung")
colnames(categorytable) = c("Mouse", "Rat", "Human", "Total")
for (colnm in colnames(categorytable)[1:3]){
  for (rownm in rownames(categorytable)){
    if (rownm != "Brain"){
      logspecies = grepl(paste0(".*",colnm,".*"), colnames(logFCmatrixregr))
      logtissue = grepl(paste0(".*",rownm,".*"), colnames(logFCmatrixregr))
      categorytable[rownm, colnm] = sum(logspecies & logtissue)
    } else {
      logspecies = grepl(paste0(".*",colnm,".*"), colnames(logFCmatrixregr))
      logtissue = grepl(paste0(".*",rownm,".*"), colnames(logFCmatrixregr)) | grepl(paste0(".*","Cerebellum",".*"), colnames(logFCmatrixregr)) | grepl(paste0(".*","Frontalcortex",".*"), colnames(logFCmatrixregr))
      categorytable[rownm, colnm] = sum(logspecies & logtissue)
    }
  }
}
categorytable[,"Total"] = categorytable[,"Mouse"] + categorytable[,"Human"] + categorytable[,"Rat"]


##### THIS IS FOR MANY SIGNATURES BUT ONE MINIMIZATION RUN FOR EACH SIGNATURE

# functions for heatmaps:
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

# prep shit:
chosencols = list()
chosencols[["Human"]] = colnames(logFCmatrixregr)[grepl(".*Human.*", colnames(logFCmatrixregr))]
chosencols[["Rat"]] = colnames(logFCmatrixregr)[grepl(".*Rat.*", colnames(logFCmatrixregr))]
chosencols[["Mouse"]] = colnames(logFCmatrixregr)[grepl(".*Mouse.*", colnames(logFCmatrixregr))]
chosencols[["Brain"]] = colnames(logFCmatrixregr)[grepl(".*Brain.*", colnames(logFCmatrixregr)) | grepl(".*Frontalcortex.*", colnames(logFCmatrixregr)) | grepl(".*Cerebellum.*", colnames(logFCmatrixregr))]
chosencols[["Muscle"]] = colnames(logFCmatrixregr)[grepl(".*Muscle.*", colnames(logFCmatrixregr))]
chosencols[["Liver"]] = colnames(logFCmatrixregr)[grepl(".*Liver.*", colnames(logFCmatrixregr))]
chosencols[["All"]] = colnames(logFCmatrixregr)

source("FUN.Deming_minimizer.R")
source("FUN.Signature_builder.R")

agingsignatures = list()

for (name in names(chosencols)){
  # filter datasets for the individual signature:
  logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
  SEmatrixchosen = SEmatrixregr[, chosencols[[name]]]
  
  # plot correlation heatmap:
  cormatrix = data.frame()
  for (i in 1:length(colnames(logFCmatrixchosen))){
    cormatrix[colnames(logFCmatrixchosen)[i], colnames(logFCmatrixchosen)[i]] = 1
  }
  for (i in 1:(length(colnames(logFCmatrixchosen))-1)){
    for (j in (i+1):length(colnames(logFCmatrixchosen))){
      cormatrix[colnames(logFCmatrixchosen)[i], colnames(logFCmatrixchosen)[j]] = cor(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[i]]][[colnames(logFCmatrixchosen)[j]]],i], logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[i]]][[colnames(logFCmatrixchosen)[j]]],j], method = "spearman", use = "complete.obs")
      cormatrix[colnames(logFCmatrixchosen)[j], colnames(logFCmatrixchosen)[i]] = cormatrix[colnames(logFCmatrixchosen)[i], colnames(logFCmatrixchosen)[j]]
    }
  }
  cormatrix_2 <- reorder_cormat(cormatrix, method="average")
  cormatrix_2 = apply(cormatrix_2, 2, rev)
  upper_tri <- get_upper_tri(cormatrix_2)
  melted_cormat <- melt(cormatrix_2, na.rm = TRUE)
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  ggheatmap
  pdf(paste0("./plots/signatureplots/", name, "/initialheatmap", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # minimize once:
  kres = deming_minimizer(logFCmatrixchosen)$coefs
  # plot an example:
  plot1 = deming(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],1] ~ logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],2] - 1)
  plot(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],1], logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],2], main = paste0("Correlation: ", as.character(cormatrixsign[colnames(logFCmatrixchosen)[1], colnames(logFCmatrixchosen)[2]]), ", p-value: ", as.character(coradjpvalsign[colnames(logFCmatrixchosen)[1], colnames(logFCmatrixchosen)[2]])))
  abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
  abline(0, kres[2]/kres[1], col = "red", lwd = 2)
  
  ggheatmap = ggplot(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],], aes_string(x = colnames(logFCmatrixchosen)[1], y = colnames(logFCmatrixchosen)[2])) + geom_point() +
    geom_abline(slope = plot1$coefficients[2], intercept = 0, colour = "blue", size = 1) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_abline(slope = kres[2]/kres[1], intercept = 0, colour = "red", size = 1)
  ggheatmap
  pdf(paste0("./plots/signatureplots/", name, "/demingexample", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # normalize by deming coefficients:
  
  for (i in 1:length(colnames(logFCmatrixchosen))){
    SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
    logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
  }
  
  # discard bad boys:
  logFCmatrixchosen$NACount = rowSums(is.na(logFCmatrixchosen))
  if (name != "Liver"){
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = floor(length(colnames(logFCmatrixchosen))/2))
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrix$NACount < floor(length(colnames(logFCmatrixchosen))/2))
  } else {
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = 4)
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrix$NACount < 4)
  }
  logFCmatrixchosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %in% goodboys)
  SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %in% goodboys)
  logFCmatrixchosen$NACount = NULL
  
  # run mixed-effect model:
  agingsignatures[[name]] = signature_builder(logFCmatrixchosen)
  # plot an example:
  helpertable = as.data.frame(t(logFCmatrixchosen[1,]))
  rownames(helpertable) = colnames(logFCmatrixchosen)
  colnames(helpertable) = c("logFC")
  helpertable$SE = t(SEmatrixchosen[1,])
  helpertable$source = as.factor(sourcedata[rownames(helpertable),"dataset"])
  helpertable$dataset = rownames(helpertable)
  helpertable = na.omit(helpertable)
  ggheatmap = ggplot(helpertable, aes(x = dataset, y = logFC, color = source)) + geom_pointrange(aes(ymin = logFC - SE, ymax = logFC + SE)) + geom_hline(yintercept = agingsignatures[[name]][rownames(logFCmatrixchosen)[1], "logFC"], colour = "red") + geom_hline(yintercept = 0)
  ggheatmap
  pdf(paste0("./plots/signatureplots/", name, "/mixedmodelexample", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # plot verification cor heatmap:
  logFCmatrixchosen = merge(logFCmatrixchosen, agingsignatures[[name]]["logFC"], by = "row.names", all = TRUE)
  colnames(logFCmatrixchosen)[which(colnames(logFCmatrixchosen) == "logFC")] = paste0(name, "_signature")
  logFCmatrixchosen = logFCmatrixchosen %>% column_to_rownames(var = "Row.names")
  significantgenematrix = logFCmatrixchosen[rownames(subset(agingsignatures[[name]], adj_pval < 0.05)),]
  cormatrix = cor(significantgenematrix, use = "pairwise.complete.obs", method = "pearson")
  
  
  
  #cormatrix = data.frame()
  #for (i in 1:length(colnames(logFCmatrixchosen))){
  #  cormatrix[colnames(logFCmatrixchosen)[i], colnames(logFCmatrixchosen)[i]] = 1
  #}
  #for (i in 1:(length(colnames(logFCmatrixchosen))-1)){
  #  for (j in (i+1):length(colnames(logFCmatrixchosen))){
  #    cormatrix[colnames(logFCmatrixchosen)[i], colnames(logFCmatrixchosen)[j]] = cor(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[i]]][[colnames(logFCmatrixchosen)[j]]],i], logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[i]]][[colnames(logFCmatrixchosen)[j]]],j], method = "spearman", use = "complete.obs")
  #    cormatrix[colnames(logFCmatrixchosen)[j], colnames(logFCmatrixchosen)[i]] = cormatrix[colnames(logFCmatrixchosen)[i], colnames(logFCmatrixchosen)[j]]
  #  }
  #}
  cormatrix_2 <- reorder_cormat(cormatrix, method="average")
  cormatrix_2 = apply(cormatrix_2, 2, rev)
  upper_tri <- get_upper_tri(cormatrix_2)
  melted_cormat <- melt(cormatrix_2, na.rm = TRUE)
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  ggheatmap
  
  pdf(paste0("./plots/signatureplots/", name, "/verificationheatmap", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  significantgenematrix[is.na(significantgenematrix)] = 0
  significantgenematrix = significantgenematrix %>% rownames_to_column(var = "id")
  meltedshit = gather(significantgenematrix, dataset, logFC, -id)
  ggheatmap <- ggplot(meltedshit, aes(dataset, id, fill = logFC))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                         midpoint = 0, space = "Lab", 
                         name="logFC") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  ggheatmap
  pdf(paste0("./plots/signatureplots/", name, "/heatmapbygene", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  print(paste0("I'm done with the ", name, " signature."))
}

save(agingsignatures, file = "agingsignatures.RData")


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
