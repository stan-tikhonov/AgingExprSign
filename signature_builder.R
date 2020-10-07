# signature builder v. 2.0

library(biomaRt)
library(tidyverse)
library(deming)
library(reshape2)
library(metafor)

load("logFClist.RData")

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

pvalmatrix = matrix(nrow = length(totalrownames), ncol = length(logFCunlisted))
rownames(pvalmatrix) = totalrownames
colnames(pvalmatrix) = names(logFCunlisted)

for (name in names(logFCunlisted)){
  pvalmatrix[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$adj.P.Val
  logFCmatrix[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$logFC
  SEmatrixregr[rownames(logFCunlisted[[name]]), name] = logFCunlisted[[name]]$SE
}

pvalmatrix = as.data.frame(pvalmatrix)
logFCmatrix = as.data.frame(logFCmatrix)
SEmatrixregr = as.data.frame(SEmatrixregr)
logFCmatrixregr = logFCmatrix

# fix spaces in logFCmatrixregr colnames:
colnames(logFCmatrixregr) = sub(" ", "", colnames(logFCmatrixregr))
colnames(SEmatrixregr) = sub(" ", "", colnames(SEmatrixregr))

names(logFCunlisted) = sub(" ", "", names(logFCunlisted))

# determine the optimal threshold:
source("FUN_Correlation_loop_top_genes.r")
res = Correlation_loop_top_genes(logFCmatrix, pvalmatrix, seq(25,1000,by=25))

# make cortestsign matrix:

source("FUN.Cormatricesmaker.R")
cormethod = "pearson"
thres = "250"
res = cormatricesmaker(logFCunlisted, cormethod, signifgenesthres = thres)
cortestsign = res$cortestsign
cormatrixsign = res$cormatrixsign
coradjpvalsign = res$coradjpvalsign
corpvalsign = res$corpvalsign

# normalize by sd:

for (i in 1:length(colnames(logFCmatrixregr))){
  SEmatrixregr[,i] = SEmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
  logFCmatrixregr[,i] = logFCmatrixregr[,i] / sd(na.omit(logFCmatrixregr[,i]))
}

# make totalrownamematrix:
source("FUN.Totalrownamemaker.R")
totalrownamematrix = totalrownamemaker(logFCunlisted, thres)

# get source table:
sourcedata = as.data.frame(colnames(logFCmatrixregr))
rownames(sourcedata) = colnames(logFCmatrixregr)
colnames(sourcedata) = "kekkekkek"
sourcedata = sourcedata %>% separate(kekkekkek, c(NA, "dataset", NA), sep = "_")

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

# assessing how many genes are differentially expressed per signature:
corgenes = list()
for (name in names(chosencols)){
  lol = as.data.frame(rep(0, length(rownames(logFCmatrixregr))))
  lol = cbind(lol, rep(0, length(rownames(logFCmatrixregr))))
  rownames(lol) = rownames(logFCmatrixregr)
  colnames(lol) = c("countspositive", "countsnegative")
  for (i in 1:length(chosencols[[name]])){
    positiverownames = rownames(subset(logFCunlisted[[chosencols[[name]][i]]], logFCunlisted[[chosencols[[name]][i]]]$adj.P.Val < 0.05 & logFCunlisted[[chosencols[[name]][i]]]$logFC > 0))
    negativerownames = rownames(subset(logFCunlisted[[chosencols[[name]][i]]], logFCunlisted[[chosencols[[name]][i]]]$adj.P.Val < 0.05 & logFCunlisted[[chosencols[[name]][i]]]$logFC < 0))
    lol[positiverownames,]$countspositive = lol[positiverownames,]$countspositive + 1
    lol[negativerownames,]$countsnegative = lol[negativerownames,]$countsnegative + 1
  }
  corgenes[[name]] = rownames(subset(lol, lol$countspositive > length(chosencols[[name]]) * 0.2 | lol$countsnegative > length(chosencols[[name]]) * 0.2))
}

##### THIS IS FOR MANY SIGNATURES BUT ONE MINIMIZATION RUN FOR EACH SIGNATURE

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
source("FUN.Corheatmapper.R")

agingsignatures = list()
deminglist = list()
deminglist[["Human"]] = list()
deminglist[["Mouse"]] = list()
deminglist[["Rat"]] = list()
deminglist[["Muscle"]] = list()
deminglist[["Brain"]] = list()
deminglist[["Liver"]] = list()
deminglist[["All"]] = list()

# minimization loop:
for (name in names(chosencols)){
  # filter datasets for the individual signature:
  logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
  # minimize:
  minimums = c()
  for (i in 1:10){
    deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen, totalrownamematrix)
    minimums = c(minimums, deminglist[[name]][[i]]$minimum)
  }
  print(paste0("I'm done with ", name))
} 

# or load its results:
load("deminglist.RData")

# visualisation of coef distributions:
for (name in names(deminglist)){
  # here a is coefs and b is minimums
  a = as.data.frame(deminglist[[name]][[1]]$coefs)
  for (i in 2:length(deminglist[[name]])){
    a = cbind(a, deminglist[[name]][[i]]$coefs)
  }
  colnames(a) = 1:10
  a$group <- row.names(a)
  a.m <- melt(a, id.vars = "group")
  b = data.frame()
  for (i in 1:10){
    b = rbind(b, deminglist[[name]][[i]]$minimum)
  }
  b = cbind(b, 1:10)
  colnames(b) = c("minimum", "variable")
  temp = cortestsign[chosencols[[name]], chosencols[[name]]]
  normcoef = sum(temp[upper.tri(temp, diag = F)] == 1)
  b$minimum = b$minimum / normcoef
  b$variable = as.factor(b$variable)
  a.m$minimum = left_join(a.m, b, by = "variable")
  print(ggplot(a.m$minimum, aes(group, value)) + geom_boxplot() + geom_jitter(aes(color = as.factor(round(minimum, 4)))))
}

##### MAIN LOOP #####

# main loop:
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
  ggheatmap = corheatmapper(cormatrix, cormethod = "spearman")
  print(ggheatmap)
  pdf(paste0("./newplots/signatureplots1/", name, "/initialheatmap", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # get the deming coefs:
  minimums = c()
  for (i in 1:10){
    #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
    minimums = c(minimums, deminglist[[name]][[i]]$minimum)
  }
  kres = deminglist[[name]][[which.min(minimums)]]$coefs
  
  # plot exapmles:
  kekmatrix = cormatrixsign[chosencols[[name]], chosencols[[name]]] %>% rownames_to_column("datasetid1")
  kekmatrix = gather(kekmatrix, datasetid2, corvalue, -datasetid1)
  kekmatrix$corvalue = as.numeric(as.character(kekmatrix$corvalue))
  kekmatrix = kekmatrix %>% filter(corvalue != 1) %>% top_n(10, corvalue) %>% distinct(corvalue, .keep_all = T)
  
  for (i in 1:length(rownames(kekmatrix))){
    plot1 = deming(logFCmatrixchosen[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid2"]] ~ logFCmatrixchosen[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid1"]] - 1)
    #plot(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],1], logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],2], main = paste0("Correlation: ", as.character(cormatrixsign[colnames(logFCmatrixchosen)[1], colnames(logFCmatrixchosen)[2]]), ", p-value: ", as.character(coradjpvalsign[colnames(logFCmatrixchosen)[1], colnames(logFCmatrixchosen)[2]])))
    #abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
    #abline(0, kres[2]/kres[1], col = "red", lwd = 2)
    
    ggheatmap = ggplot(logFCmatrixchosen[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],], aes_string(x = kekmatrix[i,"datasetid1"], y = kekmatrix[i,"datasetid2"])) + geom_point() +
      geom_abline(slope = plot1$coefficients[2], intercept = 0, colour = "blue", size = 1) + 
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      geom_abline(slope = kres[which(colnames(logFCmatrixchosen) == kekmatrix[i, "datasetid2"])]/kres[which(colnames(logFCmatrixchosen) == kekmatrix[i, "datasetid1"])], intercept = 0, colour = "red", size = 1)
    print(ggheatmap)
    pdf(paste0("./newplots/signatureplots1/", name, "/demingexample", i, ".pdf"))
    print(ggheatmap)
    dev.off()
  }
  
  # normalize by deming coefficients:
  
  for (i in 1:length(colnames(logFCmatrixchosen))){
    SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
    logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
  }
  
  # discard bad boys:
  logFCmatrixchosen$NACount = rowSums(is.na(logFCmatrixchosen))
  if (name != "Liver"){
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = floor(length(colnames(logFCmatrixchosen))/2))
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < floor(length(colnames(logFCmatrixchosen))/2))
  } else {
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = 4)
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < 4)
  }
  logFCmatrixchosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %in% goodboys)
  SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %in% goodboys)
  logFCmatrixchosen$NACount = NULL
  
  # run mixed-effect model:
  agingsignatures[[name]] = signature_builder(logFCmatrixchosen, SEmatrixchosen)
  # plot examples:
  geneids = agingsignatures[[name]] %>% rownames_to_column("Row.names") %>% top_n(-5, adj_pval) %>% column_to_rownames("Row.names")
  geneids = rownames(geneids)
  for (i in 1:length(geneids)){
    helpertable = as.data.frame(t(logFCmatrixchosen[geneids[i],]))
    rownames(helpertable) = colnames(logFCmatrixchosen)
    colnames(helpertable) = c("logFC")
    helpertable$SE = t(SEmatrixchosen[geneids[i],])
    helpertable$source = as.factor(sourcedata[rownames(helpertable),"dataset"])
    helpertable$dataset = rownames(helpertable)
    helpertable = na.omit(helpertable)
    border = max(abs(helpertable$logFC)) + max(helpertable$SE)
    ggheatmap = ggplot(helpertable, aes(x = dataset, y = logFC, color = source)) + geom_pointrange(aes(ymin = logFC - SE, ymax = logFC + SE)) + geom_hline(yintercept = agingsignatures[[name]][geneids[i], "logFC"], colour = "red") +
      geom_hline(yintercept = 0) + ylim(-border, border)
    ggheatmap
    pdf(paste0("./newplots/signatureplots1/", name, "/mixedmodelexample", i, ".pdf"))
    print(ggheatmap)
    dev.off()
  }
  
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
  ggheatmap = corheatmapper(cormatrix, cormethod = "spearman")
  print(ggheatmap)
  pdf(paste0("./newplots/signatureplots1/", name, "/verificationheatmap", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # plot a verification heatmap by gene
  significantgenematrix[is.na(significantgenematrix)] = 0
  #clusteredshit = reorder_cormat(significantgenematrix)
  #dd <- as.dist((1-cor(t(significantgenematrix), method = "spearman"))/2)
  dd = dist(significantgenematrix, method = "manhattan")
  hc <- hclust(dd,method = "average")
  clusteredshit <-significantgenematrix[hc$order,]
  clusteredshit = clusteredshit %>% rownames_to_column(var = "id")
  meltedshit = gather(clusteredshit, dataset, logFC, -id, factor_key = T)
  meltedshit$id = factor(meltedshit$id, levels = as.character(meltedshit$id[1:length(rownames(significantgenematrix))]))
  ggheatmap <- ggplot(meltedshit, aes(dataset, id, fill = logFC))+
    geom_tile()+
    scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                         midpoint = 0, space = "Lab", 
                         name="LogFC")
  ggheatmap
  pdf(paste0("./newplots/signatureplots1/", name, "/heatmapbygene", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  print(paste0("I'm done with the ", name, " signature."))
}
save(agingsignatures, file = "agingsignatures.RData")

# this is for LOO and robust signature (no plots):

source("FUN.Signature_builder_4.R")
agingsignatures_v3 = list()
# main loop:
for (name in names(chosencols)){
  # filter datasets for the individual signature:
  logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
  SEmatrixchosen = SEmatrixregr[, chosencols[[name]]]
  
  # get the deming coefs:
  minimums = c()
  for (i in 1:10){
    #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
    minimums = c(minimums, deminglist[[name]][[i]]$minimum)
  }
  kres = deminglist[[name]][[which.min(minimums)]]$coefs
  
  # normalize by deming coefficients:
  
  for (i in 1:length(colnames(logFCmatrixchosen))){
    SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
    logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
  }
  
  # discard bad boys:
  logFCmatrixchosen$NACount = rowSums(is.na(logFCmatrixchosen))
  if (name != "Liver"){
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = floor(length(colnames(logFCmatrixchosen))/2))
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < floor(length(colnames(logFCmatrixchosen))/2))
  } else {
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = 4)
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < 4)
  }
  logFCmatrixchosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %in% goodboys)
  SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %in% goodboys)
  logFCmatrixchosen$NACount = NULL
  
  # run mixed-effect model:
  agingsignatures_v3[[name]] = signature_builder(logFCmatrixchosen, SEmatrixchosen, name)
  print(paste0("I'm done with ", name))
}
save(agingsignatures_v3, file = "agingsignatures_v4.RData")




# this is only for plots:
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
  ggheatmap = corheatmapper(cormatrix, cormethod = "spearman")
  print(ggheatmap)
  pdf(paste0("./newplots/signatureplots1/", name, "/initialheatmap", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # get the deming coefs:
  minimums = c()
  for (i in 1:10){
    #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
    minimums = c(minimums, deminglist[[name]][[i]]$minimum)
  }
  kres = deminglist[[name]][[which.min(minimums)]]$coefs
  
  # plot exapmles:
  kekmatrix = cormatrixsign[chosencols[[name]], chosencols[[name]]] %>% rownames_to_column("datasetid1")
  kekmatrix = gather(kekmatrix, datasetid2, corvalue, -datasetid1)
  kekmatrix$corvalue = as.numeric(as.character(kekmatrix$corvalue))
  kekmatrix = kekmatrix %>% filter(corvalue != 1) %>% top_n(10, corvalue) %>% distinct(corvalue, .keep_all = T)
  
  for (i in 1:length(rownames(kekmatrix))){
    plot1 = deming(logFCmatrixchosen[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid2"]] ~ logFCmatrixchosen[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],kekmatrix[i,"datasetid1"]] - 1)
    #plot(logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],1], logFCmatrixchosen[totalrownamematrix[[colnames(logFCmatrixchosen)[1]]][[colnames(logFCmatrixchosen)[2]]],2], main = paste0("Correlation: ", as.character(cormatrixsign[colnames(logFCmatrixchosen)[1], colnames(logFCmatrixchosen)[2]]), ", p-value: ", as.character(coradjpvalsign[colnames(logFCmatrixchosen)[1], colnames(logFCmatrixchosen)[2]])))
    #abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
    #abline(0, kres[2]/kres[1], col = "red", lwd = 2)
    
    ggheatmap = ggplot(logFCmatrixchosen[totalrownamematrix[[kekmatrix[i,"datasetid1"]]][[kekmatrix[i,"datasetid2"]]],], aes_string(x = kekmatrix[i,"datasetid1"], y = kekmatrix[i,"datasetid2"])) + geom_point() +
      geom_abline(slope = plot1$coefficients[2], intercept = 0, colour = "blue", size = 1) + 
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
      geom_abline(slope = kres[which(colnames(logFCmatrixchosen) == kekmatrix[i, "datasetid2"])]/kres[which(colnames(logFCmatrixchosen) == kekmatrix[i, "datasetid1"])], intercept = 0, colour = "red", size = 1)
    print(ggheatmap)
    pdf(paste0("./newplots/signatureplots1/", name, "/demingexample", i, ".pdf"))
    print(ggheatmap)
    dev.off()
  }
  
  # normalize by deming coefficients:
  
  for (i in 1:length(colnames(logFCmatrixchosen))){
    SEmatrixchosen[,i] = SEmatrixchosen[,i] / kres[i]
    logFCmatrixchosen[,i] = logFCmatrixchosen[,i] / kres[i]
  }
  
  # discard bad boys:
  logFCmatrixchosen$NACount = rowSums(is.na(logFCmatrixchosen))
  if (name != "Liver"){
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = floor(length(colnames(logFCmatrixchosen))/2))
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < floor(length(colnames(logFCmatrixchosen))/2))
  } else {
    ggplot(logFCmatrixchosen, aes(x = NACount)) + geom_density() + geom_vline(xintercept = 4)
    goodboys = subset(rownames(logFCmatrixchosen), logFCmatrixchosen$NACount < 4)
  }
  logFCmatrixchosen = subset(logFCmatrixchosen, rownames(logFCmatrixchosen) %in% goodboys)
  SEmatrixchosen = subset(SEmatrixchosen, rownames(SEmatrixchosen) %in% goodboys)
  logFCmatrixchosen$NACount = NULL
  
  helpertablelist = list()
  # plot examples:
  geneids = agingsignatures_v3[[name]] %>% rownames_to_column("Row.names") %>% top_n(-5, adj_pval) %>% column_to_rownames("Row.names")
  geneids = rownames(geneids)
  for (i in 1:length(geneids)){
    helpertable = as.data.frame(t(logFCmatrixchosen[geneids[i],]))
    rownames(helpertable) = colnames(logFCmatrixchosen)
    colnames(helpertable) = c("logFC")
    helpertable$SE = t(SEmatrixchosen[geneids[i],])
    helpertable$source = as.factor(sourcedata[rownames(helpertable),"dataset"])
    helpertable$dataset = rownames(helpertable)
    helpertable$species = sub("^([^_]+)_(.*)", "\\1", rownames(helpertable))
    helpertable$tissue = sub("^([^_]+)_([^_]+)_([^_]+)_(.*)", "\\3", rownames(helpertable))
    helpertable$tissue = sub("Frontalcortex", "Brain", helpertable$tissue)
    helpertable$tissue = sub("Cerebellum", "Brain", helpertable$tissue)
    helpertable = na.omit(helpertable)
    helpertable = helpertable %>% arrange(dataset)
    helpertablelist[[geneids[i]]] = helpertable
    border = max(abs(helpertable$logFC)) + max(helpertable$SE)
    ggheatmap = ggplot(helpertable, aes(x = dataset, y = logFC, color = tissue, shape = species, label = source)) + geom_pointrange(aes(ymin = logFC - SE, ymax = logFC + SE)) + geom_hline(yintercept = agingsignatures_v3[[name]][geneids[i], "logFC"], colour = "red") +
      geom_hline(yintercept = 0) + ylim(-border, border) + theme_minimal()# + geom_text_repel(aes(label=source))
    print(ggheatmap)
    pdf(paste0("./newplots/signatureplots1/", name, "/mixedmodelexample_shapes", geneids[i], ".pdf"))
    print(ggheatmap)
    dev.off()
  }
  
  # plot verification cor heatmap:
  logFCmatrixchosen = merge(logFCmatrixchosen, agingsignatures_v3[[name]]["logFC"], by = "row.names", all = TRUE)
  colnames(logFCmatrixchosen)[which(colnames(logFCmatrixchosen) == "logFC")] = paste0(name, "_signature")
  logFCmatrixchosen = logFCmatrixchosen %>% column_to_rownames(var = "Row.names")
  significantgenematrix = logFCmatrixchosen[rownames(subset(agingsignatures_v3[[name]], adj_pval < 0.05)),]
  #significantgenematrix = logFCmatrixchosen[rownames(subset(agingsignatures_v3[[name]], adj_pval_LOO < 0.05)),]
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
  ggheatmap = corheatmapper(cormatrix, cormethod = "spearman")
  print(ggheatmap)
  pdf(paste0("./newplots/signatureplots1/", name, "/verificationheatmap", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  # plot a verification heatmap by gene
  significantgenematrix[is.na(significantgenematrix)] = 0
  #clusteredshit = reorder_cormat(significantgenematrix)
  #dd <- as.dist((1-cor(t(significantgenematrix), method = "spearman"))/2)
  dd = dist(significantgenematrix, method = "manhattan")
  hc <- hclust(dd,method = "average")
  clusteredshit <-significantgenematrix[hc$order,]
  clusteredshit = clusteredshit %>% rownames_to_column(var = "id")
  meltedshit = gather(clusteredshit, dataset, logFC, -id, factor_key = T)
  meltedshit$id = factor(meltedshit$id, levels = as.character(meltedshit$id[1:length(rownames(significantgenematrix))]))
  ggheatmap <- ggplot(meltedshit, aes(dataset, id, fill = logFC))+
    geom_tile()+
    scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                         midpoint = 0, space = "Lab", 
                         name="LogFC")
  ggheatmap
  pdf(paste0("./newplots/signatureplots1/", name, "/heatmapbygene", ".pdf"))
  print(ggheatmap)
  dev.off()
  
  print(paste0("I'm done with the ", name, " signature."))
}




################################################ L E G A C Y   C O D E ##########################################

# metafor no nothing:
agingmetafornonothing = list()

signature_metafornonothing = function(logFCmatrixregr){
  goodgenes = c()
  signature = data.frame()
  genenumber = 0
  for (genename in rownames(logFCmatrixregr)){
    genenumber = genenumber + 1
    percentready = (genenumber/length(rownames(logFCmatrixregr))) * 100
    if (genenumber %% 1000 == 0){
      print(paste0("I'm on gene No. ", genenumber, " (", round(percentready, 2), "% done)"))
    }
    logFC = logFCmatrixregr[genename,]
    logFC = logFC[!is.na(logFC)]
    SE = SEmatrixregr[genename,]
    SE = SE[!is.na(SE)]
    sourcevec = as.factor(sourcedata[colnames(logFCmatrixregr)[!is.na(logFCmatrixregr[genename,])],])
    #SE = rep(1, length(logFC))
    
    tryCatch(
      {
        mixedeffres = rma.mv(yi = logFC, V = SE ^ 2, method = "REML", random = list(~ 1 | sourcevec))
        signature = rbind(signature, c(mixedeffres$b[1], mixedeffres$pval))
        goodgenes = c(goodgenes, genename)
      },
      error=function(cond) {
        message("Fucked up")
        message("Here's the original error message:")
        message(cond)
      }
    )
  }
  #rownames(signature) = totalgenes[-which(totalgenes %in% badgenes)]
  rownames(signature) = goodgenes
  colnames(signature) = c("logFC", "pval")
  
  signature$adj_pval = p.adjust(signature$pval, method = "BH")
  return(signature)
}

for (name in names(chosencols)){
  # filter datasets for the individual signature:
  logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
  SEmatrixchosen = SEmatrixregr[, chosencols[[name]]]
  
  # get the deming coefs:
  minimums = c()
  for (i in 1:10){
    #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
    minimums = c(minimums, deminglist[[name]][[i]]$minimum)
  }
  kres = deminglist[[name]][[which.min(minimums)]]$coefs
  
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
  
  agingmetafornonothing[[name]] = signature_metafornonothing(logFCmatrixchosen)
}

# z-test for each gene:
agingztests = list()

signature_ttester = function(logFCmatrixregr){
  goodgenes = c()
  signature = data.frame()
  genenumber = 0
  for (genename in rownames(logFCmatrixregr)){
    genenumber = genenumber + 1
    percentready = (genenumber/length(rownames(logFCmatrixregr))) * 100
    if (genenumber %% 1000 == 0){
      print(paste0("I'm on gene No. ", genenumber, " (", round(percentready, 2), "% done)"))
    }
    logFC = logFCmatrixregr[genename,]
    logFC = logFC[!is.na(logFC)]
    ttestres = t.test(logFC)
    signature = rbind(signature, c(ttestres$estimate, ttestres$p.value))
    goodgenes = c(goodgenes, genename)
  }
  #rownames(signature) = totalgenes[-which(totalgenes %in% badgenes)]
  rownames(signature) = goodgenes
  colnames(signature) = c("logFC", "pval")
  
  signature$adj_pval = p.adjust(signature$pval, method = "BH")
  return(signature)
}

for (name in names(chosencols)){
  # filter datasets for the individual signature:
  logFCmatrixchosen = logFCmatrixregr[, chosencols[[name]]]
  SEmatrixchosen = SEmatrixregr[, chosencols[[name]]]
  
  # get the deming coefs:
  minimums = c()
  for (i in 1:10){
    #deminglist[[name]][[i]] = deming_minimizer(logFCmatrixchosen)
    minimums = c(minimums, deminglist[[name]][[i]]$minimum)
  }
  kres = deminglist[[name]][[which.min(minimums)]]$coefs
  
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
  
  agingztests[[name]] = signature_ttester(logFCmatrixchosen)
}

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
