##### collecting bad boys with NAs:

library(biomaRt)
library(tidyverse)
library(deming)
library(reshape2)
library(sgd)
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



#correlation matrix and heatmap
logFCmatrix = logFClist$Mouse$GSE6591$Lung$Male$DBA2J["logFC"]
logFCmatrix = logFCmatrix %>% rename(Mouse_GSE6581_Lung_Male_DBA2J = logFC)
logFCmatrix = merge(logFCmatrix, logFClist$Mouse$GSE6591$Lung$Male$C57BL6J["logFC"], by=0)
logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
logFCmatrix = logFCmatrix %>% rename(Mouse_GSE6581_Lung_Male_C57BL6J = logFC)
logFClist1 = logFClist
logFClist1$Mouse$GSE6591 = NULL
totalrownames = rownames(logFCmatrix)

#remove bad datasets:
logFClist1$Mouse$GSE53959 = NULL
logFClist1$Human$GSE40645 = NULL
logFClist1$Human$GSE5086 = NULL

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

badboys = subset(rownames(logFCmatrix), logFCmatrix$NACount >=35)

##### assembling the cortestsign matrix

cormethod = "spearman"
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
      if (cormatrixsign[rowname, colname] > 0.2){
        cortestsign[rowname, colname] = 1
      } else if (cormatrixsign[rowname, colname] < -0.2){
        cortestsign[rowname, colname] = -1
      } else {
        cortestsign[rowname, colname] = 0
      }
    } else{
      cortestsign[rowname, colname] = 0
    }
  }
}

# stacked barplot showing how many ones and minus ones are there for each data point (as proportions), and then names of bad data points: 
new.df<-melt(cortestsign)
new.df<-new.df[complete.cases(new.df),]
new.df$Score<-factor(new.df$value, levels = c("1", "-1", "0"))
new.df = new.df %>% filter(Score != 0)
ggplot(new.df, aes(x=variable, fill=Score)) + geom_bar(aes(y = (..count..)/58))

anothadf = cortestsign
anothadf$ones = rowSums(anothadf == 1)
anothadf$minusones = rowSums(anothadf == -1)
rownames(anothadf)[which(anothadf$ones <= anothadf$minusones)]
# density plot of number of ones across datasets:
ggplot(anothadf, aes(x = ones)) + geom_density()
# based on the visually established threshold, filter these bitches out:
rownames(anothadf)[which(anothadf$ones <= 10)]


##### heatmap, looking for bad datasets:
# Reorder the correlation matrix
cormatrix_2 <- reorder_cormat(cormatrixsign,method="average")
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

##### applying Deming regression:

for (el in logFCunlisted){ 
  el$SE = el$SE / sd(el$logFC)
  el$logFC = el$logFC / sd(el$logFC)
}

#totalrownamematrix = matrix(nrow = length(logFCunlisted), ncol = length(logFCunlisted))
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
    tempdata = matrix(nrow = length(totalrownames), ncol = 2)
    rownames(tempdata) = totalrownames
    tempdata[, 1] = logFCunlisted[[i]][totalrownames,]$logFC
    tempdata[, 2] = logFCunlisted[[j]][totalrownames,]$logFC
    tempdata = na.omit(tempdata)
    totalrownames = rownames(tempdata)
    totalrownamematrix[[names(logFCunlisted)[i]]][[names(logFCunlisted)[j]]] = totalrownames
  }
}

  fn = function(k_no_first){
    k = c()
    k[1] = 1
    k[2:length(logFCunlisted)] = k_no_first
  #  k[2:10] = k_no_first
    res = 0
    for (i in 1:(length(logFCunlisted)-1)){
  #  for (i in 1:9){
      namei = names(logFCunlisted)[i]
      for (j in (i + 1):length(logFCunlisted)){
  #    for (j in (i + 1):10){
        namej = names(logFCunlisted)[j]
        if (cortestsign[namei, namej] != 1){
          next
        }
        totalrownames = totalrownamematrix[[namei]][[namej]]
        ai = logFCunlisted[[i]][totalrownames,]$logFC
        aj = logFCunlisted[[j]][totalrownames,]$logFC
        res = res + sum(
          (((aj - (k[j]/k[i])*ai)^2)*((ai - (k[i]/k[j])*aj)^2))/
            (((aj - (k[j]/k[i])*ai)^2)+((ai - (k[i]/k[j])*aj)^2)))/length(totalrownames)
      }
    }
    return(res)
  }
  kvec = rnorm(length(logFCunlisted) - 1, 1, 1)
  #kvec = rnorm(9, 1, 1)
  ptm <- proc.time()
  optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B")
  proc.time() - ptm

kres = c(1, optimized$par)



# compare to deming:

dem_coefs <- c()
for (i in 1:9){
  for (j in (i + 1):10){
    if (cortestsign[names(logFCunlisted)[i], names(logFCunlisted)[j]] != 1){
      next
    }
    topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
    topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topA = topA %>% column_to_rownames(var = "row.names")
    topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
    topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
    topB = topB %>% column_to_rownames(var = "row.names")
    totalrownames = union(rownames(topA), rownames(topB))
    tempdata = matrix(nrow = length(totalrownames), ncol = 2)
    rownames(tempdata) = totalrownames
    tempdata[, 1] = logFCunlisted[[i]][totalrownames,]$logFC
    tempdata[, 2] = logFCunlisted[[j]][totalrownames,]$logFC
    tempdata = na.omit(tempdata)
    totalrownames = rownames(tempdata)
    plot1 = deming(logFCunlisted[[j]][totalrownames,]$logFC ~ logFCunlisted[[i]][totalrownames,]$logFC - 1)
    dem_coefs <- c(dem_coefs,plot1$coefficients[2])
    plot(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, main = paste0("Correlation: ", as.character(cormatrixsign[names(logFCunlisted)[i], names(logFCunlisted)[j]]), ", p-value: ", as.character(coradjpvalsign[names(logFCunlisted)[i], names(logFCunlisted)[j]])))
    abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
    abline(0, kres[j]/kres[i], col = "red", lwd = 2)
    abline(v = 0)
    abline(h = 0)
    plot3 = lm(logFCunlisted[[j]][totalrownames,]$logFC ~ logFCunlisted[[i]][totalrownames,]$logFC - 1)
    plot27 = lm(logFCunlisted[[i]][totalrownames,]$logFC ~ logFCunlisted[[j]][totalrownames,]$logFC - 1)
    abline(0, plot3$coef, col = "green", lwd = 2)
    abline(0, 1/plot27$coef, col = "darkgreen", lwd = 2)
  }
}

# normalize:

optimcoefs = optimized$par / optimized$par[1]

for (i in 1:length(logFCunlisted)){ 
  logFCunlisted[[i]]$SE = logFCunlisted[[i]]$SE / optimcoefs[i]
  logFCunlisted[[i]]$logFC = logFCunlisted[[i]]$logFC / optimcoefs[i]
}

# discard bad boys:
for (el in logFCunlisted){
  for (i in 1:length(rownames(el))){
    if (rownames(el)[i] %in% badboys){
      el = el[-i,]
    }
  }
}

# run mixed-effect model:

sourcedata = as.data.frame(names(logFCunlisted))
rownames(sourcedata) = names(logFCunlisted)
colnames(sourcedata) = "kekkekkek"
sourcedata = sourcedata %>% separate(kekkekkek, c(NA, "dataset", NA), sep = "_")

totalgenes = c()
for (el in logFCunlisted){
  totalgenes = union(totalgenes, rownames(el))
}
for (genename in totalgenes){
  logFC = c()
  SE = c()
  for (el in logFCunlisted){
    if (genename %in% rownames(el)){
      logFC = c(logFC, el$logFC[which(rownames(el) == genename)])
      SE = c(SE, el$SE[which(rownames(el) == genename)])
    }
  }
  mixedeffres = rma.mv(yi = logFC, V = SE, mods = ~ 1, method = "REML", random = list(~ 1 | as.factor(sourcedata$dataset)))
  
}


for(i in 1:5){
  for(j in (i+1):6){
    plot1 = deming(logFCunlisted[[j]][] ~ d[[i]] - 1)
    plot(d[[i]], d[[j]], xlim = c(-1, 1), ylim = c(-1, 1))
    plot(d[[i]], d[[j]])
    abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
    abline(0, plot2$par[j]/plot2$par[i], col = "red", lwd = 2)
    abline(v = 0)
    abline(h = 0)
    plot3 = lm(d[[j]] ~ d[[i]] - 1)
    plot27 = lm(d[[i]] ~ d[[j]] - 1)
    #abline(0, plot3$coef, col = "green", lwd = 2)
    #abline(0, 1/plot27$coef, col = "darkgreen", lwd = 2)
    dem_coefs <- c(dem_coefs,plot1$coefficients[2])
  }
}











d = list()
a = rnorm(1000, 0, 1)
b = 2*a+rnorm(1000, 0, 1)
b = b
b <- b/sd(b)
c = 3*a + rnorm(1000, 0, 1)
d = list(a, b, c)
d = list(a, b)
fn = function(k){
  res = 0
  for (i in 1:length(d)){
    for (j in i:length(d)){
      if (i != j){
        res = res + sum((((d[[j]] - (k[j]/k[i])*d[[i]])^2)*((d[[i]] - (k[i]/k[j])*d[[j]])^2))/(((d[[j]] - (k[j]/k[i])*d[[i]])^2)+((d[[i]] - (k[i]/k[j])*d[[j]])^2)))/1000
      }  
    }
  }
  return(res)
}
fn1 = function(k){
  sum((((a - k*b)^2)+((b-(1/k)*a)^2)))
}

plot2 = optim(rnorm(length(d),0,1), fn)
dem_coefs <- c()
for(i in 1:(length(d)-1)){
  for(j in (i+1):length(d)){
    plot1 = deming(d[[j]] ~ d[[i]] - 1)
    plot(d[[i]], d[[j]], xlim = c(-1, 1), ylim = c(-1, 1))
    plot(d[[i]], d[[j]])
    abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
    abline(0, plot2$par[j]/plot2$par[i], col = "red", lwd = 2)
    abline(v = 0)
    abline(h = 0)
    plot3 = lm(d[[j]] ~ d[[i]] - 1)
    plot27 = lm(d[[i]] ~ d[[j]] - 1)
    #abline(0, plot3$coef, col = "green", lwd = 2)
    #abline(0, 1/plot27$coef, col = "darkgreen", lwd = 2)
    dem_coefs <- c(dem_coefs,plot1$coefficients[2])
  }
}
plot2$par/plot2$par[1]
plot2$par[3]/plot2$par[2]
dem_coefs[2]/dem_coefs[1]

dem_coefs/plot2$par[2]*plot2$par[1]
plot1 = deming(d[[1]] ~ b - 1)
plot2 = optim(c(1, 1, 1), fn)
plot3 = lm(a ~ b - 1)
plot27 = lm(b ~ a - 1)
plot4 = optim(c(1, 1), fn1)
plot(d[[2]], d[[1]])
abline()
plot(b, a,xlim=c(-3,3),ylim=c(-3,3))
abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
abline(0, plot2$par, col = "red", lwd = 2)
#abline(0, plot3$coef, col = "green", lwd = 2)
abline(v = 0)
abline(h = 0)
abline(0, plot4$par, col = "darkred", lwd = 2)
abline(0, 1/plot27$coef, col = "darkgreen", lwd = 2)