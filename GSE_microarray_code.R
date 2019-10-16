#import libs
library(tidyverse)
library(dplyr)
library(ggplot2)
library(GEOquery)
#obtain phenodata and filter the liver tissue
mygse = getGEO("GSE11845")
mygsedf = data.frame(pData(mygse[[1]]))
filteredphenodata = subset(mygsedf, grepl("^Liver tissue.*", source_name_ch1))
#obtain expression data and filter the liver tissue
exprdata = data.frame(exprs(mygse[[1]]))
filteredexprdata = exprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)
#visualise the expression data
filteredexprstack = stack(filteredexprdata)
ggplot(filteredexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
#now scale the distribution of each gene expression in all the samples and visualize
scaledexprdata = data.frame(scale(filteredexprdata))
scaledexprstack = stack(scaledexprdata)
ggplot(scaledexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
#quantile normalize for a better dataset
library(preprocessCore)
scaledexprdatam = normalize.quantiles(as.matrix(scaledexprdata))
scaledexprdataqn = data.frame(scaledexprdatam)
colnames(scaledexprdataqn) = colnames(scaledexprdata)
rownames(scaledexprdataqn) = rownames(scaledexprdata)
scaledexprstackqn = stack(scaledexprdataqn)
ggplot(scaledexprstackqn, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
#the data has already been logarithmized, and inadequately expressed genes filtered out prior to that
#transform microarray IDs into Entrez IDs by taking mean values of log expression
exd <- scaledexprdataqn
lookup = fData(mygse[[1]])[, c('ID', 'Entrez_Gene_ID')]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=Entrez_Gene_ID)
new_exd <- exd %>% group_by(Entrez_Gene_ID) %>% summarize_each(mean, matches("^GSM.*"))
entrezexprdata <- new_exd %>% remove_rownames() %>% column_to_rownames(var = "Entrez_Gene_ID")
entrezexprstack = stack(entrezexprdata)
ggplot(entrezexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
#next I just checked that I did what I wanted by comparing the numbers:
#View(distinct(subset(filteredfdata, ID %in% rownames(filteredexprdata)), Entrez_Gene_ID, .keep_all = TRUE))
#View(exd[order(exd$Entrez_Gene_ID),])
#edit phenodata table for PCA
editedphenodata <- filteredphenodata
editedphenodata$source_name_ch1 = sub("^([^ ]+) ([^ ]+) ([^ ]+) (.*)", "\\4", editedphenodata$source_name_ch1)
library(xlsx)
library(readxl)
write.xlsx(editedphenodata, "editedphenodata.xlsx")
#...edit in Excel...
#***
#old version of the library (don't use it):
editedphenodata = read.xlsx("editedphenodata.xlsx", sheetName = "Sheet1", colnames = TRUE, rownames = TRUE)
editedphenodata <- editedphenodata %>% remove_rownames() %>% column_to_rownames(var = "NA.")
#on Windows and Linux:
editedphenodata = read_xlsx("editedphenodata.xlsx", sheet = "Sheet1")
editedphenodata <- editedphenodata %>% remove_rownames() %>% column_to_rownames(var = "...1")
#***
editedphenodataforpca = editedphenodata[editedphenodata$age != 18,]
entrezexprdataforpca = entrezexprdata[, rownames(editedphenodataforpca)]
#plot PCA
exprforpca = data.frame(t(scale(t(entrezexprdataforpca)))) # this way every gene will contribute evenly to pca
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = editedphenodataforpca[, "diet"]
shape_pca = editedphenodataforpca[, "resveratrol"]
cluster_plot + geom_point(aes(color = color_pca, shape = shape_pca))
#checking for other effective pcs
screeplot(pcamodel, npcs = 6, type = "barplot")
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
cluster_plot + geom_point(aes(color = color_pca, shape = shape_pca))
#analysing diff expression with limma
library(limma)
#before running limma, filter out genes with low standard dev across samples
library(genefilter)
sdevs = apply(entrezexprdata, 1, sd)
hist(sdevs) #check
filteredexprforlimma = entrezexprdata[sdevs > shorth(sdevs), ]
#limma analysis in a loop
top_genes1 = list()
design_matrices = list()
control <- editedphenodata$limma_indices[1]
for (limmaindex in subset(unique(editedphenodata$limma_indices), unique(editedphenodata$limma_indices) != control)){
  currentphenodata = subset(editedphenodata, editedphenodata$limma_indices == control | editedphenodata$limma_indices == limmaindex)
  currentexprdata = filteredexprforlimma[, rownames(currentphenodata)]
  currentfactor <- factor(currentphenodata$limma_indices)
  currentfactor = relevel(currentfactor, control)
  design_matrix = model.matrix(~ currentfactor)
  rownames(design_matrix) = rownames(currentphenodata)
  colnames(design_matrix) <- c("Intercept", limmaindex)
  design_matrices[[paste(limmaindex, "_vs_control", sep = "")]] = design_matrix
  contrast_matrix <- makeContrasts(paste(paste(limmaindex, "_vs_control", sep = ""), "=", limmaindex), levels=design_matrix)
  diffexprfit <- lmFit(currentexprdata, design_matrix)
  # diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
  diffexprfitcontrast <- eBayes(diffexprfit)
  top_genes1[[paste(limmaindex, "_vs_control", sep = "")]] <- topTable(diffexprfitcontrast, coef = limmaindex, number = Inf, adjust = "BH")
}
#temporary solution that worked, but had wrong total variance due to a large design matrix (hence the wrong number of small p-values)
age_diet_resveratrol_factor = factor(editedphenodata$Age_diet_resveratrol)
age_diet_resveratrol_factor = relevel(age_diet_resveratrol_factor, "27_standard_0")
design_matrix = model.matrix(~ age_diet_resveratrol_factor)
contrast_matrix <- makeContrasts(diet_EOD_vs_standard = age_diet_resveratrol_factor27_EOD_0, levels=design_matrix)
diffexprfit <- lmFit(qnentrezexprdata, design_matrix) #qnentrezexprdata is the same as entrezexprdata
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
top_genes <- topTable(diffexprfitcontrast, number = Inf, adjust = "BH")

age_diet_resveratrol_factor = factor(editedphenodata$limma_indices)
age_diet_resveratrol_factor = relevel(age_diet_resveratrol_factor, "a27_standard_0")
design_matrix = model.matrix(~ age_diet_resveratrol_factor)
rownames(design_matrix) = rownames(editedphenodata)
contrast_matrix <- makeContrasts(diet_EOD_vs_standard = age_diet_resveratrol_factora27_EOD_0, levels=design_matrix)
diffexprfit <- lmFit(entrezexprdata, design_matrix) #qnentrezexprdata is the same as entrezexprdata
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
topEOD1 <- topTable(diffexprfitcontrast, number = Inf, adjust = "BH")

#legacy code for limma


#first, create a design matrix
age_factor = factor(editedphenodata$age)
diet_factor = factor(editedphenodata$diet)
resveratrol_factor = factor(editedphenodata$resveratrol)
age_factor = relevel(age_factor, "27")
diet_factor = relevel(diet_factor, "standard diet")
resveratrol_factor = relevel(resveratrol_factor, "-")
design_matrix = model.matrix(~ age_factor + diet_factor + resveratrol_factor)
colnames(design_matrix) <- c("intercept", "age18", "dietEOD", "diethighcaloric", "resveratrol0.01", "resveratrol0.04")
#second, create a contrast matrix
contrast_matrix <- makeContrasts(age_18_vs_27 = age18, diet_EOD_vs_standard = dietEOD,
                                 diet_highcaloric_vs_standard = diethighcaloric, 
                                 diet_highcaloric_vs_EOD = diethighcaloric - dietEOD, 
                                 resveratrol_0.01_vs_0 = resveratrol0.01, 
                                 resveratrol_0.04_vs_0 = resveratrol0.04, 
                                 resveratrol_0.04_vs_0.01 = resveratrol0.04 - resveratrol0.01, 
                                 levels=design_matrix)
#finally, run limma's multiple component ANOVA
diffexprfit <- lmFit(newscaledexprdataqn, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
volcanoplotot(diffexprfitcontrast, coef = "age_18_vs_27") #you can substitute age18vs27 with any comparison group
top_genes <- topTable(diffexprfitcontrast, number = Inf, adjust = "BH") #you can add coef here too
diffexprresult <- decideTests(diffexprfitcontrast) #the threshold is 0.05 for adjusted p-value by default 
summary(diffexprresult)