#master script for analyzing rnaseq data
library(tidyverse)
library(edgeR)
library(GEOquery)
library(limma)
library(ggrepel)

#>>>>> NEW TISSUE (for a new tissue, start here)
#$$$$$ get filtered phenodata and exprdata
# get phenodata
gse = getGEO("GSE#####") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))
# if you have the dataframe already, just name it filteredphenodata

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata = subset(filteredphenodata, subset = grepl('.*sedentary.*', filteredphenodata$title))
# takes only rows containing "sedentary" in the column "title"

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

# separate by tissue
# NONALPHANUMERIC
# (if everything is in one column separated with underscores or sth similar)
filteredphenodata = filteredphenodata %>% separate(V2, c("Tissue", "Age", NA))
filteredphenodata = subset(filteredphenodata, Tissue == "liver")
# SEPARATED ALREADY
# (if you have a separate column)
filteredphenodata = subset(filteredphenodata, Tissue == "liver")
# rename the column if necessary:
# filteredphenodata = filteredphenodata %>% dplyr::rename(Tissue=tissue.ch1)
# GREP
# (if it is bad and you need to grep)
filteredphenodata = subset(filteredphenodata, grepl("^Liver tissue.*", source_name_ch1))
# >>>>> OUTLIERS (remove them here)

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)


#$$$$$ filter genes with low expression
# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 10) # take only genes having more than
#                                                                10 reads...
filteredexprdata1 = filteredexprdata[filteredexprdata$rowsum > 3,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata1, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))


#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
source("FUN.Ensembl_mouse_dictionary_create.R")
dic = Ensembl_mouse_dictionary_create(filteredexprdata)
# CONVERT TRANSCRIPTS
source("FUN.Ensembl_mouse_dictionary_create_for_trans.R") #if rat or human, use the corresponding function
dic = Ensembl_mouse_dictionary_create_for_trans(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
source("FUN.RLE_normalization.R")
normdata = RLE_normalization(normdata)
normdata  = log2(normdata + 1)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# LOOKUP (entrez is already there, just take the sum)
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(sum, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))


#$$$$$ PCA
exprforpca = data.frame(t(scale(t(normdata))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = filteredphenodata[, "Age"]
cluster_plot + geom_point(aes(color = color_pca))
# OUTLIERS
# if boxes obscure something, use this:
# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))


#>>>>> NEW SEX (for a new sex, start here)
#$$$$$ Filter sex, scale before limma
# separate by sex
# SEPARATED ALREADY
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex == "Male")
# GREP
# (if it is bad and you need to grep)
sexyphenodata = subset(filteredphenodata, grepl("^Male.*", source_name_ch1))

# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

#!!!
# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata
#!!!

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))


#$$$$$ Limma
# CONTINUOUS (design matrix for continuous data)
# sexyphenodata$Age = sub("m$", "", sexyphenodata$Age) # delete nonnumeric characters
sexyphenodata$Age = as.numeric(as.character(sexyphenodata$Age))
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age")
print(design_matrix)

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "3m") # adjust for your age entries
design_matrix = model.matrix(~ currentfactor)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age") #the "Age" ones should be older
print(design_matrix)

# Launch limma
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
diffexprfit <- lmFit(sexyexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
target[["LogFC_table"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH", confint = TRUE)
# calculate standard error
target[["LogFC_table"]]$SE = (target[["LogFC_table"]]$CI.R - target[["LogFC_table"]]$CI.L) / 3.92


#$$$$$ Plot expression for top 1 diff expressed gene
# EXPR PLOT
topgenes = target[["LogFC_table"]]
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), as.numeric(sexyphenodata$Age)))
copps$V3 = as.numeric(as.character(copps$V3))
copps$V2 = as.double(as.character(copps$V2))
#copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE####"]][["Liver"]][["Male"]] <- target$LogFC_table
name = "Human_GSE5086_Muscle"
gender = "Male"
pdf(paste0("./plots/", name, "_", gender, "_processeddensity", ".pdf"))
print(target$Processed_density)
dev.off()
pdf(paste0("./plots/", name, "_rawdensity", ".pdf"))
print(target$Raw_density)
dev.off()
pdf(paste0("./plots/", name, "_pca", ".pdf"))
print(target$PCA)
dev.off()
pdf(paste0("./plots/", name, "_", gender, "_topgeneexpression", ".pdf"))
print(target$Topgene_expression)
dev.off()

# save progress
save(logFClist, file = "logFClist.RData")