#master script for analyzing rnaseq data
library(tidyverse)
library(edgeR)
library(GEOquery)
library(limma)


#>>>>> NEW TISSUE (for a new tissue, start here)
#$$$$$ get filtered phenodata and exprdata
# get phenodata
gse = getGEO("GSE11845") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))
# if you have the dataframe already, just name it filteredphenodata

# separate by tissue
# NONALPHANUMERIC
# (if everything is in one column separated with underscores or sth similar)
filteredphenodata = filteredphenodata %>% separate(V2, c("Tissue", "Age", NA))
filteredphenodata = subset(filteredphenodata, Tissue == "liver")
# SEPARATED ALREADY
# (if you have a separate column)
filteredphenodata = subset(filteredphenodata, Tissue == "liver")
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
# CHECK PEAK
# (first look if there is a peak at low values of expression)
logdata = log2(filteredexprdata1 + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
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
source("FUN.Ensembl_mouse_dictionary_create.R")
dic = Ensembl_mouse_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
source("FUN.RLE_normalization.R")
normdata = RLE_normalization(normdata)
normdata  = log2(normdata + 1)
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
cluster_plot + geom_point(aes(color = color_pca)) + geom_text(aes(label=rownames(cluster_values)),hjust=0, vjust=0)


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
scaledexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(scaledexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))


#$$$$$ Limma
# CONTINUOUS (design matrix for continuous data)
sexyphenodata$Age = sub("m$", "", sexyphenodata$Age) # delete nonnumeric characters
sexyphenodata$Age = as.integer(sexyphenodata$Age)
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("(Intercept)", "Age")

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "3m") # adjust for your age entries
design_matrix = model.matrix(~ currentfactor)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("(Intercept)", "Age") #the "Age" ones should be older

# Launch limma
contrast_matrix <- makeContrasts(paste("Age", "=", "Age"), levels=design_matrix)
diffexprfit <- lmFit(currentexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
logFClist[["GSE####"]][["Liver"]][["Male"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH")



