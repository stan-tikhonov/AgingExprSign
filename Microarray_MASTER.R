#master script for analyzing microarray data
library(tidyverse)
library(preprocessCore)
library(GEOquery)
library(limma)


#>>>>> NEW TISSUE (for a new tissue, start here)
#$$$$$ get filtered phenodata and exprdata and fixed featuredata
# get phenodata
gse = getGEO("GSE11845") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))
# if you have the dataframe already, just name it filteredphenodata

# fix featuredata:
featuredata = fData(gse[[1]])
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata = subset(filteredphenodata, subset = grepl('.*sedentary.*', filteredphenodata$title))
# takes only rows containing "sedentary" in the column "title"

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
# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
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


#$$$$$ Scale
# SCALE (ordinary scaling)
scaledexprdata = data.frame(scale(logdata))
visualstack = stack(scaledexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# QUANTILE (quantile normalization, if needed)
normdata = normalize.quantiles(as.matrix(scaledexprdata))
normdata = data.frame(normdata)
colnames(normdata) = colnames(scaledexprdata)
rownames(normdata) = rownames(scaledexprdata)
rm(scaledexprdata)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# if it is not needed: normdata <- scaledexprdata


#$$$$$ Convert to Entrez (from Ensembl)(take mean of fluorescence)
# CONVERT (convert and take means)
source("FUN.Ensembl_mouse_dictionary_create.R") #if rat or human, use the corresponding function
dic = Ensembl_mouse_dictionary_create(normdata)
source("FUN.Ensembl_to_entrez_for_microarray.R")
normdata = Ensembl_to_entrez(normdata, dic)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata$ENTREZ_GENE_ID)) # get rid of multiple entrez per chip ID
exd <- normdata
lookup = featuredata[, c('ID', 'ENTREZ_GENE_ID')]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=ENTREZ_GENE_ID)
exd <- exd %>% group_by(ENTREZ_GENE_ID) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = "ENTREZ_GENE_ID")
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
colnames(design_matrix) <- c("Intercept", "Age")

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "3m") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age") #the "Age" ones should be older

# Launch limma
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
diffexprfit <- lmFit(sexyexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
logFClist[["Mouse"]][["GSE####"]][["Liver"]][["Male"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH", confint = TRUE)
# calculate standard error
logFClist$Mouse$GSE123981$Liver$Male$SE = (logFClist$Mouse$GSE123981$Liver$Male$CI.R - logFClist$Mouse$GSE123981$Liver$Male$CI.L) / 3.92





