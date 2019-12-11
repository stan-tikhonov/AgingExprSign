##### GSE5086 (BEGIN)

# get phenodata
gse = getGEO("GSE5086") #here, type your gse id
# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

load("filteredphenodataforGSE5086.RData") #Here I filtered out those who took medications or had any disease

filteredphenodata = filteredphenodata1

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#PCA
exprforpca = data.frame(t(scale(t(normdata))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = filteredphenodata[, "Age"]
cluster_plot + geom_point(aes(color = color_pca))
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

#analyze males:
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "Male")
# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Limma
# CONTINUOUS (design matrix for continuous data)
sexyphenodata$Age = as.numeric(as.character(sexyphenodata$Age))
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age")
print(design_matrix)

# Launch limma
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
diffexprfit <- lmFit(sexyexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
target[["LogFC_table"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH", confint = TRUE)
# calculate standard error
target[["LogFC_table"]]$SE = (target[["LogFC_table"]]$CI.R - target[["LogFC_table"]]$CI.L) / 3.92

# EXPR PLOT
topgenes = target[["LogFC_table"]]
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), as.integer(sexyphenodata$Age)))
copps$V3 = as.integer(as.character(copps$V3))
copps$V2 = as.double(as.character(copps$V2))
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in years", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

logFClist[["Human"]][["GSE5086"]][["Muscle"]][["Male"]] <- target$LogFC_table
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

#analyze females:
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "Female")
# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Limma
# CONTINUOUS (design matrix for continuous data)
sexyphenodata$Age = as.numeric(as.character(sexyphenodata$Age))
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age")
print(design_matrix)

# Launch limma
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
diffexprfit <- lmFit(sexyexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
target[["LogFC_table"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH", confint = TRUE)
# calculate standard error
target[["LogFC_table"]]$SE = (target[["LogFC_table"]]$CI.R - target[["LogFC_table"]]$CI.L) / 3.92

# EXPR PLOT
topgenes = target[["LogFC_table"]]
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), as.integer(sexyphenodata$Age)))
copps$V3 = as.integer(as.character(copps$V3))
copps$V2 = as.double(as.character(copps$V2))
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in years", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

logFClist[["Human"]][["GSE5086"]][["Muscle"]][["Female"]] <- target$LogFC_table
name = "Human_GSE5086_Muscle"
gender = "Female"
pdf(paste0("./plots/", name, "_", gender, "_processeddensity", ".pdf"))
print(target$Processed_density)
dev.off()
pdf(paste0("./plots/", name, "_", gender, "_topgeneexpression", ".pdf"))
print(target$Topgene_expression)
dev.off()
# save progress
save(logFClist, file = "logFClist.RData")

##### GSE5086 (END)

##### GSE9103 (BEGIN)

gse = getGEO("GSE9103") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

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
filteredphenodata$Age = sub("", "", filteredphenodata$title)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

filteredphenodata$Age = c("old", "old", "old", "old", "old", "old", "old", "old", "old", "old", "young", "young", "young", "young", "young", "young", "young", "young", "young", "young")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Limma
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "young") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), sexyphenodata$Age))
copps$V2 = as.double(as.character(copps$V2))
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in years (18-30 and 59-76)", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Human"]][["GSE9103"]][["Muscle"]][["NoSex"]] <- target$LogFC_table
name = "Human_GSE9103_Muscle"
gender = "NoSex"
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

##### GSE9103 (END)

##### GSE3150 (BEGIN)

# get phenodata
gse = getGEO("GSE3150") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata = subset(filteredphenodata, subset = grepl('.*Wild Type.*', filteredphenodata$title))

filteredphenodata$Age = c(22, 22, 22, 22, 22, 22, 10, 10, 10, 10, 10, 4, 4, 4, 4, 4)

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ PCA
exprforpca = data.frame(t(scale(t(normdata))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = as.factor(filteredphenodata[, "Age"])
cluster_plot + geom_point(aes(color = color_pca))
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Limma
# CONTINUOUS (design matrix for continuous data)
sexyphenodata$Age = as.numeric(as.character(sexyphenodata$Age))
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age")
print(design_matrix)

# Launch limma
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
diffexprfit <- lmFit(sexyexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
target[["LogFC_table"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH", confint = TRUE)
# calculate standard error
target[["LogFC_table"]]$SE = (target[["LogFC_table"]]$CI.R - target[["LogFC_table"]]$CI.L) / 3.92

# EXPR PLOT
topgenes = target[["LogFC_table"]]
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), as.integer(sexyphenodata$Age)))
copps$V3 = as.integer(as.character(copps$V3))
copps$V2 = as.double(as.character(copps$V2))
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE3150"]][["Liver"]][["Male"]] <- target$LogFC_table
name = "Mouse_GSE3150_Liver"
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

##### GSE3150 (END)

##### GSE53959 (BEGIN)

# get phenodata
gse = getGEO("GSE53959") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENSEMBL_ID"
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CONVERT TRANSCRIPTS
source("FUN.Ensembl_mouse_dictionary_create_for_trans.R") #if rat or human, use the corresponding function
dic = Ensembl_mouse_dictionary_create_for_trans(normdata)
# dic done here
source("FUN.Ensembl_to_entrez_for_microarray.R")
normdata = Ensembl_to_entrez_for_microarray(normdata, dic)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# One gene has sd == 0 across samples
# Its IDs: A_51_P205779, ENSMUST00000015998, and entrez 11801
# how do you know that? code:
which(apply(normdata, 1, sd) == 0)
normdata[584,]
xx <- as.list(org.Mm.egENSEMBLTRANS)
xx[["11801"]]
which(featuredata$ENSEMBL_ID == "ENSMUST00000015998")
featuredata[18009,]
filteredexprdata["A_51_P205779",]

# take care of this motherfucker:
filteredexprdata = filteredexprdata[-18009,]

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENSEMBL_ID"
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CONVERT TRANSCRIPTS
source("FUN.Ensembl_mouse_dictionary_create_for_trans.R") #if rat or human, use the corresponding function
dic = Ensembl_mouse_dictionary_create_for_trans(normdata)
# dic done here
source("FUN.Ensembl_to_entrez_for_microarray.R")
normdata = Ensembl_to_entrez_for_microarray(normdata, dic)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

exprforpca = data.frame(t(scale(t(normdata))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = filteredphenodata[, "Age"]
cluster_plot + geom_point(aes(color = color_pca))
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Limma
# CONTINUOUS (design matrix for continuous data)
sexyphenodata$Age = as.numeric(as.character(sexyphenodata$Age))
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age")
print(design_matrix)

# Launch limma
contrast_matrix <- makeContrasts("Age", levels=design_matrix)
diffexprfit <- lmFit(sexyexprdata, design_matrix)
diffexprfitcontrast <- contrasts.fit(diffexprfit, contrast_matrix)
diffexprfitcontrast <- eBayes(diffexprfitcontrast)
target[["LogFC_table"]] <- topTable(diffexprfitcontrast, coef = "Age", number = Inf, adjust = "BH", confint = TRUE)
# calculate standard error
target[["LogFC_table"]]$SE = (target[["LogFC_table"]]$CI.R - target[["LogFC_table"]]$CI.L) / 3.92

# EXPR PLOT
topgenes = target[["LogFC_table"]]
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), as.numeric(sexyphenodata$Age)))
copps$V3 = as.numeric(as.character(copps$V3))
copps$V2 = as.double(as.character(copps$V2))
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE53959"]][["Muscle"]][["Female"]] <- target$LogFC_table
gender = "Female"
name = "Mouse_GSE53959_Muscle"

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

##### GSE53959 (END)

##### E-MTAB-3374 (BEGIN)

# obtain the dataframes:
filteredexprdata = read.table("E-MTAB-3374-expr.txt", header = T)
filteredexprdata = filteredexprdata %>% remove_rownames() %>% column_to_rownames(var = "PROBE_ID")
for (i in colnames(filteredexprdata)){
  filteredexprdata[,i] = as.vector(filteredexprdata[,i])
}
filteredphenodata = read.table("E-MTAB-3374.sdrf.txt", sep = "\t", header = T)
filteredphenodata = filteredphenodata %>% remove_rownames() %>% column_to_rownames(var = "Source.Name")
for (i in colnames(filteredphenodata)){
  filteredphenodata[,i] = as.vector(filteredphenodata[,i])
}
featuredata = read.csv("E-MTAB-3374-featuredata.csv")
for (i in colnames(featuredata)){
  featuredata[,i] = as.vector(featuredata[,i])
}

# set up target list element:
target = list()

# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata = subset(filteredphenodata, subset = grepl('.*no.*', filteredphenodata$Characteristics.treatment.))

# get age:
filteredphenodata$Age = as.character(filteredphenodata$Factor.Value.age.)

# filter exprdata:
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]

#$$$$$ filter genes with low expression
# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

#entrez
entrezcolname = "Reporter.Database.Entry.entrez."
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
featuredata = featuredata %>% dplyr::rename(ID = Comment.AEReporterName.)
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_at(1:6, mean)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#limma
# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "2.5") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["E-MTAB-3374"]][["Liver"]][["Female"]] <- target$LogFC_table
name = "Human_E-MTAB-3374_Liver"
gender = "Female"
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

##### E-MTAB-3374 (END)

##### GSE6591 (BEGIN)

# get phenodata
gse = getGEO("GSE6591") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

filteredphenodata$Age = as.character(c(2, 2, 2, 18, 18, 18, 26, 26, 26, 2, 2, 2, 18, 18, 18))

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

sexyphenodata = subset(filteredphenodata, subset = grepl('.*DBA.*', filteredphenodata$characteristics_ch1))
# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "2") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE6591"]][["Lung"]][["Male"]][["DBA2J"]] <- target$LogFC_table
name = "Mouse_GSE6591_Lung"
gender = "Male_DBA2J"
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

sexyphenodata = subset(filteredphenodata, subset = grepl('.*C57.*', filteredphenodata$characteristics_ch1))
# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CONTINUOUS (design matrix for continuous data)
# sexyphenodata$Age = sub("m$", "", sexyphenodata$Age) # delete nonnumeric characters
sexyphenodata$Age = as.numeric(as.character(sexyphenodata$Age))
design_matrix = model.matrix(~ sexyphenodata$Age)
rownames(design_matrix) = rownames(sexyphenodata)
colnames(design_matrix) <- c("Intercept", "Age")
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE6591"]][["Lung"]][["Male"]][["C57BL6J"]] <- target$LogFC_table
name = "Mouse_GSE6591_Lung"
gender = "Male_C57BL6J"
pdf(paste0("./plots/", name, "_", gender, "_processeddensity", ".pdf"))
print(target$Processed_density)
dev.off()
pdf(paste0("./plots/", name, "_", gender, "_topgeneexpression", ".pdf"))
print(target$Topgene_expression)
dev.off()

# save progress
save(logFClist, file = "logFClist.RData")

##### GSE6591 (END)

##### GSE74463 (BEGIN)

# get phenodata
gse = getGEO("GSE74463") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])

m = gregexpr("ENSMUST[[:digit:]]+", featuredata$gene_assignment, perl = TRUE)
ensembls = regmatches(featuredata$gene_assignment, m)
for (i in 1:length(ensembls)){
  if (length(ensembls[[i]]) == 0){
    ensembls[[i]] = NA
  } else if (length(ensembls[[i]]) > 1){
    ensembls[[i]] = NA
  }
}
featuredata$ENSEMBL_ID = unlist(ensembls)
featuredata[,"ID"] = as.character(featuredata[,"ID"])

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

rownames(filteredexprdata) = as.character(rownames(filteredexprdata))

logdata = filteredexprdata
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENSEMBL_ID"
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CONVERT TRANSCRIPTS
source("FUN.Ensembl_mouse_dictionary_create_for_trans.R") #if rat or human, use the corresponding function
dic = Ensembl_mouse_dictionary_create_for_trans(normdata)
# dic done here
source("FUN.Ensembl_to_entrez_for_microarray.R")
normdata = Ensembl_to_entrez_for_microarray(normdata, dic)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ PCA
exprforpca = data.frame(t(scale(t(normdata))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = filteredphenodata[, "Age"]
cluster_plot + geom_point(aes(color = color_pca))
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "14") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE74463"]][["Kidney"]][["Female"]] <- target$LogFC_table
name = "Mouse_GSE74463_Kidney"
gender = "Female"
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

##### GSE74463 (END)

##### GSE11291 (BEGIN)

# get phenodata
gse = getGEO("GSE11291") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata = subset(filteredphenodata, subset = grepl('.*control.*', filteredphenodata$source_name_ch1))
# takes only rows containing "sedentary" in the column "title"

# separate by tissue
# NONALPHANUMERIC
# (if everything is in one column separated with underscores or sth similar)
filteredphenodata = filteredphenodata %>% separate(source_name_ch1, c("Tissue", "Age", NA))

filteredphenodata1 = filteredphenodata

# analyze heart:
filteredphenodata = subset(filteredphenodata1, Tissue == "heart")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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

# discarding outliers: GSM284967 and GSM284975
filteredphenodata = filteredphenodata[-c(1,9),]

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
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

entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "5") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE11291"]][["Heart"]][["Male"]] <- target$LogFC_table
name = "Mouse_GSE11291_Heart"
gender = "Male"

# analyze muscle:
filteredphenodata = subset(filteredphenodata1, Tissue == "gastrocnemius")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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

# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# discard outliers: GSM284990 and GSM284992
filteredphenodata = filteredphenodata[-c(4,6),]

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
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

entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "5") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE11291"]][["Muscle"]][["Male"]] <- target$LogFC_table
name = "Mouse_GSE11291_Muscle"
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



# analyze brain:
filteredphenodata = subset(filteredphenodata1, Tissue == "neocortex")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "5") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE11291"]][["Brain"]][["Male"]] <- target$LogFC_table
name = "Mouse_GSE11291_Brain"
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

##### GSE11291 (END)

##### GSE34378 (BEGIN)

# get phenodata
gse = getGEO("GSE34378") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

filteredphenodata1 = filteredphenodata

#analyze liver:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "liver")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "geneENTREZ"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE34378"]][["Liver"]][["Female"]] <- target$LogFC_table
name = "Mouse_GSE34378_Liver"
gender = "Female"
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

#analyze kidney:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "kidney")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "geneENTREZ"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE34378"]][["Kidney"]][["Female"]] <- target$LogFC_table
name = "Mouse_GSE34378_Kidney"
gender = "Female"
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

#analyze lung:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "lung")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "geneENTREZ"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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

# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE34378"]][["Lung"]][["Female"]] <- target$LogFC_table
name = "Mouse_GSE34378_Lung"
gender = "Female"
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

#analyze brain:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "brain")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "geneENTREZ"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# discard outliers: GSM847865, GSM847851, and GSM847857
filteredphenodata = filteredphenodata[-c(3, 9, 17),]

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "geneENTREZ"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE34378"]][["Brain"]][["Female"]] <- target$LogFC_table
name = "Mouse_GSE34378_Brain"
gender = "Female"
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

##### GSE34378 (END)

##### GSE40645 (BEGIN)

# get phenodata
gse = getGEO("GSE40645") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

featuredata = fData(gse[[1]])

# set up target list element:
target = list()

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "Entrez_Gene_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# discard outlier: GSM998808
filteredphenodata = filteredphenodata[-12,]

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "Entrez_Gene_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

#analyze males:
sexyphenodata = subset(filteredphenodata, gender.ch1 == "male")

# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in years", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Human"]][["GSE40645"]][["Muscle"]][["Male"]] <- target$LogFC_table
name = "Human_GSE40645_Muscle"
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

#analyze females:
sexyphenodata = subset(filteredphenodata, gender.ch1 == "female")

# filter expression data
sexyexprdata = normdata[, rownames(sexyphenodata)]

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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in years", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Human"]][["GSE40645"]][["Muscle"]][["Female"]] <- target$LogFC_table
name = "Human_GSE40645_Muscle"
gender = "Female"
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

##### GSE40465 (END)

##### GSE12480 (BEGIN)

# get phenodata
gse = getGEO("GSE12480") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

filteredphenodata$Age = c("old", "old", "old", "old", "old", "old", "old", "old", "old", "old", "young", "young", "young", "young", "young", "young", "young", "young", "young", "young")

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
normdata <- exd %>% remove_rownames() %>% column_to_rownames(var = entrezcolname)
visualstack = stack(normdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# One gene has sd == 0 across samples
# Its IDs: 1443692_at, and entrez 433940
# how do you know that? code:
which(apply(normdata, 1, sd) == 0)
normdata[11896,]
which(featuredata$ENTREZ_GENE_ID == "433940")
rownames(featuredata)[25119]
which(rownames(filteredexprdata) == "1443692_at")

# take care of this motherfucker:
filteredexprdata = filteredexprdata[-27998,]

logdata = filteredexprdata
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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# discard outlier: GSM313564
filteredphenodata = filteredphenodata[-8,]

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# take care of this motherfucker:
filteredexprdata = filteredexprdata[-27998,]

logdata = filteredexprdata
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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "young") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
copps = as.data.frame(cbind(colnames(sexyexprdata), as.numeric(sexyexprdata[rownames(topgenes)[1],]), sexyphenodata$Age))
#copps$V3 = as.numeric(as.character(copps$V3))
copps$V2 = as.double(as.character(copps$V2))
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Mouse"]][["GSE12480"]][["Heart"]][["Male"]] <- target$LogFC_table
name = "Mouse_GSE12480_Heart"
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

##### GSE12480 (END)

##### GSE27625 (BEGIN)

# get phenodata
gse = getGEO("GSE27625") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata = subset(filteredphenodata, subset = grepl('.*control.*', filteredphenodata$diet.ch1))

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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

# discrad outlier: GSM684654
filteredphenodata = filteredphenodata[-13,]

# get filtered exprdata
filteredexprdata = data.frame(exprs(gse[[1]]))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

logdata = filteredexprdata
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

# LOOKUP (entrez is already there, just take the means)
# featuredata = fData(gse[[1]])
entrezcolname = "ENTREZ_GENE_ID"
featuredata = subset(featuredata, subset = grepl("^[0-9]+$", featuredata[, entrezcolname])) # get rid of multiple entrez per read ID
exd <- normdata
lookup = featuredata[, c('ID', entrezcolname)]
exd$ID = rownames(exd)
exd <- left_join(exd, lookup)
exd = na.omit(exd, cols=as.name(entrezcolname))
exd <- exd %>% group_by(!!as.name(entrezcolname)) %>% summarize_each(mean, matches("^GSM.*"))
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "3") # adjust for your age entries (here 3m
design_matrix = model.matrix(~ currentfactor) #                     is the youngest)
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
# copps$V2 = as.numeric(sub(",", ".", copps$V2))
copps = copps %>% arrange(V2)
ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))
# change months to weeks or years if needed in the next line:
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in months", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE27625"]][["Liver"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE27625_Liver"
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

##### GSE27625 (END)

##### GSE53960 (BEGIN)

# get phenodata
gse = getGEO("GSE53960") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# get expr counts
# code from exprfilemerger.R
setwd("./GSE53960")
filelist = list.files(pattern = ".*.txt")
datalist = lapply(filelist, function(x)read.table(x, header=T))
for (i in 1:length(datalist)){datalist[[i]] = datalist[[i]] %>% remove_rownames() %>% column_to_rownames(var = "AceVeiwGeneSymbol")}
filteredexprdata = datalist[[1]]
for (i in 2:length(datalist)){
  filteredexprdata = merge(filteredexprdata, datalist[[i]], by=0)
  filteredexprdata = filteredexprdata %>% column_to_rownames(var = "Row.names")
}
setwd("../")

filteredphenodata1 = filteredphenodata

filteredexprdata1 = filteredexprdata
filteredphenodata1$Age = as.character(filteredphenodata1$developmental.stage..week..ch1)

# set up target list element:
target = list()

# analyze brain:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "Brain")

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
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
# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# discard outlier: Brn_M_104_2
filteredphenodata = filteredphenodata[-30,]

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# analyze males:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "M")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Brain"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE53960_Brain"
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

# analyze females:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "F")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Brain"]][["Female"]] <- target$LogFC_table
name = "Rat_GSE53960_Brain"
gender = "Female"
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

# analyze heart:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "Heart")

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# analyze males:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "M")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Heart"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE53960_Heart"
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

# analyze females:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "F")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Heart"]][["Female"]] <- target$LogFC_table
name = "Rat_GSE53960_Heart"
gender = "Female"
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

# analyze kidney:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "Kidney")

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# analyze males:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "M")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Kidney"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE53960_Kidney"
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

# analyze females:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "F")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Kidney"]][["Female"]] <- target$LogFC_table
name = "Rat_GSE53960_Kidney"
gender = "Female"
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

# analyze liver:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "Liver")

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# analyze males:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "M")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Liver"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE53960_Liver"
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

# analyze females:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "F")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Liver"]][["Female"]] <- target$LogFC_table
name = "Rat_GSE53960_Liver"
gender = "Female"
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

# analyze lung:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "Lung")

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# analyze males:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "M")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Lung"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE53960_Lung"
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

# analyze females:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "F")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Lung"]][["Female"]] <- target$LogFC_table
name = "Rat_GSE53960_Lung"
gender = "Female"
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

# analyze muscle:
filteredphenodata = subset(filteredphenodata1, tissue.ch1 == "Muscle")

filteredexprdata = filteredexprdata1[, filteredphenodata$title]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# SMASH IT
# (if there is one, filter it out)
filteredexprdata$rowsum = rowSums(filteredexprdata > 20) # take only genes having more than
#                                                                10 reads...
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 10,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

source("FUN.Genesymbol_rat_dictionary_create.R")
dic = Genesymbol_rat_dictionary_create(filteredexprdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(filteredexprdata, dic)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# analyze males:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "M")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Muscle"]][["Male"]] <- target$LogFC_table
name = "Rat_GSE53960_Muscle"
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

# analyze females:
# (if you have a separate column)
sexyphenodata = subset(filteredphenodata, Sex.ch1 == "F")

# filter expression data
sexyexprdata = normdata[, sub("SEQC_", "", sexyphenodata$title)]

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
target[["Topgene_expression"]] <- ggplot(copps, aes(x = copps$V3, y = copps$V2, color = copps$V3)) + geom_point() + labs(colour = "age in weeks", x = "age", y = "Expression",title = paste("entrez ID", rownames(topgenes)[1], sep = " "))

# dump the data:
logFClist[["Rat"]][["GSE53960"]][["Muscle"]][["Female"]] <- target$LogFC_table
name = "Rat_GSE53960_Muscle"
gender = "Female"
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

##### GSE53960 (END)

##### GSE66715 (BEGIN)

#get exprdata
setwd("./GSE66715")
filelist = list.files(pattern = ".*.txt")
datalist = lapply(filelist, function(x)read.table(x, header=F))
for (i in 1:length(datalist)){
  rownames(datalist[[i]]) = datalist[[i]][,1]
  datalist[[i]] = datalist[[i]][3]
  colnames(datalist[[i]]) = c(filelist[i])
  }
filteredexprdata = datalist[[1]]
for (i in 2:length(datalist)){
  filteredexprdata = merge(filteredexprdata, datalist[[i]], by=0)
  filteredexprdata = filteredexprdata %>% column_to_rownames(var = "Row.names")
}
setwd("../")

# get phenodata
gse = getGEO("GSE66715") #here, type your gse id
filteredphenodata = data.frame(pData(gse[[1]]))

# set up target list element:
target = list()

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$age.ch1)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

colnames(filteredexprdata) = sub("_.*_.*_.*$", "", colnames(filteredexprdata))

filteredphenodata1 = filteredphenodata
# FILTER CONTROL REGEX (filter out noncontrol groups in phenodata)
filteredphenodata1 = subset(filteredphenodata1, subset = grepl('.*mRNA.*', filteredphenodata1$title))

# analyze liver:
filteredphenodata = subset(filteredphenodata1, source_name_ch1 == "liver")

filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

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
filteredexprdata1 = filteredexprdata[filteredexprdata$rowsum > 2,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata1, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CONVERT TRANSCRIPTS
source("FUN.Ensembl_rat_dictionary_create_for_trans.R") #if rat or human, use the corresponding function
dic = Ensembl_rat_dictionary_create_for_trans(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "6") # adjust for your age entries
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
logFClist[["Rat"]][["GSE66715"]][["Liver"]][["NoSex"]] <- target$LogFC_table
name = "Rat_GSE66715_Liver"
gender = "NoSex"
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

# analyze brain:
filteredphenodata = subset(filteredphenodata1, source_name_ch1 == "brain")

filteredexprdata = filteredexprdata[, rownames(filteredphenodata)]
filteredexprdata = na.omit(filteredexprdata)

# LOG2 (logarithmize if needed)
logdata = log2(filteredexprdata + 1)
# CHECK PEAK
# (first look if there is a peak at low values of expression)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Raw_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CONVERT TRANSCRIPTS
source("FUN.Ensembl_rat_dictionary_create_for_trans.R") #if rat or human, use the corresponding function
dic = Ensembl_rat_dictionary_create_for_trans(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

# Scale
sexyexprdata = data.frame(scale(sexyexprdata))
visualstack = stack(sexyexprdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
# save it:
target[["Processed_density"]] <- ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

# CATEGORIAL (design matrix for categorial data)
currentfactor <- factor(sexyphenodata$Age)
currentfactor = relevel(currentfactor, "6") # adjust for your age entries
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
logFClist[["Rat"]][["GSE66715"]][["Brain"]][["NoSex"]] <- target$LogFC_table
name = "Rat_GSE66715_Brain"
gender = "NoSex"
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

##### GSE66715 (END)

##### PRJNA281127 (BEGIN)

#first, obtain the annotation table using the tools on globe
#then, make it a phenodata table in R:
srapheno = read.csv("annotation_281127.csv")
srapheno = subset(srapheno, srapheno$LibraryStrategy == "RNA-Seq")
srapheno = srapheno %>% remove_rownames() %>% column_to_rownames(var = "Run")

#now, get the expression counts table:
setwd("./count_result_star")
source("../FUN.Download_raw_reads.R")
sraexpr = Download_raw_reads("./", srapheno)
setwd("../")

#after that, analyse the percentage of assigned reads...
setwd("./count_result_star")
assigneddata = data.frame()
totaldata = data.frame()
file.names <- dir("./", pattern =".count.summary")
for(i in 1:length(file.names)){
  countstable = read.table(file.names[i], sep = "")
  countstable = countstable[-1, ]
  assigneddata[str_remove(file.names[i], ".count.summary"), 1] = as.integer(as.character(countstable$V2[1]))
  totaldata[str_remove(file.names[i], ".count.summary"), 1] = sum(as.integer(as.character(countstable$V2)))}
colnames(assigneddata) = c("Assigned")
colnames(totaldata) = c("Total")
tableforplot = merge(assigneddata, totaldata, by=0)
tableforplot = tableforplot %>% column_to_rownames(var = "Row.names")
plot(tableforplot$Total, tableforplot$Assigned)
setwd("../")

# ...and filter out samples of poor quality:
tableforplot$Ratio = tableforplot$Assigned / tableforplot$Total
tableforplot = tableforplot %>% rownames_to_column("ID") %>% filter(Ratio >= 0.5) %>% column_to_rownames("ID")
ggplot(tableforplot, aes(x = Total, y = Assigned)) + geom_point()

# get filtered exprdata
filteredexprdata = data.frame(exprs(sraexpr))
# if you have the dataframe already, just name it filteredexprdata
filteredexprdata = filteredexprdata[, rownames(tableforplot)]
filteredexprdata = na.omit(filteredexprdata)
filteredexprdata1 = filteredexprdata

# get filtered phenodata
filteredphenodata = subset(srapheno, rownames(srapheno) %in% rownames(tableforplot))

# get age:
filteredphenodata$Age = as.character(c(12, 12, 12, 29, 29, 3, 3, 3, 3, 3, 12, 12, 29, 29, 3, 3, 3, 12, 12, 12, 29, 29, 29, 3, 29))

# get tissue:
filteredphenodata$Tissue = c("heart", "heart", "heart", "heart", "heart", "heart", "heart", "heart", "cerebellum", "cerebellum", "OB", "OB", "OB", "OB", "OB", "OB", "OB", "cerebellum", "cerebellum", "cerebellum", "cerebellum", "cerebellum", "cerebellum", "cerebellum", "OB")

filteredphenodata1 = filteredphenodata

# analyze heart:
filteredphenodata = subset(filteredphenodata1, Tissue == "heart")
filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]


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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 3,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Mouse"]][["PRJNA281127"]][["Heart"]][["Male"]] <- target$LogFC_table
name = "Mouse_PRJNA281127_Heart"
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

# analyze cerebellum:
filteredphenodata = subset(filteredphenodata1, Tissue == "cerebellum")
filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 3,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
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
# cluster_plot + geom_point(aes(color = color_pca)) + geom_text_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")
cluster_plot + geom_point(aes(color = color_pca)) + geom_label_repel(aes(label=rownames(cluster_values)),hjust="inward", vjust="inward")

# discard outliers: SRR5642536, SRR5642568
filteredphenodata = filteredphenodata[-c(1, 8),]
filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 2,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
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

# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Mouse"]][["PRJNA281127"]][["Cerebellum"]][["Male"]] <- target$LogFC_table
name = "Mouse_PRJNA281127_Cerebellum"
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

# analyze olfactory bulb:
filteredphenodata = subset(filteredphenodata1, Tissue == "OB")
filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 2,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Mouse"]][["PRJNA281127"]][["OB"]][["Male"]] <- target$LogFC_table
name = "Mouse_PRJNA281127_OB"
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

##### PRJNA281127 (END)

##### PRJNA516151 (BEGIN)

#first, obtain the annotation table using the tools on globe
#then, make it a phenodata table in R:
srapheno = read.csv("annotation_PRJNA516151.csv")
srapheno = subset(srapheno, srapheno$LibraryStrategy == "RNA-Seq")
srapheno = srapheno %>% remove_rownames() %>% column_to_rownames(var = "Run")

#now, get the expression counts table:
setwd("./ratsra")
source("../FUN.Download_raw_reads.R")
sraexpr = Download_raw_reads("./", srapheno)
setwd("../")

#after that, analyse the percentage of assigned reads...
setwd("./ratsra")
assigneddata = data.frame()
totaldata = data.frame()
file.names <- dir("./", pattern =".count.summary")
for(i in 1:length(file.names)){
  countstable = read.table(file.names[i], sep = "")
  countstable = countstable[-1, ]
  assigneddata[str_remove(file.names[i], ".count.summary"), 1] = as.integer(as.character(countstable$V2[1]))
  totaldata[str_remove(file.names[i], ".count.summary"), 1] = sum(as.integer(as.character(countstable$V2)))}
colnames(assigneddata) = c("Assigned")
colnames(totaldata) = c("Total")
tableforplot = merge(assigneddata, totaldata, by=0)
tableforplot = tableforplot %>% column_to_rownames(var = "Row.names")
plot(tableforplot$Total, tableforplot$Assigned)
setwd("../")

# ...and filter out samples of poor quality:
tableforplot$Ratio = tableforplot$Assigned / tableforplot$Total
#tableforplot = tableforplot %>% rownames_to_column("ID") %>% filter(Ratio >= 0.5) %>% column_to_rownames("ID")
#ggplot(tableforplot, aes(x = Total, y = Assigned)) + geom_point()

# get filtered exprdata
filteredexprdata = data.frame(exprs(sraexpr))
# if you have the dataframe already, just name it filteredexprdata
#filteredexprdata = filteredexprdata[, rownames(tableforplot)]
filteredexprdata = na.omit(filteredexprdata)
filteredexprdata1 = filteredexprdata

# get filtered phenodata
#filteredphenodata = subset(srapheno, rownames(srapheno) %in% rownames(tableforplot))

filteredphenodata = srapheno %>% separate(SampleName, c("Tissue", "Age", NA))

# GET AGE (make column "Age" with age integers)
filteredphenodata$Age = sub("^[^[:digit:]]*", "", filteredphenodata$Age)
filteredphenodata$Age = sub("[^[:digit:]]*$", "", filteredphenodata$Age)

filteredphenodata1 = filteredphenodata

# analyze liver:
filteredphenodata = subset(filteredphenodata1, Tissue == "Kidney")

filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 15,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
source("FUN.Ensembl_rat_dictionary_create.R")
dic = Ensembl_rat_dictionary_create(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Rat"]][["PRJNA516151"]][["Kidney"]][["Male"]] <- target$LogFC_table
name = "Rat_PRJNA516151_Kidney"
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

# analyze liver:
filteredphenodata = subset(filteredphenodata1, Tissue == "Liver")

filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 15,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
source("FUN.Ensembl_rat_dictionary_create.R")
dic = Ensembl_rat_dictionary_create(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Rat"]][["PRJNA516151"]][["Liver"]][["Male"]] <- target$LogFC_table
name = "Rat_PRJNA516151_Liver"
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

# analyze muscle:
filteredphenodata = subset(filteredphenodata1, Tissue == "Gastrocnemius")

filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 15,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
source("FUN.Ensembl_rat_dictionary_create.R")
dic = Ensembl_rat_dictionary_create(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Rat"]][["PRJNA516151"]][["Muscle"]][["Male"]] <- target$LogFC_table
name = "Rat_PRJNA516151_Muscle"
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

# analyze brain:
filteredphenodata = subset(filteredphenodata1, Tissue == "Hippocampus")

filteredexprdata = filteredexprdata1[, rownames(filteredphenodata)]

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
filteredexprdata2 = filteredexprdata[filteredexprdata$rowsum > 15,] # ...in at least
#                                                  one third of samples (here it is 3 samples)
filteredexprdata = subset(filteredexprdata2, select = -c(rowsum))
logdata = log2(filteredexprdata + 1)
visualstack = stack(logdata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

#$$$$$ Convert to Entrez (from Ensembl)(sum reads), and normalize with RLE
# CONVERT GENES (convert and take means)
source("FUN.Ensembl_rat_dictionary_create.R")
dic = Ensembl_rat_dictionary_create(filteredexprdata)
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
# save it:
target[["PCA"]] <- cluster_plot + geom_point(aes(color = color_pca))

# IF THERE IS ONLY ONE SEX:
sexyphenodata = filteredphenodata
sexyexprdata = normdata

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
logFClist[["Rat"]][["PRJNA516151"]][["Brain"]][["Male"]] <- target$LogFC_table
name = "Rat_PRJNA516151_Brain"
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





