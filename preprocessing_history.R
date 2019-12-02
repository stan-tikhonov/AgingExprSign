##### GSE5086 (BEGIN)

# get phenodata
gse = getGEO("GSE11845") #here, type your gse id
# set up target list element:
target = list()

# get featuredata!:
featuredata = fData(gse[[1]])
# fix it:
for (row in 1:nrow(problems(fData(gse[[1]])))) {
  featuredata[as.integer(problems(fData(gse[[1]]))[row, "row"]), as.character(problems(fData(gse[[1]]))[row, "col"])] = as.character(problems(fData(gse[[1]]))[row, "actual"])
}

load("filteredphenodataforGSE5086.RData") #Here I filtered out those who took medications or had any disease

#analyse male:
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





















