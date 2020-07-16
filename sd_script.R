
library(tidyverse)
library(preprocessCore)
library(GEOquery)
library(limma)
library(ggrepel)
library(edgeR)
library(onewaytests)
library(lawstat)
source("FUN.Sdtablemaker.R")

sdlist = list()

brown_forsythe_pval = function(y, group){
  n <- length(y)
  x.levels <- levels(factor(group))
  y.vars <- y.means <- m <- y.n <- NULL
  y.mean = mean(y)
  for (i in x.levels) {
    y.vars[i] <- var(y[group == i])
    y.means[i] <- mean(y[group == i])
    y.n[i] <- length(y[group == i])
  }
  for (j in x.levels) {
    m[j] <- (1 - y.n[j]/n) * (y.vars[j])/sum((1 - y.n/n) * 
                                               (y.vars))
  }
  SSb = sum(y.n * ((y.means - y.mean)^2))
  denom = sum((1 - y.n/n) * (y.vars))
  Ftest = SSb/denom
  df1 = length(x.levels) - 1
  df2 = 1/(sum(m^2/(y.n - 1)))
  p.value = pf(Ftest, df1, df2, lower.tail = F)
  return(p.value)
}

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

# FUNC

sdtable = sdtablemaker(normdata, filteredphenodata)

sdlist[["Mouse_GSE3150_Liver_Male"]] = sdtable

# CODE

# calculate SDs and perform Levene's test

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
maxage = as.character(max(as.numeric(colnames(sddata))))
minage = as.character(min(as.numeric(colnames(sddata))))

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
if (length(unique(tnormdata$Age)) > 2){
  tnormdata = tnormdata %>% filter((Age %in% c(minage, maxage)))
}
tnormdata$Age = as.factor(tnormdata$Age)
pvalues = c()

for (rowname in rownames(sddata)){
  #out = brown_forsythe_pval(tnormdata[[rowname]], tnormdata[["Age"]])
  out = levene.test(tnormdata[, rowname], tnormdata[,"Age"], location = "mean")
  pvalues = c(pvalues, out$p.value)
}
sddata$pvalue = pvalues
sddata$adjpval = p.adjust(sddata$pvalue, method = "BH")
sddata$ratio = sddata[, maxage] / sddata[, minage]

#proof of principle:
temp = sddata %>% filter(adjpval < 0.05) %>% mutate(logratio = log(ratio))
temp = temp %>% group_by(logratio) %>% summarise(count = n())
#print(ggplot(temp %>% filter(logratio > -0.005 & logratio < 0.005), aes(x = logratio, y = count)) + geom_point())
#print(ggplot(temp, aes(x = logratio)) + geom_point(aes(x = logratio, y = count())))
print(ggplot(temp, aes(x = logratio)) + geom_density() + theme_minimal())

# perform the "mean" approach

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(mean)
tnormdata = tnormdata %>% column_to_rownames("Age")
meandata = as.data.frame(t(tnormdata))

metrictable = as.data.frame(cbind(rownames(filteredphenodata), filteredphenodata$Age))
metrictable = metrictable %>% column_to_rownames("V1")
metrictable = cbind(metrictable, rep("kok", length(rownames(metrictable))), rep("kok", length(rownames(metrictable))))
colnames(metrictable) = c("Age", "Euclidean", "Pearson")
metrictable$Euclidean = as.character(metrictable$Euclidean)
metrictable$Pearson = as.character(metrictable$Pearson)

for(object in rownames(metrictable)){
  age = as.character(filteredphenodata[object, "Age"])
  metrictable[object, "Pearson"] = cor(normdata[, object], meandata[, age], method = "pearson")
  metrictable[object, "Euclidean"] = dist(rbind(normdata[, object], meandata[, age]))
}
metrictable$Euclidean = as.numeric(metrictable$Euclidean)
metrictable$Pearson = as.numeric(metrictable$Pearson)

print(ggplot(metrictable, aes(x = Age, y = Euclidean)) + geom_boxplot() + theme_minimal())
print(ggplot(metrictable, aes(x = Age, y = Pearson)) + geom_boxplot() + theme_minimal())





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

#discard an outlier:
filteredphenodata = filteredphenodata[-13,]

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

sdtable = sdtablemaker(normdata, filteredphenodata)

#proof of principle:
temp = sdtable %>% filter(adjpval < 0.05) %>% mutate(logratio = log(ratio))
temp = temp %>% group_by(logratio) %>% summarise(count = n())
print(ggplot(temp %>% filter(logratio > -0.005 & logratio < 0.005), aes(x = logratio, y = count)) + geom_point())
#print(ggplot(temp, aes(x = logratio)) + geom_point(aes(x = logratio, y = count())))

sdlist[["Rat_GSE27625_Liver_Male"]] = sdtable

# calculate SDs and perform Levene's test

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
maxage = as.character(max(as.numeric(colnames(sddata))))
minage = as.character(min(as.numeric(colnames(sddata))))

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
if (length(unique(tnormdata$Age)) > 2){
  tnormdata = tnormdata %>% filter((Age %in% c(minage, maxage)))
}
tnormdata$Age = as.factor(tnormdata$Age)
pvalues = c()

for (rowname in rownames(sddata)){
  #out = brown_forsythe_pval(tnormdata[[rowname]], tnormdata[["Age"]])
  out = levene.test(tnormdata[, rowname], tnormdata[,"Age"], location = "mean")
  pvalues = c(pvalues, out$p.value)
}
sddata$pvalue = pvalues
sddata$adjpval = p.adjust(sddata$pvalue, method = "BH")
sddata$ratio = sddata[, maxage] / sddata[, minage]

#proof of principle:
temp = sddata %>% filter(adjpval < 0.05) %>% mutate(logratio = log(ratio))
temp = temp %>% group_by(logratio) %>% summarise(count = n())
#print(ggplot(temp %>% filter(logratio > -0.005 & logratio < 0.005), aes(x = logratio, y = count)) + geom_point())
#print(ggplot(temp, aes(x = logratio)) + geom_point(aes(x = logratio, y = count())))
print(ggplot(temp, aes(x = logratio)) + geom_density() + theme_minimal())

# perform the "mean" approach

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(mean)
tnormdata = tnormdata %>% column_to_rownames("Age")
meandata = as.data.frame(t(tnormdata))

metrictable = as.data.frame(cbind(rownames(filteredphenodata), filteredphenodata$Age))
metrictable = metrictable %>% column_to_rownames("V1")
metrictable = cbind(metrictable, rep("kok", length(rownames(metrictable))), rep("kok", length(rownames(metrictable))))
colnames(metrictable) = c("Age", "Euclidean", "Pearson")
metrictable$Euclidean = as.character(metrictable$Euclidean)
metrictable$Pearson = as.character(metrictable$Pearson)

for(object in rownames(metrictable)){
  age = as.character(filteredphenodata[object, "Age"])
  metrictable[object, "Pearson"] = cor(normdata[, object], meandata[, age], method = "pearson")
  metrictable[object, "Euclidean"] = dist(rbind(normdata[, object], meandata[, age]))
}
metrictable$Euclidean = as.numeric(metrictable$Euclidean)
metrictable$Pearson = as.numeric(metrictable$Pearson)

print(ggplot(metrictable, aes(x = Age, y = Euclidean)) + geom_boxplot() + theme_minimal())
print(ggplot(metrictable, aes(x = Age, y = Pearson)) + geom_boxplot() + theme_minimal())

# general SD distribution
tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
visualstack = stack(sddata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + theme_minimal()

##### GSE27625 (END)

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

# general SD distribution

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
visualstack = stack(sddata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + theme_minimal()

# calculate SDs and perform Levene's test

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
maxage = as.character(max(as.numeric(colnames(sddata))))
minage = as.character(min(as.numeric(colnames(sddata))))

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
if (length(unique(tnormdata$Age)) > 2){
  tnormdata = tnormdata %>% filter((Age %in% c(minage, maxage)))
}
tnormdata$Age = as.factor(tnormdata$Age)
pvalues = c()

for (rowname in rownames(sddata)){
  #out = brown_forsythe_pval(tnormdata[[rowname]], tnormdata[["Age"]])
  out = levene.test(tnormdata[, rowname], tnormdata[,"Age"], location = "mean")
  pvalues = c(pvalues, out$p.value)
}
sddata$pvalue = pvalues
sddata$adjpval = p.adjust(sddata$pvalue, method = "BH")
sddata$ratio = sddata[, maxage] / sddata[, minage]

print(paste("Number of significant genes:", as.character(sum(sddata$adjpval < 0.05))))

#proof of principle:
temp = sddata %>% filter(adjpval < 0.05) %>% mutate(logratio = log(ratio))
temp = temp %>% group_by(logratio) %>% summarise(count = n())
#print(ggplot(temp %>% filter(logratio > -0.005 & logratio < 0.005), aes(x = logratio, y = count)) + geom_point())
#print(ggplot(temp, aes(x = logratio)) + geom_point(aes(x = logratio, y = count())))
print(ggplot(temp, aes(x = logratio)) + geom_density() + theme_minimal())

# perform the "mean" approach

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(mean)
tnormdata = tnormdata %>% column_to_rownames("Age")
meandata = as.data.frame(t(tnormdata))

metrictable = as.data.frame(cbind(rownames(filteredphenodata), filteredphenodata$Age))
metrictable = metrictable %>% column_to_rownames("V1")
metrictable = cbind(metrictable, rep("kok", length(rownames(metrictable))), rep("kok", length(rownames(metrictable))))
colnames(metrictable) = c("Age", "Euclidean", "Pearson")
metrictable$Euclidean = as.character(metrictable$Euclidean)
metrictable$Pearson = as.character(metrictable$Pearson)

for(object in rownames(metrictable)){
  age = as.character(filteredphenodata[object, "Age"])
  metrictable[object, "Pearson"] = cor(normdata[, object], meandata[, age], method = "pearson")
  metrictable[object, "Euclidean"] = dist(rbind(normdata[, object], meandata[, age]))
}
metrictable$Euclidean = as.numeric(metrictable$Euclidean)
metrictable$Pearson = as.numeric(metrictable$Pearson)

print(ggplot(metrictable, aes(x = Age, y = Euclidean)) + geom_boxplot() + theme_minimal())
print(ggplot(metrictable, aes(x = Age, y = Pearson)) + geom_boxplot() + theme_minimal())





# analyze muscle:
filteredphenodata = subset(filteredphenodata1, Tissue == "gastrocnemius")

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

# general SD distribution

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
visualstack = stack(sddata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + theme_minimal()

# calculate SDs and perform Levene's test

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
maxage = as.character(max(as.numeric(colnames(sddata))))
minage = as.character(min(as.numeric(colnames(sddata))))

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
if (length(unique(tnormdata$Age)) > 2){
  tnormdata = tnormdata %>% filter((Age %in% c(minage, maxage)))
}
tnormdata$Age = as.factor(tnormdata$Age)
pvalues = c()

for (rowname in rownames(sddata)){
  #out = brown_forsythe_pval(tnormdata[[rowname]], tnormdata[["Age"]])
  out = levene.test(tnormdata[, rowname], tnormdata[,"Age"], location = "mean")
  pvalues = c(pvalues, out$p.value)
}
sddata$pvalue = pvalues
sddata$adjpval = p.adjust(sddata$pvalue, method = "BH")
sddata$ratio = sddata[, maxage] / sddata[, minage]

print(paste("Number of significant genes:", as.character(sum(sddata$adjpval < 0.05))))

#proof of principle:
temp = sddata %>% filter(adjpval < 0.05) %>% mutate(logratio = log(ratio))
temp = temp %>% group_by(logratio) %>% summarise(count = n())
#print(ggplot(temp %>% filter(logratio > -0.005 & logratio < 0.005), aes(x = logratio, y = count)) + geom_point())
#print(ggplot(temp, aes(x = logratio)) + geom_point(aes(x = logratio, y = count())))
print(ggplot(temp, aes(x = logratio)) + geom_density() + theme_minimal())

# perform the "mean" approach

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(mean)
tnormdata = tnormdata %>% column_to_rownames("Age")
meandata = as.data.frame(t(tnormdata))

metrictable = as.data.frame(cbind(rownames(filteredphenodata), filteredphenodata$Age))
metrictable = metrictable %>% column_to_rownames("V1")
metrictable = cbind(metrictable, rep("kok", length(rownames(metrictable))), rep("kok", length(rownames(metrictable))))
colnames(metrictable) = c("Age", "Euclidean", "Pearson")
metrictable$Euclidean = as.character(metrictable$Euclidean)
metrictable$Pearson = as.character(metrictable$Pearson)

for(object in rownames(metrictable)){
  age = as.character(filteredphenodata[object, "Age"])
  metrictable[object, "Pearson"] = cor(normdata[, object], meandata[, age], method = "pearson")
  metrictable[object, "Euclidean"] = dist(rbind(normdata[, object], meandata[, age]))
}
metrictable$Euclidean = as.numeric(metrictable$Euclidean)
metrictable$Pearson = as.numeric(metrictable$Pearson)

print(ggplot(metrictable, aes(x = Age, y = Euclidean)) + geom_boxplot() + theme_minimal())
print(ggplot(metrictable, aes(x = Age, y = Pearson)) + geom_boxplot() + theme_minimal())







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

# general SD distribution

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
visualstack = stack(sddata)
ggplot(visualstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind)) + theme_minimal()

# calculate SDs and perform Levene's test

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
tnormdata = tnormdata %>% column_to_rownames("Age")
sddata = as.data.frame(t(tnormdata))
maxage = as.character(max(as.numeric(colnames(sddata))))
minage = as.character(min(as.numeric(colnames(sddata))))

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
if (length(unique(tnormdata$Age)) > 2){
  tnormdata = tnormdata %>% filter((Age %in% c(minage, maxage)))
}
tnormdata$Age = as.factor(tnormdata$Age)
pvalues = c()

for (rowname in rownames(sddata)){
  #out = brown_forsythe_pval(tnormdata[[rowname]], tnormdata[["Age"]])
  out = levene.test(tnormdata[, rowname], tnormdata[,"Age"], location = "mean")
  pvalues = c(pvalues, out$p.value)
}
sddata$pvalue = pvalues
sddata$adjpval = p.adjust(sddata$pvalue, method = "BH")
sddata$ratio = sddata[, maxage] / sddata[, minage]

print(paste("Number of significant genes:", as.character(sum(sddata$adjpval < 0.05))))

#proof of principle:
temp = sddata %>% filter(adjpval < 0.05) %>% mutate(logratio = log(ratio))
temp = temp %>% group_by(logratio) %>% summarise(count = n())
#print(ggplot(temp %>% filter(logratio > -0.005 & logratio < 0.005), aes(x = logratio, y = count)) + geom_point())
#print(ggplot(temp, aes(x = logratio)) + geom_point(aes(x = logratio, y = count())))
print(ggplot(temp, aes(x = logratio)) + geom_density() + theme_minimal())

# perform the "mean" approach

tnormdata = t(normdata)
tnormdata = as.data.frame(tnormdata)
tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(mean)
tnormdata = tnormdata %>% column_to_rownames("Age")
meandata = as.data.frame(t(tnormdata))

metrictable = as.data.frame(cbind(rownames(filteredphenodata), filteredphenodata$Age))
metrictable = metrictable %>% column_to_rownames("V1")
metrictable = cbind(metrictable, rep("kok", length(rownames(metrictable))), rep("kok", length(rownames(metrictable))))
colnames(metrictable) = c("Age", "Euclidean", "Pearson")
metrictable$Euclidean = as.character(metrictable$Euclidean)
metrictable$Pearson = as.character(metrictable$Pearson)

for(object in rownames(metrictable)){
  age = as.character(filteredphenodata[object, "Age"])
  metrictable[object, "Pearson"] = cor(normdata[, object], meandata[, age], method = "pearson")
  metrictable[object, "Euclidean"] = dist(rbind(normdata[, object], meandata[, age]))
}
metrictable$Euclidean = as.numeric(metrictable$Euclidean)
metrictable$Pearson = as.numeric(metrictable$Pearson)
metrictable$Age = as.numeric(as.character(metrictable$Age))

metrictable = metrictable %>% arrange(Age)

print(ggplot(metrictable, aes(x = Age, y = Euclidean, group = Age)) + geom_boxplot() + theme_minimal())
print(ggplot(metrictable, aes(x = Age, y = Pearson, group = Age)) + geom_boxplot() + theme_minimal())

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



