library(tidyverse)
sanya_expr = read.csv("RNA_combined_stranded_counts.csv")
sanya_expr = sanya_expr %>% remove_rownames %>% column_to_rownames(var="ID")
sanya_expr = subset(sanya_expr, select = -c(X))
sanya_pheno = read.table("Sanya_rnaseq_pheno.txt", sep = "\t")
sanya_pheno = subset(sanya_pheno, grepl("^kidney.*", V2))
#sanya_pheno = sanya_pheno[-5,] #delete liver 32m_1
sanya_expr = sanya_expr[, sanya_pheno$V1]
sanya_expr = na.omit(sanya_expr)
filteredexprdata = sanya_expr
filteredphenodata = sanya_pheno
filteredphenodata = filteredphenodata[-10,]
filteredphenodata = filteredphenodata[-c(5, 9),]
filteredphenodata$V1 = trimws(filteredphenodata$V1)
filteredexprdata = filteredexprdata[, filteredphenodata$V1]
filteredphenodata = filteredphenodata %>% separate(V2, c("Tissue", "Age", NA))
#logdata = log2(filteredexprdata + 1)
#filteredexprstack = stack(logdata)
#ggplot(filteredexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
filteredexprdata$rowsum = rowSums(filteredexprdata > 10)
filteredexprdata1 = filteredexprdata[filteredexprdata$rowsum > 3,]
filteredexprdata1 = subset(filteredexprdata1, select = -c(rowsum))
logdata = log2(filteredexprdata1 + 1)
filteredexprstack = stack(logdata)
ggplot(filteredexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
source("FUN.RLE_normalization.R")
library(edgeR)
normdata = RLE_normalization(logdata)
keknames = sub("_(.*)$", "", rownames(normdata))
rownames(normdata) = keknames
source("FUN.Ensembl_mouse_dictionary_create.R")
dic = Ensembl_mouse_dictionary_create(normdata)
source("FUN.Ensembl_to_entrez.R")
normdata = Ensembl_to_entrez(normdata, dic)
#plotting PCA
exprforpca = data.frame(t(scale(t(normdata))))
pcamodel = prcomp(t(exprforpca))
cluster_values = as.data.frame(pcamodel[['x']])
cluster_plot = ggplot(cluster_values, aes(x = PC1, y = PC2))
color_pca = filteredphenodata[, "Age"]
cluster_plot + geom_point(aes(color = color_pca))
#edgeR diff expression
keknames = sub("_(.*)$", "", rownames(filteredexprdata1))
rownames(filteredexprdata1) = keknames
dic = Ensembl_mouse_dictionary_create(filteredexprdata1)
filteredexprdata1 = Ensembl_to_entrez(filteredexprdata1, dic)
