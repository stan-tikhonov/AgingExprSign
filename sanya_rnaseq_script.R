library(tidyverse)
sanya_expr = read.csv("RNA_combined_stranded_counts.csv")
sanya_expr = sanya_expr %>% remove_rownames %>% column_to_rownames(var="ID")
sanya_expr = subset(sanya_expr, select = -c(X))
sanya_pheno = read.table("Sanya_rnaseq_pheno.txt", sep = "\t")
sanya_pheno = subset(sanya_pheno, grepl("^kidney.*", V2))
#sanya_pheno = sanya_pheno[-5,] #delete liver 32m_1
sanya_expr = sanya_expr[, as.character(sanya_pheno$V1)]
sanya_expr = na.omit(sanya_expr)
filteredexprdata = sanya_expr
filteredphenodata = sanya_pheno
filteredphenodata = filteredphenodata[-18,] #32m_1 liver
filteredphenodata = filteredphenodata[-5,] #32m_1 kidney
filteredphenodata$V1 = trimws(as.character(filteredphenodata$V1))
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
#for outlier removal:
#cluster_plot + geom_point(aes(color = color_pca)) + geom_text(aes(label=rownames(cluster_values)),hjust=0, vjust=0)
#edgeR diff expression for categorial data
keknames = sub("_(.*)$", "", rownames(filteredexprdata1))
rownames(filteredexprdata1) = keknames
dic = Ensembl_mouse_dictionary_create(filteredexprdata1)
filteredexprdata1 = Ensembl_to_entrez(filteredexprdata1, dic)
currentfactor <- factor(filteredphenodata$Age)
currentfactor = relevel(currentfactor, "3m")
design_matrix = model.matrix(~ currentfactor)
rownames(design_matrix) = filteredphenodata$V1
colnames(design_matrix) <- c("(Intercept)", "20m")
edger_matrix = DGEList(filteredexprdata1)
edger_matrix = calcNormFactors(edger_matrix, method = "RLE")
edger_matrix = estimateDisp(edger_matrix, design_matrix, robust = TRUE)
fit = glmFit(edger_matrix, design_matrix)
result <- topTags(glmLRT(fit,coef="20m"),n=Inf,adjust.method = "BH")$table
#edgeR diff expression for continuous data
keknames = sub("_(.*)$", "", rownames(filteredexprdata1))
rownames(filteredexprdata1) = keknames
dic = Ensembl_mouse_dictionary_create(filteredexprdata1)
filteredexprdata1 = Ensembl_to_entrez(filteredexprdata1, dic)
filteredphenodata$Age = sub("m$", "", filteredphenodata$Age)
filteredphenodata$Age = as.integer(filteredphenodata$Age)
design_matrix = model.matrix(~ filteredphenodata$Age)
rownames(design_matrix) = filteredphenodata$V1
colnames(design_matrix) <- c("(Intercept)", "Age")
edger_matrix = DGEList(filteredexprdata1)
edger_matrix = calcNormFactors(edger_matrix, method = "RLE")
edger_matrix = estimateDisp(edger_matrix, design_matrix, robust = TRUE)
fit = glmFit(edger_matrix, design_matrix)
result <- topTags(glmLRT(fit,coef="Age"),n=Inf,adjust.method = "BH")$table
#correlation matrix and heatmap
alex_signatures = dget("Signatures_mouse_genes.R")
logFCmatrix = sanya_kidney_3ages["logFC"]
logFCmatrix = logFCmatrix %>% rename(sanya_kidney = logFC)
logFCmatrix = merge(logFCmatrix, sanya_liver_6ages["logFC"], by=0)
logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
logFCmatrix = logFCmatrix %>% rename(sanya_liver = logFC)
i = 1
for (el in alex_signatures$Species){
  logFCmatrix = merge(logFCmatrix, el["logFC"], by=0)
  logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
  colnames(logFCmatrix)[i + 2] = names(alex_signatures$Species)[i]
  i = i + 1
}
i = 1
for (el in alex_signatures$Interventions){
  logFCmatrix = merge(logFCmatrix, el["logFC"], by=0)
  logFCmatrix = logFCmatrix %>% column_to_rownames(var = "Row.names")
  colnames(logFCmatrix)[i + 5] = names(alex_signatures$Interventions)[i]
  i = i + 1
}
cormatrix <- round(cor(logFCmatrix, method = "spearman"),2)
library(reshape2)
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
upper_tri <- get_upper_tri(cormatrix)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()