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
filteredphenodata$V1 = trimws(filteredphenodata$V1)
filteredexprdata = filteredexprdata[, filteredphenodata$V1]
#logdata = log2(filteredexprdata + 1)
#filteredexprstack = stack(logdata)
#ggplot(filteredexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))
filteredexprdata$rowsum = rowSums(filteredexprdata > 10)
filteredexprdata1 = filteredexprdata[filteredexprdata$rowsum > 3,]
filteredexprdata1 = subset(filteredexprdata1, select = -c(rowsum))
logdata = log2(filteredexprdata1 + 1)
filteredexprstack = stack(logdata)
ggplot(filteredexprstack, aes=(x=values)) + geom_density(aes(x=values, group=ind, color=ind))

