library(reshape2)
library(psych)
library(BlandAltmanLeh)

##### correlation of correlation
corofcor = matrix(nrow = length(names(cormatrixsign)), ncol = length(names(cormatrixsign)))
colnames(corofcor) = names(cormatrixsign)
rownames(corofcor) = names(cormatrixsign)
for (thres1 in names(cormatrixsign)){
  # vectorize the upper triangles of matrices (diagonals excluded)
  vec1 = as.vector(as.matrix(cormatrixsign[[thres1]])[upper.tri(as.matrix(cormatrixsign[[thres1]]), diag = F)])
  for (thres2 in names(cormatrixsign)){
    vec2 = as.vector(as.matrix(cormatrixsign[[thres2]])[upper.tri(as.matrix(cormatrixsign[[thres2]]), diag = F)])
    corofcor[thres1, thres2] = cor(vec1, vec2, method = "spearman")
  }
}
# View the result
View(corofcor)

# Shapiro-Wilk test and QQplots:
for (thres in names(cormatrixsign)){
  print(shapiro.test(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F)]))
  qqnorm(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F)])
  qqline(cormatrixsign[[thres]][upper.tri(cormatrixsign[[thres]], diag = F)])
}
# thereby we proved that the correlations are not distributed normally and we cannot use pearson correlations with them

# plots of correlations of different methods against each other:
for (thres1 in names(cormatrixsign)){
  tempmatrix = matrix(0, length(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]), 2)
  tempmatrix[,1] = cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]
  for (thres2 in names(cormatrixsign)){
    if (thres1 == thres2) {
      break
    }
    tempmatrix[,2] = cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]
    tempmatrix = as.data.frame(tempmatrix)
    colnames(tempmatrix) = c(paste0("threshold_", thres1), paste0("threshold_", thres2))
    a = ggplot(tempmatrix, aes_string(x = colnames(tempmatrix)[1], y = colnames(tempmatrix)[2])) + geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + geom_abline(intercept = 0)
    print(a)
  }
}

##### chi-squared test, #const/#total, kappa coef and stacked barplot
# generate the 1, -1 and 0 tables
cortestsign = list()
cortestsign[["100"]] = data.frame()
cortestsign[["200"]] = data.frame()
cortestsign[["300"]] = data.frame()
cortestsign[["all"]] = data.frame()
for (thres in names(cormatrixsign)){
  for (colname in colnames(cormatrixsign[[thres]])){
    for (rowname in rownames(cormatrixsign[[thres]])){
      if (coradjpvalsign[[thres]][rowname, colname] < 0.1){
        if (cormatrixsign[[thres]][rowname, colname] > 0){
          cortestsign[[thres]][rowname, colname] = 1
        } else {
          cortestsign[[thres]][rowname, colname] = -1
        } 
      } else{
        cortestsign[[thres]][rowname, colname] = 0
      }
    }
  }
}

# stacked barplot showing how many ones and minus ones are there for each data point, and then names of bad data points: 
new.df<-melt(cortestsign[["200"]])
new.df<-new.df[complete.cases(new.df),]
new.df$Score<-factor(new.df$value, levels = c("1", "-1", "0"))
new.df = new.df %>% filter(Score != 0)
ggplot(new.df, aes(x=variable, fill=Score)) + geom_bar()
anothadf = cortestsign[["200"]]
anothadf$ones = rowSums(anothadf == 1)
anothadf$minusones = rowSums(anothadf == -1)
rownames(anothadf)[which(anothadf$ones <= anothadf$minusones)]

# run chi-squared
chisqtable = matrix(nrow = length(names(cortestsign)), ncol = length(names(cortestsign)))
colnames(chisqtable) = names(cortestsign)
rownames(chisqtable) = names(cortestsign)
for (thres1 in names(cortestsign)){
  # vectorize the upper triangles of matrices (diagonals excluded)
  vec1 = as.vector(as.matrix(cortestsign[[thres1]])[upper.tri(as.matrix(cortestsign[[thres1]]), diag = F)])
  for (thres2 in names(cortestsign)){
    vec2 = as.vector(as.matrix(cortestsign[[thres2]])[upper.tri(as.matrix(cortestsign[[thres2]]), diag = F)])
    chisqtable[thres1, thres2] = as.numeric(chisq.test(vec1, vec2)$p.value)
  }
}
# View the result
View(chisqtable)

# calculate #const/#total:
portiontable = matrix(nrow = length(names(cortestsign)), ncol = length(names(cortestsign)))
colnames(portiontable) = names(cortestsign)
rownames(portiontable) = names(cortestsign)
for (thres1 in names(cortestsign)){
  for (thres2 in names(cortestsign)){
    consistentones = sum(cortestsign[[thres1]][upper.tri(cortestsign[[thres1]], diag = F)] == cortestsign[[thres2]][upper.tri(cortestsign[[thres2]], diag = F)])
    total = length(cortestsign[[thres1]][upper.tri(cortestsign[[thres1]], diag = F)])
    portiontable[thres1, thres2] = consistentones / total
    portiontable[thres2, thres1] = consistentones / total
  }
}
# check out the result:
View(portiontable)

# calculate Cohen's kappa coef:
kappatable = matrix(nrow = length(names(cortestsign)), ncol = length(names(cortestsign)))
colnames(kappatable) = names(cortestsign)
rownames(kappatable) = names(cortestsign)
for (thres1 in names(cortestsign)){
  for (thres2 in names(cortestsign)){
    tempmatrix = matrix(0, length(cortestsign[[thres1]][upper.tri(cortestsign[[thres1]], diag = F)]), 2)
    tempmatrix[,1] = as.factor(cortestsign[[thres1]][upper.tri(cortestsign[[thres1]], diag = F)])
    tempmatrix[,2] = as.factor(cortestsign[[thres2]][upper.tri(cortestsign[[thres2]], diag = F)])
    #tempmatrix = tempmatrix + 2
    kappatable[thres1, thres2] = cohen.kappa(tempmatrix)$kappa
    kappatable[thres2, thres1] = paste(as.vector(cohen.kappa(tempmatrix)$confid["unweighted kappa",]), collapse = "<")
  }
}
# check out the result:
View(kappatable)

# Bland-Altman plots for signed value of correlation:
for (thres1 in names(cormatrixsign)){
  for (thres2 in names(cormatrixsign)){
    if (thres1 == thres2){
      break
    }
    meandif = round(mean(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)] - cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]), 5)
    a = as.numeric(wilcox.test(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)], cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)])$p.value)
    bland.altman.plot(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)], cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)], main=paste0("Bland Altman Plot of ", thres1, " against ", thres2, " (Mann-Wintey U test p-value: ", a, "; mean difference: ", meandif, ")"), xlab="Means", ylab=paste0("Differences (", thres1, " - ", thres2, ")"))
  }
}

# Bland-Altman plots for absolute value of correlation:
for (thres1 in names(cormatrixsign)){
  for (thres2 in names(cormatrixsign)){
    if (thres1 == thres2){
      break
    }
    meandif = round(mean(abs(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]) - abs(cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)])), 5)
    a = as.numeric(wilcox.test(abs(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]), abs(cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]))$p.value)
    bland.altman.plot(abs(cormatrixsign[[thres1]][upper.tri(cormatrixsign[[thres1]], diag = F)]), abs(cormatrixsign[[thres2]][upper.tri(cormatrixsign[[thres2]], diag = F)]), main=paste0("Bland Altman Plot of ", thres1, " against ", thres2, " (Mann-Wintey U test p-value: ", a, "; mean difference: ", meandif, ")"), xlab="Means", ylab=paste0("Differences (", thres1, " - ", thres2, ")"))
  }
}

##### datasets consistency:

# distribution of correlation coefficients:



