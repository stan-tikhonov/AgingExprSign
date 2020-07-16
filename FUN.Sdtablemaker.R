# set up the Brown-Forsythe test function
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

sdtablemaker = function(normdata, filteredphenodata){
  # Get SDs for all of the genes
  tnormdata = t(normdata)
  tnormdata = as.data.frame(tnormdata)
  tnormdata$Age = filteredphenodata[rownames(tnormdata),]$Age
  tnormdata = tnormdata %>% group_by(Age) %>% summarise_all(sd)
  tnormdata = tnormdata %>% column_to_rownames("Age")
  sddata = as.data.frame(t(tnormdata))
  maxage = as.character(max(as.numeric(colnames(sddata))))
  minage = as.character(min(as.numeric(colnames(sddata))))
  
  # Perform BF test between the oldest and the youngest; add sd ratios to the table
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
  return(sddata)
}