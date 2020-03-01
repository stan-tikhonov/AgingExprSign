# make cortestsign:
cormatricesmaker = function(logFCunlisted, cormethod, signifgenesthres, adjmethod = "BH", corthres = 0.1){
  corpvalsign = data.frame()
  cormatrixsign = data.frame()
  for (i in 1:length(logFCunlisted)){
    for (j in i:length(logFCunlisted)){
      topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
      topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topA = topA %>% column_to_rownames(var = "row.names")
      topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
      topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topB = topB %>% column_to_rownames(var = "row.names")
      totalrownames = union(rownames(topA), rownames(topB))
      cormatrixsign[names(logFCunlisted)[i], names(logFCunlisted)[j]] = cor(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
      cormatrixsign[names(logFCunlisted)[j], names(logFCunlisted)[i]] = cor(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod, use = "complete.obs")
      corpvalsign[names(logFCunlisted)[i], names(logFCunlisted)[j]] = as.numeric(cor.test(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
      corpvalsign[names(logFCunlisted)[j], names(logFCunlisted)[i]] = as.numeric(cor.test(logFCunlisted[[i]][totalrownames,]$logFC, logFCunlisted[[j]][totalrownames,]$logFC, method = cormethod)$p.value)
      #    mergedmatrix = logFCunlisted[[i]]["logFC"]
      #    mergedmatrix = mergedmatrix %>% dplyr::rename(logFCi = logFC)
      #    mergedmatrix = merge(mergedmatrix, logFCunlisted[[j]]["logFC"], by=0, all=TRUE)
      #    mergedmatrix = mergedmatrix %>% column_to_rownames("Row.names")
      #    cormatrixdenoised[names(logFCunlisted)[i], names(logFCunlisted)[j]] = round(cor(mergedmatrix[union(rownames(topA), rownames(topB)),], method = "spearman", use = "complete.obs"),2)[2,1]
    }
  }
  coradjpvalsign = data.frame()
  vec = as.vector(corpvalsign[upper.tri(corpvalsign, diag = F)])
  vec = p.adjust(vec, method = adjmethod)
  tempmatrix = matrix(0, length(logFCunlisted), length(logFCunlisted))
  tempmatrix[upper.tri(tempmatrix, diag = F)] = vec
  tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
  #tempmatrix[lower.tri(tempmatrix, diag = F)] = vec
  #tempmatrix = t(tempmatrix)
  #tempmatrix[lower.tri(tempmatrix, diag = F)] = t(tempmatrix)[lower.tri(t(tempmatrix), diag = F)]
  coradjpvalsign = tempmatrix
  colnames(coradjpvalsign) = colnames(corpvalsign)
  rownames(coradjpvalsign) = rownames(corpvalsign)
  cortestsign = data.frame()
  for (colname in colnames(cormatrixsign)){
    for (rowname in rownames(cormatrixsign)){
      if (coradjpvalsign[rowname, colname] < 0.05){
        if (cormatrixsign[rowname, colname] > corthres){
          cortestsign[rowname, colname] = 1
        } else if (cormatrixsign[rowname, colname] < -corthres){
          cortestsign[rowname, colname] = -1
        } else {
          cortestsign[rowname, colname] = 0
        }
      } else{
        cortestsign[rowname, colname] = 0
      }
    }
  }
  res = list(cortestsign, cormatrixsign, corpvalsign, coradjpvalsign)
  names(res) = c("cortestsign", "cormatrixsign", "corpvalsign", "coradjpvalsign")
  return(res)
}