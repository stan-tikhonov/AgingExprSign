##### correlation of correlation
corofcor = matrix(nrow = length(names(cormatrixsign)), ncol = length(names(cormatrixsign)))
colnames(corofcor) = names(cormatrixsign)
rownames(corofcor) = names(cormatrixsign)
for (thres1 in names(cormatrixsign)){
  # vectorize the upper triangles of matrices (diagonals excluded)
  vec1 = as.vector(as.matrix(cormatrixsign[[thres1]])[upper.tri(as.matrix(cormatrixsign[[thres1]]), diag = F)])
  for (thres2 in names(cormatrixsign)){
    vec2 = as.vector(as.matrix(cormatrixsign[[thres2]])[upper.tri(as.matrix(cormatrixsign[[thres2]]), diag = F)])
    corofcor[thres1, thres2] = cor(vec1, vec2)
  }
}
# View the result
View(corofcor)

##### chi-squared test and #const/#total
# generate the 1, -1 and 0 tables
cortestsign = list()
cortestsign[["100"]] = data.frame()
cortestsign[["200"]] = data.frame()
cortestsign[["300"]] = data.frame()
cortestsign[["all"]] = data.frame()
for (thres in names(cormatrixsign)){
  for (colname in colnames(cormatrixsign[[thres]])){
    for (rowname in rownames(cormatrixsign[[thres]])){
      if (corpvalsign[[thres]][rowname, colname] < 0.05){
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

