# mixed effect model signature builder 2.0 (slow)

signature_builder = function(logFCmatrixregr, SEmatrixregr, name){
  goodgenes = c()
  signature = data.frame()
  genenumber = 0
  for (genename in rownames(logFCmatrixregr)){
    genenumber = genenumber + 1
    percentready = (genenumber/length(rownames(logFCmatrixregr))) * 100
    if (genenumber %% 1000 == 0){
      print(paste0("I'm on gene No. ", genenumber, " (", round(percentready, 2), "% done)"))
    }
    logFC = logFCmatrixregr[genename,]
    logFC = logFC[!is.na(logFC)]
    SE = SEmatrixregr[genename,]
    SE = SE[!is.na(SE)]
    sourcevec = as.factor(sourcedata[colnames(logFCmatrixregr)[!is.na(logFCmatrixregr[genename,])],])
    tissuevec = as.factor(sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\3", colnames(logFCmatrixregr)))
    speciesvec = as.factor(sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\1", colnames(logFCmatrixregr)))
    
    
    tryCatch(
      {
        if (name %in% c("Brain", "Muscle", "Liver")){
          mixedeffresfull = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec, ~ 1 | speciesvec))
        } else if(name %in% c("Human", "Rat", "Mouse")){
          mixedeffresfull = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec, ~ 1 | tissuevec))
        } else {
          mixedeffresfull = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec, ~ 1 | speciesvec, ~ 1 | tissuevec))
        }
        pvalsrob = c()
        logFCsrob = c()
        for (i in 1:length(colnames(logFCmatrixregr))){
          logFCmatrixoneout = logFCmatrixregr[,-i]
          SEmatrixoneout = SEmatrixregr[,-i]
          
          logFC = logFCmatrixoneout[genename,]
          logFC = logFC[!is.na(logFC)]
          SE = SEmatrixoneout[genename,]
          SE = SE[!is.na(SE)]
          sourcevec = as.factor(sourcedata[colnames(logFCmatrixoneout)[!is.na(logFCmatrixoneout[genename,])],])
          tissuevec = as.factor(sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\3", colnames(logFCmatrixoneout)))
          speciesvec = as.factor(sub("([^_]+)_([^_]+)_([^_]+)_.*", "\\1", colnames(logFCmatrixoneout)))
          if (name %in% c("Brain", "Muscle", "Liver")){
            mixedeffresrob = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec, ~ 1 | speciesvec))
          } else if(name %in% c("Human", "Rat", "Mouse")){
            mixedeffresrob = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec, ~ 1 | tissuevec))
          } else {
            mixedeffresrob = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec, ~ 1 | speciesvec, ~ 1 | tissuevec))
          }
          pvalsrob = c(pvalsrob, mixedeffresrob$pval)
          logFCsrob = c(logFCsrob, mixedeffresrob$b[1])
        }
        signature = rbind(signature, c(mixedeffresfull$b[1], mixedeffresfull$pval, logFCsrob[which.max(pvalsrob)], pvalsrob[which.max(pvalsrob)], logFCsrob[which.min(pvalsrob)], pvalsrob[which.min(pvalsrob)]))
        goodgenes = c(goodgenes, genename)
      },
      error=function(cond) {
        message("Fucked up")
        message("Here's the original error message:")
        message(cond)
      }
    )
    
    
    
  }
  #rownames(signature) = totalgenes[-which(totalgenes %in% badgenes)]
  rownames(signature) = goodgenes
  colnames(signature) = c("logFC", "pval", "logFC_robust", "pval_robust", "logFC_LOO", "pval_LOO")
  
  signature$adj_pval = p.adjust(signature$pval, method = "BH")
  signature$adj_pval_robust = p.adjust(signature$pval_robust, method = "BH")
  signature$adj_pval_LOO = p.adjust(signature$pval_LOO, method = "BH")
  return(signature)
}