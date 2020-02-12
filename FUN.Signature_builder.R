# mixed effect model signature builder

signature_builder = function(logFCmatrixregr){
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
    
    tryCatch(
      {
        mixedeffres = rma.mv(yi = logFC, V = SE^2, method = "REML", random = list(~ 1 | sourcevec))
        signature = rbind(signature, c(mixedeffres$b[1], mixedeffres$pval))
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
  colnames(signature) = c("logFC", "pval")
  
  signature$adj_pval = p.adjust(signature$pval, method = "BH")
  return(signature)
}