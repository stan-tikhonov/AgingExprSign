# Deming regression minimizer

deming_minimizer = function(logFCmatrixregr, totalrownamematrix){
  fn = function(k_no_first){
    k = c()
    k[1] = 1
    k[2:length(colnames(logFCmatrixregr))] = k_no_first
    #k[2:10] = k_no_first
    res = 0
    for (i in 1:(length(colnames(logFCmatrixregr))-1)){
      #for (i in 1:9){
      namei = colnames(logFCmatrixregr)[i]
      for (j in (i + 1):length(colnames(logFCmatrixregr))){
        #for (j in (i + 1):10){
        namej = colnames(logFCmatrixregr)[j]
        if (cortestsign[namei, namej] != 1){
          next
        }
        totalrownames = totalrownamematrix[[namei]][[namej]]
        ai = logFCmatrixregr[totalrownames, namei]
        aj = logFCmatrixregr[totalrownames, namej]
        res = res + sum(
          (((aj - (k[j]/k[i])*ai)^2)*((ai - (k[i]/k[j])*aj)^2))/
            (((aj - (k[j]/k[i])*ai)^2)+((ai - (k[i]/k[j])*aj)^2)))/length(totalrownames)
      }
    }
    return(res)
  }
  
  kvec = rnorm(length(colnames(logFCmatrixregr)) - 1, 1, 1)
  #kvec = rnorm(9, 1, 1)
  ptm <- proc.time()
  #optimized = optim(kvec, fn)
  #optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B", control = list(factr = 1e3))
  optimized = optim(kvec, fn, lower = 0.01, upper = 100, method = "L-BFGS-B")
  proc.time() - ptm
  
  kres = c(1, optimized$par)
  minimum = optimized$value
  bigres = list(kres, minimum)
  names(bigres) = c("coefs", "minimum")
  return(bigres)
}