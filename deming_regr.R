d = list()
a = rnorm(1000, 0, 1)
b = 2*a+rnorm(1000, 0, 1)
b = b
b <- b/sd(b)
c = 3*a + rnorm(1000, 0, 1)
d = list(a, b, c)
d = list(a, b)
fn = function(k){
  res = 0
  for (i in 1:length(d)){
    for (j in i:length(d)){
      if (i != j){
        res = res + sum((((d[[j]] - (k[j]/k[i])*d[[i]])^2)*((d[[i]] - (k[i]/k[j])*d[[j]])^2))/(((d[[j]] - (k[j]/k[i])*d[[i]])^2)+((d[[i]] - (k[i]/k[j])*d[[j]])^2)))/1000
      }  
    }
  }
  return(res)
}
fn1 = function(k){
  sum((((a - k*b)^2)+((b-(1/k)*a)^2)))
}

plot2 = optim(rnorm(length(d),0,1), fn)
dem_coefs <- c()
for(i in 1:(length(d)-1)){
  for(j in (i+1):length(d)){
    plot1 = deming(d[[j]] ~ d[[i]] - 1)
    plot(d[[i]], d[[j]], xlim = c(-1, 1), ylim = c(-1, 1))
    plot(d[[i]], d[[j]])
    abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
    abline(0, plot2$par[j]/plot2$par[i], col = "red", lwd = 2)
    abline(v = 0)
    abline(h = 0)
    plot3 = lm(d[[j]] ~ d[[i]] - 1)
    plot27 = lm(d[[i]] ~ d[[j]] - 1)
    #abline(0, plot3$coef, col = "green", lwd = 2)
    #abline(0, 1/plot27$coef, col = "darkgreen", lwd = 2)
    dem_coefs <- c(dem_coefs,plot1$coefficients[2])
  }
}
plot2$par/plot2$par[1]
plot2$par[3]/plot2$par[2]
dem_coefs[2]/dem_coefs[1]

dem_coefs/plot2$par[2]*plot2$par[1]
plot1 = deming(d[[1]] ~ b - 1)
plot2 = optim(c(1, 1, 1), fn)
plot3 = lm(a ~ b - 1)
plot27 = lm(b ~ a - 1)
plot4 = optim(c(1, 1), fn1)
plot(d[[2]], d[[1]])
abline()
plot(b, a,xlim=c(-3,3),ylim=c(-3,3))
abline(0, plot1$coefficients[2], col = "blue", lwd = 2)
abline(0, plot2$par, col = "red", lwd = 2)
#abline(0, plot3$coef, col = "green", lwd = 2)
abline(v = 0)
abline(h = 0)
abline(0, plot4$par, col = "darkred", lwd = 2)
abline(0, 1/plot27$coef, col = "darkgreen", lwd = 2)