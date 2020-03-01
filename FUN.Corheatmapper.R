# corheatmapper

library(R.utils)
library(reshape2)

# functions for heatmaps:
reorder_cormat <- function(cormat,method="complete"){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd,method = method)
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

corheatmapper = function(cormatrix, cormethod, corlimit = 1){
  cormatrix_2 <- reorder_cormat(cormatrix, method="average")
  cormatrix_2 = apply(cormatrix_2, 2, rev)
  upper_tri <- get_upper_tri(cormatrix_2)
  melted_cormat <- melt(cormatrix_2, na.rm = TRUE)
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue4", high = "red4", mid = "white", 
                         midpoint = 0, limit = c(-corlimit,corlimit), space = "Lab", 
                         name=paste0(capitalize(cormethod), "\nCorrelation")) +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
  return(ggheatmap)
}