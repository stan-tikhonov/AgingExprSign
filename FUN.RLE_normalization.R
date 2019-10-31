###############################################
#This script performs RLE normalization and returns normalized dataset

RLE_normalization  <- function(original_dataset){
         library(edgeR)
        #First, we need to get library sizes vector
        lib.size = apply(original_dataset, 2, sum)
        
        #Now we should obtain normalization factors for each samples using RLE technique
        edger.rle = lib.size * calcNormFactors(original_dataset, method="RLE")
        
        #Finally we should divide expression values on normalization factors
        RLE_dataset <- original_dataset
        RLE_dataset = sweep(original_dataset, 2, edger.rle, "/")*10^7
        
        return(RLE_dataset)
}


