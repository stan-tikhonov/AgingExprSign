#################################################################################3
#This script performs p-value and FDR calculation for RNAseq data

adjusted_pvalue_RNAseq_calculator <- function(exprs_data,gene_list,samplename_matrix){
  #make some corrections (if needed)
  samplename_matrix = as.matrix(samplename_matrix)
  data_for_analysis = exprs(exprs_data)
  
  #create output matrices
  phenotype_pairs <- levels(factor(samplename_matrix[,1]))
  output_matrix <- matrix(nrow=length(gene_list),ncol=length(phenotype_pairs))
  rownames(output_matrix) <- gene_list
  colnames(output_matrix) <- phenotype_pairs
  output_matrix_pvalue <- output_matrix
  output_matrix_FDR <- output_matrix

  #Time for linear model; we will build it for now without considering batch effects
  edger_matrix = DGEList(data_for_analysis)
  edger_matrix = calcNormFactors(edger_matrix, method = "RLE")
  
  #Create design matrix
  control = rep(1,ncol(data_for_analysis))
  design <- control
  for (i in seq_along(phenotype_pairs)){
    int_names <- samplename_matrix[samplename_matrix[,1]==phenotype_pairs[i],2]
    temp <- as.numeric(colnames(data_for_analysis) %in% int_names)
    design <- cbind(design,temp)
  }
  rownames(design) <- colnames(data_for_analysis)
  colnames(design) <- c("control",phenotype_pairs)
  
  #Fit
  edger_matrix = estimateDisp(edger_matrix, design,robust = TRUE)
  #plotBCV(edger_matrix)
  fit = glmFit(edger_matrix, design)
  for (i in seq_along(phenotype_pairs)){
    result <- topTags(glmLRT(fit,coef=phenotype_pairs[i]),n=Inf,adjust.method = "BH")$table
    output_matrix_pvalue[rownames(result),phenotype_pairs[i]] <- result$PValue
    output_matrix_FDR[rownames(result),phenotype_pairs[i]] <- result$FDR
  }
  
 
  return(list(output_matrix_pvalue,output_matrix_FDR))
}