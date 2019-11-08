#Download_raw_reads is a function for downloading RNA-seq reads from featureCounts output
#phenodata rownames should be the same as names of featureCounts count files

Download_raw_reads <- function(featurecounts_dir,phenodata){
  setwd(featurecounts_dir)
  
  temp_data <- read.csv(dir()[!grepl(".summary$",dir())][1],header=T,sep='\t',skip = 1)
  counts_star <- data.frame(ID=temp_data$Geneid)
  rownames(counts_star) <- counts_star$ID
  for (i in dir()[!grepl(".summary$",dir())]){
    temp_name <- strsplit(i,".count")[[1]][1]
    temp_name <- gsub("-","_",temp_name)
    temp_data <- read.csv(i,header=T,sep='\t',skip = 1)
    temp_counts <- temp_data[,7]
    counts_star[temp_name] <- temp_counts
  }
  rm(temp_data)
  counts_star$ID <- rownames(counts_star)
  counts_star <- counts_star[,-1]
  inter_samples <- intersect(rownames(phenodata),colnames(counts_star))
  counts_star <- counts_star[,inter_samples]

  #Create expression set object
  library("lumi")
  phenodata <- phenodata[inter_samples,]
  meta.info <- data.frame(labelDescription=colnames(phenodata))
  pheno <- new("AnnotatedDataFrame", data = phenodata, varMetadata = meta.info)
  RNAseq_counts_star <- new("ExpressionSet", exprs=as.matrix(counts_star),phenoData=pheno)
  
  return(RNAseq_counts_star)
}