#script for assembling subread's featureCounts results into an expression counts table

#first, obtain the annotation table using the tools on globe
#then, make it a phenodata table in R:
mousesrapheno = read.csv("annotation_281127.csv")
mousesrapheno = subset(mousesrapheno, mousesrapheno$LibraryStrategy == "RNA-Seq")
mousesrapheno = mousesrapheno %>% remove_rownames() %>% column_to_rownames(var = "Run")

#now, get the expression counts table:
setwd("./count_result_star")
source("../../FUN.Download_raw_reads.R")
mousesraexpr = Download_raw_reads("./", mousesrapheno)

#after that, analyse the percentage of assigned reads and filter out samples of poor quality:
assigneddata = data.frame()
totaldata = data.frame()
file.names <- dir("./", pattern =".count.summary")
for(i in 1:length(file.names)){
  countstable = read.table(file.names[i], sep = "")
  countstable = countstable[-1, ]
  assigneddata[str_remove(file.names[i], ".count.summary"), 1] = as.integer(as.character(countstable$V2[1]))
  totaldata[str_remove(file.names[i], ".count.summary"), 1] = sum(as.integer(as.character(countstable$V2)))}
colnames(assigneddata) = c("Assigned")
colnames(totaldata) = c("Total")
tableforplot = merge(assigneddata, totaldata, by=0)
tableforplot = tableforplot %>% column_to_rownames(var = "Row.names")
plot(tableforplot$Total, tableforplot$Assigned)



