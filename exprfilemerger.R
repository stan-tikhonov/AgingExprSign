setwd("./GSE53960")
filelist = list.files(pattern = ".*.txt")
datalist = lapply(filelist, function(x)read.table(x, header=T))
for (i in 1:length(datalist)){datalist[[i]] = datalist[[i]] %>% remove_rownames() %>% column_to_rownames(var = "AceVeiwGeneSymbol")}
filteredexprdata = datalist[[1]]
for (i in 2:length(datalist)){
  filteredexprdata = merge(filteredexprdata, datalist[[i]], by=0)
  filteredexprdata = filteredexprdata %>% column_to_rownames(var = "Row.names")
}