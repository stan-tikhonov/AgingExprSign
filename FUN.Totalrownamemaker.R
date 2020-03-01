totalrownamemaker = function(logFCunlisted, thres){
  # make totalrownamematrix:
  totalrownamematrix = list()
  for (i in 1:(length(logFCunlisted)-1)){
    for (j in (i + 1):length(logFCunlisted)){
      topA = logFCunlisted[[i]] %>% rownames_to_column(var = "row.names")
      topA = topA %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topA = topA %>% column_to_rownames(var = "row.names")
      topB = logFCunlisted[[j]] %>% rownames_to_column(var = "row.names")
      topB = topB %>% top_n(-1 * as.integer(as.character(thres)), adj.P.Val)
      topB = topB %>% column_to_rownames(var = "row.names")
      totalrownames = union(rownames(topA), rownames(topB))
      #tempdata = matrix(nrow = length(totalrownames), ncol = 2)
      #rownames(tempdata) = totalrownames
      #tempdata[, 1] = logFCunlisted[[i]][totalrownames,]$logFC
      #tempdata[, 2] = logFCunlisted[[j]][totalrownames,]$logFC
      #tempdata = na.omit(tempdata)
      #totalrownames = rownames(tempdata)
      rownamesA = rownames(subset(logFCunlisted[[i]], rownames(logFCunlisted[[i]]) %in% totalrownames))
      rownamesB = rownames(subset(logFCunlisted[[j]], rownames(logFCunlisted[[j]]) %in% totalrownames))
      totalrownamematrix[[names(logFCunlisted)[i]]][[names(logFCunlisted)[j]]] = intersect(rownamesA, rownamesB)
      totalrownamematrix[[names(logFCunlisted)[j]]][[names(logFCunlisted)[i]]] = intersect(rownamesA, rownamesB)
    }
  }
  return(totalrownamematrix)
}