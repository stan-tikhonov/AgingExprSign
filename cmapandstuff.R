# this script prints lists of top significant genes for the artyomovlab tool
library(tidyverse)
library(annotate)
library(org.Hs.eg.db)


for (name in names(agingsignatures_v3)){
  if (!(name %in% c("All", "Mouse", "Muscle"))){
    temp = agingsignatures_v3[[name]] %>% rownames_to_column("rownames") %>% filter(adj_pval < 0.05) %>% column_to_rownames("rownames")
    answer = ""
    for (gene in rownames(temp)){
      answer = paste0(answer, gene, " ")
    }
    print(paste0("Signature: ", name))
    print(answer)
  } else {
    temp = agingsignatures_v3[[name]] %>% rownames_to_column("rownames") %>% filter(adj_pval_LOO < 0.05) %>% column_to_rownames("rownames")
    answer = ""
    for (gene in rownames(temp)){
      answer = paste0(answer, gene, " ")
    }
    print(paste0("Signature: ", name))
    print(answer)
  }
}

# Association test
source("FUN_Association_test.r")

lookup = read.table(file = "GSE92742_Broad_LINCS_gene_info.txt", header = T, sep = "\t")
lookup = lookup[,1:2]
lookup$pr_gene_symbol = as.character(lookup$pr_gene_symbol)

kek = rbind(
c(91860, "CALML4"),
c(51275, "MAPKAPK5-AS1"),
c(57175, 'CORO1B'),
c(126298, 'IRGQ'),
c(158747, 'MOSPD2'),
c(64748, 'PLPPR2'),
c(78992, 'YIPF2'),
c(63897, 'HEATR6'),
c(85002, 'FAM86B1'),
c(283232, 'TMEM80'),
c(25960, 'ADGRA2'),
c(89941, 'RHOT2'),
c(79874, 'RABEP2'),
c(100289678, 'ZNF783'),
c(6376, 'CX3CL1'),
c(79716, 'NPEPL1'),
c(11033, 'ADAP1'),
c(4034, 'LRCH4'),
c(399664, 'MEX3D'),
c(54869, 'EPS8L1'),
c(90379, 'DCAF15'),
c(60, 'ACTB'))

colnames(kek) = colnames(lookup)
kek = as.data.frame(kek)

lookup = rbind(lookup, kek)
lookup$pr_gene_id = as.character(lookup$pr_gene_id)

zscoretable = read.table(file = "Inchi_GSE92742.csv", header = T, sep = "\t")
#zscoretable = read.table(file = "trt_sh_GSE92742_n1.csv", header = T, sep = "\t")
zscoretable = zscoretable %>% distinct(sample_id, .keep_all = T) %>% column_to_rownames("sample_id")
zscoretable = t(zscoretable)

zscoretable = as.data.frame(zscoretable)
zscoretable = zscoretable %>% rownames_to_column("pr_gene_id")
zscoretable$pr_gene_id = sub("X", "", zscoretable$pr_gene_id)

zscoretable = left_join(zscoretable, lookup)

zscoretable = na.omit(zscoretable, cols = pr_gene_symbol)
zscoretable = zscoretable %>% remove_rownames() %>% column_to_rownames(var = "pr_gene_symbol")
zscoretable$pr_gene_id = NULL

zscoretable = zscoretable[,grepl(":24 h:", colnames(zscoretable))]

source("FUN.Entrez_converter.R")
agingsignatures_v3_for_gsea = list()
for (name in names(agingsignatures_v3)){
  agingsignatures_v3_for_gsea[[name]] = entrez_converter(agingsignatures_v3[[name]], from = "Mouse", to = "Human")
  symbols = getSYMBOL(rownames(agingsignatures_v3_for_gsea[[name]]), data = "org.Hs.eg.db")
  agingsignatures_v3_for_gsea[[name]]$genesymbol = symbols[rownames(agingsignatures_v3_for_gsea[[name]])]
}

agingsignatures_v3_for_cmap = list()

for (name in names(agingsignatures_v3_for_gsea)){
  if (name %in% c("All", "Mouse")){
    agingsignatures_v3_for_cmap[[name]] = list()
    temp = agingsignatures_v3_for_gsea[[name]] %>% rownames_to_column("row.names") %>% filter(logFC > 0) %>% arrange(adj_pval) %>% column_to_rownames("row.names")
    agingsignatures_v3_for_cmap[[name]][["Up"]] = temp[1:250,]$genesymbol
    temp = agingsignatures_v3_for_gsea[[name]] %>% rownames_to_column("row.names") %>% filter(logFC < 0) %>% arrange(adj_pval) %>% column_to_rownames("row.names")
    agingsignatures_v3_for_cmap[[name]][["Down"]] = temp[1:250,]$genesymbol
  } else {
    agingsignatures_v3_for_cmap[[name]] = list()
    temp = agingsignatures_v3_for_gsea[[name]] %>% rownames_to_column("row.names") %>% filter(logFC > 0 & adj_pval < 0.05) %>% arrange(adj_pval) %>% column_to_rownames("row.names")
    agingsignatures_v3_for_cmap[[name]][["Up"]] = temp$genesymbol
    temp = agingsignatures_v3_for_gsea[[name]] %>% rownames_to_column("row.names") %>% filter(logFC < 0 & adj_pval < 0.05) %>% arrange(adj_pval) %>% column_to_rownames("row.names")
    agingsignatures_v3_for_cmap[[name]][["Down"]] = temp$genesymbol
  }
}

cmapout = Association_test(agingsignatures_v3_for_cmap, zscoretable)

cmapout1 = cmapout

cmapout$pvalue$Signature = NULL

cmapout$adjpval = apply(cmapout$pvalue, 1, p.adjust)

cmapout$adjpval = matrix(p.adjust(as.vector(as.matrix(cmapout$pvalue)), method = "BH"), ncol = length(colnames(cmapout$pvalue)))

colnames(cmapout$adjpval) = colnames(cmapout$pvalue)
rownames(cmapout$adjpval) = rownames(cmapout$pvalue)

cmapout$adjpval = p.adjust(cmapout$pvalue)





