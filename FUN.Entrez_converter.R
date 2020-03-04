# convert entrez ids:

library(biomaRt)
tryCatch(
  {
    # default:
    ensembl <- useMart("ensembl")
  },
  error=function(cond) {
    # in case www.ensembl.org is under maintenance:
    ensembl = useEnsembl("ensembl", host = "uswest.ensembl.org")
  }
)

datasets <- listDatasets(ensembl)
rat_dataset = useDataset("rnorvegicus_gene_ensembl", mart=ensembl)
human_dataset = useDataset("hsapiens_gene_ensembl", mart=ensembl)
mouse_dataset = useDataset("mmusculus_gene_ensembl", mart=ensembl)

tryCatch(
  {
    # default:
    rat = useMart("ensembl","rnorvegicus_gene_ensembl")
    human = useMart("ensembl", "hsapiens_gene_ensembl")
    mouse = useMart("ensembl","mmusculus_gene_ensembl")
  },
  error=function(cond) {
    # in case www.ensembl.org is under maintenance:
    rat = useEnsembl("ensembl","rnorvegicus_gene_ensembl", host = "uswest.ensembl.org")
    human = useEnsembl("ensembl", "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
    mouse = useEnsembl("ensembl","mmusculus_gene_ensembl", host = "uswest.ensembl.org")
  }
)

# make entrez maps (human-mouse, and rat-mouse):
Rat_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                 mart=rat_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Rat_to_mouse_orthologs) <- c("Rat_Entrez","Mouse_Entrez")
Mouse_to_rat_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                 values=as.character(unique(Rat_to_mouse_orthologs$Mouse_Entrez)),
                                 mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=rat_dataset)
colnames(Mouse_to_rat_orthologs) <- c("Mouse_Entrez","Rat_Entrez")
Rat_to_mouse_orthologs1 = Rat_to_mouse_orthologs %>% group_by(Rat_Entrez) %>% filter(n() == 1)
Rat_to_mouse_orthologs1 = na.omit(Rat_to_mouse_orthologs1)
Mouse_to_rat_orthologs1 = subset(Mouse_to_rat_orthologs, Mouse_Entrez %in% Rat_to_mouse_orthologs1$Mouse_Entrez)
Mouse_to_rat_orthologs1 = Mouse_to_rat_orthologs1 %>% group_by(Mouse_Entrez) %>% filter(n() == 1)
mouse_rat_entrez_map = na.omit(Mouse_to_rat_orthologs1)
mouse_rat_entrez_map = as.data.frame(mouse_rat_entrez_map)
mouse_rat_entrez_map = mouse_rat_entrez_map %>% mutate_all(as.character)

Human_to_mouse_orthologs <- getLDS(attributes=c("entrezgene_id"),
                                   mart=human_dataset,attributesL=c("entrezgene_id"), martL=mouse_dataset)
colnames(Human_to_mouse_orthologs) <- c("Human_Entrez","Mouse_Entrez")
Mouse_to_human_orthologs <- getLDS(attributes=c("entrezgene_id"), filters="entrezgene_id",
                                   values=as.character(unique(Human_to_mouse_orthologs$Mouse_Entrez)),
                                   mart=mouse_dataset,attributesL=c("entrezgene_id"), martL=human_dataset)
colnames(Mouse_to_human_orthologs) <- c("Mouse_Entrez","Human_Entrez")
Human_to_mouse_orthologs1 = Human_to_mouse_orthologs %>% group_by(Human_Entrez) %>% filter(n() == 1)
Human_to_mouse_orthologs1 = na.omit(Human_to_mouse_orthologs1)
Mouse_to_human_orthologs1 = subset(Mouse_to_human_orthologs, Mouse_Entrez %in% Human_to_mouse_orthologs1$Mouse_Entrez)
Mouse_to_human_orthologs1 = Mouse_to_human_orthologs1 %>% group_by(Mouse_Entrez) %>% filter(n() == 1)
mouse_human_entrez_map = na.omit(Mouse_to_human_orthologs1)
mouse_human_entrez_map = as.data.frame(mouse_human_entrez_map)
mouse_human_entrez_map = mouse_human_entrez_map %>% mutate_all(as.character)

entrez_converter = function(dataset, from, to){
  # convert entrez:
  if (from == "Human" & to == "Mouse"){
    exd <- dataset
    rownames(exd) = as.character(rownames(exd))
    exd$Human_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_human_entrez_map)
    exd = na.omit(exd, cols=Mouse_Entrez)
    dataset <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
  }
  if (from == "Rat" & to == "Mouse"){
    exd <- dataset
    rownames(exd) = as.character(rownames(exd))
    exd$Rat_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_rat_entrez_map)
    exd = na.omit(exd, cols=Mouse_Entrez)
    dataset <- exd %>% remove_rownames() %>% column_to_rownames(var = "Mouse_Entrez")
  }
  if (from == "Mouse" & to == "Human"){
    exd <- dataset
    rownames(exd) = as.character(rownames(exd))
    exd$Mouse_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_human_entrez_map)
    exd = na.omit(exd, cols=Human_Entrez)
    dataset <- exd %>% remove_rownames() %>% column_to_rownames(var = "Human_Entrez")
  }
  if (from == "Mouse" & to == "Rat"){
    exd <- dataset
    rownames(exd) = as.character(rownames(exd))
    exd$Mouse_Entrez = rownames(exd)
    exd <- left_join(exd, mouse_rat_entrez_map)
    exd = na.omit(exd, cols=Rat_Entrez)
    dataset <- exd %>% remove_rownames() %>% column_to_rownames(var = "Rat_Entrez")
  }
  
  if (from == "Rat" & to == "Human"){
    print("Stasyan was a lazy ass and didn't write the code for this combination of 'from' and 'to'")
  }
  if (from == "Human" & to == "Rat"){
    print("Stasyan was a lazy ass and didn't write the code for this combination of 'from' and 'to'")
  }
  
  return(dataset)
}








