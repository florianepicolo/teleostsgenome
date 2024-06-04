##### LIBRARY #####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", "ape")
BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(ape)

##### FUNCTIONS #####
recup_species_list <- function(){
  tree <- read.tree("tree_v112.txt")
  node_id <- which(tree$tip.label == "lepisosteus_oculatus") 
  parent_node <- tree$edge[tree$edge[,2] == node_id, 1] +1
  subtree <- extract.clade(tree, node=parent_node)
  species_list <- subtree$tip.label
  
  salmo_id = which(tree$tip.label == "esox_lucius")
  parent_node <- tree$edge[tree$edge[,2] == salmo_id, 1]
  subtree_salmo <- extract.clade(tree, node=parent_node)
  salmo_species = subtree_salmo$tip.label
  carp_id = which(tree$tip.label == "sinocyclocheilus_grahami")
  parent_node <- tree$edge[tree$edge[,2] == carp_id, 1]
  subtree_carp <- extract.clade(tree, node=parent_node)
  carp_species = subtree_carp$tip.label

  return(list(subtree = subtree, species_list = list(teleosts = species_list, salmo = salmo_species, carp = carp_species)))
}
transform_species_name <- function(species_name) {
  parts <- unlist(strsplit(species_name, "_"))
  if (length(parts) == 2) {
    first_letter <- substring(parts[1], 1, 1)
    rest <- parts[2]
    result <- paste0(first_letter, rest, "_gene_ensembl")
  } else {
    first_letters <- sapply(parts[-length(parts)], function(x) substring(x, 1, 1))
    rest <- parts[length(parts)]
    result <- paste0(paste(first_letters, collapse = ""), rest, "_gene_ensembl")
  }
  return(result)
}
plot_tree_with_color <- function(tree, species_ids, color, species_to_remove) {
  tree <- drop.tip(tree, species_to_remove)
  plot(tree, cex = 0.8, no.margin = TRUE, adj = 0.15, main = "sous-arbre des téléostéens")
  tiplabels(pch = 19, col = ifelse(tree$tip.label %in% species_ids, color, "black"))
}
get_genes_for_species <- function(specie) {
  ensembl <- useMart("ensembl", dataset = specie)
  
  ortho <- getBM(attributes = c('ensembl_gene_id', "loculatus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_subtype"),
                 mart = ensembl)
  para <- getBM(attributes = c('ensembl_gene_id', paste0(substr(specie, 1, nchar(specie) - 13), "_paralog_ensembl_gene")),
                mart = ensembl)
  colnames(ortho)[1] <- specie
  colnames(para)[1] <- specie
  return(list(ortho = ortho, para = para))
}
count_paralogues <- function(para){
  paralog_counts <- para %>%
    group_by(sformosus_gene_ensembl) %>%
    summarise(number_of_paralogues = n_distinct(sformosus_paralog_ensembl_gene))
  paralog_distribution <- paralog_counts %>%
    group_by(number_of_paralogues) %>%
    summarise(count = n()) %>%
    ungroup()
  total_genes <- sum(paralog_distribution$count)
  paralog_distribution <- paralog_distribution %>%
    mutate(percentage = (count / total_genes) * 100)
  return(paralog_distribution)
}

##### PROGRAMM #####
results <- recup_species_list()
subtree = results[[1]]
species_list = results[[2]]

ensembl = useMart("ensembl") 
dataset = listDatasets(ensembl)

transformed_species_list <- sapply(species_list$teleosts, transform_species_name)
matches <- transformed_species_list %in% dataset$dataset
matched_species <- species_list$teleosts[matches]
unmatched_species <- species_list$teleosts[!matches]

list(matched_species = matched_species, unmatched_species = unmatched_species)
plot_tree_with_color(subtree, c(species_list$salmo, species_list$carp), "red", unmatched_species)

global_cp <- data.frame(number_of_paralogues = integer(),
                        count = integer(),
                        percentage = numeric(),
                        specie = character(),
                        stringsAsFactors = FALSE)


for(sp in transformed_species_list) {
  print(paste0("######", sp, "######"))
  if(sp %in% dataset$dataset){
    results <- get_genes_for_species(sp)
    write.table(results$ortho, paste0(sp, "_ortho.csv"), row.names = FALSE, sep = ";")
    write.table(results$para, paste0(sp, "_para.csv"), row.names = FALSE, sep = ";")
    
    cp = count_paralogues(results$para) 
    cp$specie = substr(sp, 1, nchar(sp) - 13)
    global_cp <- rbind(global_cp, cp)
    
    ## récupérer la liste des gènes unique avec leur nombre d'ortho
    
    ## chez les ortho regarder le nombre de chez homo et spotted gar ! 
    ## faire un énorme tableau avec homme > spotted gar > all species (très certainement ENORME)
    break
  }
}


