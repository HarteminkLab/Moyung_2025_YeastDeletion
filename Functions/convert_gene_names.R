## Functions to convert between gene names and IDs
## Kevin Moyung
## MacAlpine Lab

id_to_gene <- function(ID) {
  if (ID %in% gene.df$name) {
    return(gene.df[gene.df$name == ID, "alias"])
  } else {
    return(ID)
  }
 
}

gene_to_id <- function(gene) {
  return(gene.df[gene.df$alias %in% gene, "name"])
}
