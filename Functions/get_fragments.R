## Function to obtain fragment enrichments

library(fst)

get_fragments <- function(mut, chr, start, end, type) {
  # Check what type of fragments to obtain; nucleosome or subnucleosome
  if (type == "Nuc") {
    frag.df <- read.fst(paste0("/Users/kmoyung/global_nuc_enrichment/", mut, "_global_nuc_enrichment.fst"), 
                        from = coordinates$Cum_Start[match(chr, coordinates$Chr)], 
                        to = coordinates$Cum_Length[match(chr, coordinates$Chr)])
  } else {
    frag.df <- read.fst(paste0("/Users/kmoyung/global_subnuc_enrichment/", mut, "_global_subnuc_enrichment.fst"), 
                        from = coordinates$Cum_Start[match(chr, coordinates$Chr)], 
                        to = coordinates$Cum_Length[match(chr, coordinates$Chr)])
  }
  
  # Subset for the location in the chromosome
  frag.df <- subset(frag.df, between(Pos, start, end))
  
  return(frag.df)
}
