## Function to calculate nucleosome disorganization
## Kevin Moyung

# Given a list of mutant-gene pairs, calculate the nucleosome disorganization change
calc_nuc_disorg <- function(direct_indirect.df) {
  direct_indirect.df$Chr <- gene.df$chrom[match(direct_indirect.df$Locus, gene.df$alias)]
  direct_indirect.df$Strand <- gene.df$strand[match(direct_indirect.df$Locus, gene.df$alias)]
  direct_indirect.df$TSS <- gene.df$tss[match(direct_indirect.df$Locus, gene.df$alias)]
  direct_indirect.df$Left <- ifelse(direct_indirect.df$Strand == "+", direct_indirect.df$TSS, direct_indirect.df$TSS - 500)
  direct_indirect.df$Right <- ifelse(direct_indirect.df$Strand == "+", direct_indirect.df$TSS + 500, direct_indirect.df$TSS)
  
  # Calculate nucleosome entropies
  direct_indirect.df$Nuc_Entropy <- 0
  for (i in 1:nrow(direct_indirect.df)) {
    nucs.df <- get_fragments(direct_indirect.df$Mutant[i], direct_indirect.df$Chr[i], direct_indirect.df$Left[i], direct_indirect.df$Right[i], "Nuc")
    # Try a small amount of rolling mean smoothing for the fragment counts before calculating nucleosome entropy
    nucs.df$Count <- zoo::rollmean(nucs.df$Count, k = 5, align = "center", fill = nucs.df$Count)
    nucs.df$Control_Count <- zoo::rollmean(nucs.df$Control_Count, k = 5, align = "center", fill = nucs.df$Control_Count)
    
    # Calculate log2 FC in entropy between mutant and control
    direct_indirect.df$Nuc_Entropy[i] <- log2(entropy::entropy(nucs.df$Count) / entropy::entropy(nucs.df$Control_Count)) 
  }
  
  # Normalize the log2 FC in entropy
  # direct_indirect.df$Nuc_Entropy <- scale(direct_indirect.df$Nuc_Entropy)
  return(scale(direct_indirect.df$Nuc_Entropy))
}

