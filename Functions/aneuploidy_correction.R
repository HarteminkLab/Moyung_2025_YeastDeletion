## Function to correct for aneuploidy events in the mutant dataset
## Kevin Moyung
## MacAlpine Lab

correct_aneuploidy <- function(samples, mut) {
  # Mutants with duplicate chromosomes/aneuploidies
  dupe_mutants <- c("CTF4", "MED1", "MOT3", "NUP170", "PAF1", "RSC2", "SRB8", "SWI6", "XRN1")
  dupe_chrs <- c("chrI", "chrXII", "chrV", "chrVIII", "chrXI", "chrVII", "chrXII", "chrV", "chrXI")
  
  dupes.df <- data.frame(Mutant = dupe_mutants, Chr = dupe_chrs)
  
  # Check if the current mutant is in this list; if so, downsample the specific chromosome to the same coverage as chrIV
  if (mut %in% dupe_mutants) {
    cur_chr <- dupes.df[dupes.df$Mutant == mut, "Chr"]
    # Calculate ratio of chrIV reads to current chr, scaled by chr length
    scale_ratio <- (nrow(samples[["chrIV"]]) / coordinates[coordinates$V1 == "chrIV", "V3"]) / 
      (nrow(samples[[cur_chr]]) / coordinates[coordinates$V1 == cur_chr, "V3"])
    samples[[cur_chr]] <- samples[[cur_chr]] %>% sample_frac(scale_ratio)
  } 
  
  # Return list of dataframes with the corrected chrs
  return(samples)
}

