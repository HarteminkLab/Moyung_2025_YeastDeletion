## Find motif occurrences using FIMO

## Create FIMO background by only selecting promoters of genes
get_promoter_background <- function() {
  # Load genes of interest
  gene.df = read.csv("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/sacCer3_genes_for_making_schematic.csv", header = T)
  
  # Load genome
  genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
  
  # Subset for set of promoters
  prom_window <- 250
  promoter_seqs <- list()
  promoters.df <- gene.df %>%
    mutate(Prom_Start = ifelse(strand == "+", tss - prom_window, tss), 
           Prom_End = ifelse(strand == "+", tss, tss + prom_window)) %>%
    dplyr::select(alias, chrom, Prom_Start, Prom_End) %>%
    filter(Prom_Start > 0 & Prom_End > 0)
  
  for (i in 1:nrow(promoters.df)) {
    cur_name <- promoters.df$alias[i]
    cur_chr <- promoters.df$chrom[i]
    cur_start <- promoters.df$Prom_Start[i]
    cur_end <- promoters.df$Prom_End[i]
    promoter_seqs[[cur_name]] <- genome[[cur_chr]][cur_start:cur_end]
  }
  
  promoter_set <- DNAStringSet(promoter_seqs)
  
  saveRDS(promoter_set, file = "/Users/kmoyung/MacAlpine_KM/Metadata/promoter_seqs.RDS")
  
  return(promoter_set)
}

get_motifs <- function(mut, thresh = 1e-4) {
  suppressPackageStartupMessages(library(magrittr))
  suppressPackageStartupMessages(library(universalmotif))
  library(memes)
  library(MotifDb)
  
  # Get MacIsaac sites for FIMO
  # macisaac_sites <- get_macisaac()
  # 
  # # Get Rossi sites for FIMO
  # rossi_sites <- get_rossi_all()
  # rossi_sites$left_coord <- rossi_sites$peakStart - 5
  # rossi_sites$right_coord <- rossi_sites$peakEnd + 5
  # 
  # # Combine both
  # all_sites <- rbind(macisaac_sites[, c("chrom", "left_coord", "right_coord")], 
  #                    rossi_sites[, c("chrom", "left_coord", "right_coord")])
  # 
  # # Load the genome
  # genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
  # genome_seqs <- list()
  # for (i in 1:nrow(all_sites)) {
  #   cur_chr <- all_sites$chrom[i]
  #   start <- all_sites$left_coord[i]
  #   end <- all_sites$right_coord[i]
  #   genome_seqs[[paste0(cur_chr,"_",start,"_",end)]] <- genome[[cur_chr]][start:end]
  # }
  # 
  # genome_set <- DNAStringSet(genome_seqs)
  # 
  # saveRDS(genome_set, file = "/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/MacIsaac_Rossi_genome_set.RDS")

  # genome_set <- readRDS(file = "/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/MacIsaac_Rossi_genome_set.RDS")
  
  # Get all motif scans from entire genome
  # Load the genome
  # genome <- BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3
  # genome_seqs <- list()
  # for (i in 1:length(chrs)) {
  #   cur_chr <- chrs[i]
  #   genome_seqs[[cur_chr]] <- genome[[cur_chr]]
  # }
  # 
  # genome_set <- DNAStringSet(genome_seqs)
  # saveRDS(genome_set, file = "/Users/kmoyung/MacAlpine_KM/Metadata/FIMO_genome_set.RDS")
  
  genome_set <- readRDS("/Users/kmoyung/MacAlpine_KM/Metadata/FIMO_genome_set.RDS")
  
  # TEST: Create background using promoters of genes (250bp upstream of TSS)
  # promoter_set <- readRDS("/Users/kmoyung/MacAlpine_KM/Metadata/promoter_seqs.RDS")
  
  # Query database for motif of interest
  motif <- MotifDb::MotifDb %>%
    MotifDb::query(mut, andStrings=c("Scerevisiae")) %>%
    # Convert from motifdb format to universalmotif format
    universalmotif::convert_motifs() %>% 
    # The result is a list, to simplify the object, return it as a single universalmotif
    .[[1]]
  
  # Run Fimo and convert into dataframe
  fimo_results <- runFimo(genome_set, motif, thresh = thresh, text = F)
  fimo_results <- data.frame(fimo_results)
  fimo_results$seqnames <- as.character(fimo_results$seqnames)
  # split <- as.data.frame(strsplit(fimo_results$seqnames, split = "_"))
  # fimo_results$start <- as.numeric(split[2, ])
  # fimo_results$end <- as.numeric(split[3, ])
  fimo_results$midpoint <- (fimo_results$start + fimo_results$end) / 2
  # fimo_results$seqnames <- as.character(split[1, ])
  
  return(fimo_results)
}

save_all_motifs <- function() {
  if (Sys.info()[['sysname']] == "Windows") {
    source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/get_midpoints.R")
  } else {
    source("/Users/kmoyung/MacAlpine_KM/Code/get_midpoints.R")
  }
  
  for (i in 1:nrow(mutants)) {
    tryCatch({
      mut <- mutants$V1[i]
      dm_num <- as.numeric(mutants$V2[i])
      fimo_results <- get_motifs(mut, 1e-4)
      
      if (i == 1) {
        all_fimo.df <- fimo_results
      } else {
        all_fimo.df <- rbind(all_fimo.df, fimo_results)
      }
      
    }, error=function(e){})
  }
  saveRDS(all_fimo.df, file = paste0("/Users/kmoyung/fimo_motifs/all_fimo_", thresh, ".RDS"))
}
