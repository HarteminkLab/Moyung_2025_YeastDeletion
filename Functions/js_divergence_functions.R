## Js Divergence Analysis Functions
## Kevin Moyung
## MacAlpine Lab

library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(matrixStats)
library(data.table)


## Function: get_scores
## Reads in the scores from a specified chromosome, filters out transposable element regions
## and mutants on the same chromosome
get_scores <- function(chr) {
  # Detect OS and read in files from the correct paths
  if (Sys.info()[['sysname']] == "Windows") {
    # List the JS scores
    files <- list.files(paste0("D:/js_scores/", chr,"_scores/"), full.names = TRUE)
    
    # Load the chromosome coordinates
    coordinates <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = FALSE, sep = "\t")
    
    # Load genome features
    features <- read.csv("C:/Users/Kevin/Desktop/MacAlpine_KM/feature_files/sacCer3_genes_for_making_schematic.csv", header = TRUE)
    
    # Load mutants
    mutants_meta <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/Deletion_Mutants_Primers.txt", header = TRUE, sep = "\t")
    
    # Source the typhoon plot code
    source("D:/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
  } else {
    # List the JS scores
    files <- list.files(paste0("/Users/kmoyung/js_scores/", chr,"_scores_sliding/"), full.names = TRUE)
    
    # Load the chromosome coordinates
    coordinates <- read.table("/Users/kmoyung/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = FALSE, sep = "\t")
    
    # Load genome features
    features <- read.csv("/Users/kmoyung/MacAlpine_KM/feature_files/sacCer3_genes_for_making_schematic.csv", header = TRUE)
    
    # Load mutants
    mutants_meta <- read.table("/Users/kmoyung/MacAlpine_KM/Deletion_Mutants_Primers.txt", header = TRUE, sep = "\t")
    
    # Source the typhoon plot code
    source("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
    source('/Users/kmoyung/MacAlpine_KM/Code/get_midpoint_dataframe.R')
  }
  
  start <- coordinates[match(chr, coordinates$V1), 2]
  end <- coordinates[match(chr, coordinates$V1), 3]
  window_size <- 500
  locus_window <- end - start
  
  scores <- read.table(files[1], header = TRUE)
  for (i in 2:length(files)) {
    cur_file <- fread(files[i], header = TRUE)
    scores <- cbind(scores, cur_file)
  }
  
  ## Filter out scores from regions w/transposable elements (automatic, removes locations w median > 0.048)
  medians <- apply(scores, 1, median, na.rm = TRUE)
  scores <- cbind(scores, medians)
  scores[!is.na(scores$medians) & scores$medians > 0.045, ] <- NA
  scores <- scores[, !(colnames(scores) %in% "medians")]
  
  # Optional: Filter out mutant at the site of deletion
  chr_number <- as.numeric(as.roman(substring(chr, 4, nchar(chr))))
  same_chr <- as.character(mutants_meta[mutants_meta$Chromosome == chr_number, ]$Gene_name)
  
  steps <- seq(0, (locus_window - window_size), by = 50)
  for (i in 1:length(same_chr)) {
    sample <- same_chr[i]
    sample_start <- mutants_meta$Start[mutants_meta$Gene_name == sample]
    sample_end <- mutants_meta$End[mutants_meta$Gene_name == sample]
    # Shift up by 500bp to account for the window
    sample_start <- sample_start - 500
    # Round coordinates to the nearest 50bp
    sample_start <- sample_start - (sample_start %% 50)
    sample_end <- sample_end - (sample_end %% 50)
    # Match the position of the scores
    pos1 <- match(sample_start, steps)
    pos2 <- match(sample_end, steps)
    
    scores[pos1:pos2, sample] <- NA
  }
  
  return(scores)
}

## Function: get_promoter_coords
## Obtains the list of features for the chromosome based on coordinates of each gene
get_promoter_coords <- function(chr, upstream, downstream) {
  # Detect OS and read in files from the correct paths
  if (Sys.info()[['sysname']] == "Windows") {
    # List the JS scores
    files <- list.files(paste0("D:/js_scores/", chr,"_scores/"), full.names = TRUE)
    
    # Load the chromosome coordinates
    coordinates <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = FALSE, sep = "\t")
    
    # Load genome features
    features <- read.csv("C:/Users/Kevin/Desktop/MacAlpine_KM/feature_files/sacCer3_genes_for_making_schematic.csv", header = TRUE)
    
    # Load mutants
    mutants_meta <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/Deletion_Mutants_Primers.txt", header = TRUE, sep = "\t")
    
    # Source the typhoon plot code
    source("D:/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
  } else {
    # List the JS scores
    files <- list.files(paste0("/Users/kmoyung/js_scores/", chr,"_scores_sliding/"), full.names = TRUE)
    
    # Load the chromosome coordinates
    coordinates <- read.table("/Users/kmoyung/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = FALSE, sep = "\t")
    
    # Load genome features
    features <- read.csv("/Users/kmoyung/MacAlpine_KM/feature_files/sacCer3_genes_for_making_schematic.csv", header = TRUE)
    
    # Load mutants
    mutants_meta <- read.table("/Users/kmoyung/MacAlpine_KM/Deletion_Mutants_Primers.txt", header = TRUE, sep = "\t")
    
    # Source the typhoon plot code
    source("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
    source('/Users/kmoyung/MacAlpine_KM/Code/get_midpoint_dataframe.R')
  }
  
  start <- coordinates[match(chr, coordinates$V1), 2]
  end <- coordinates[match(chr, coordinates$V1), 3]
  
  #### Get JS scores at promoter regions (upstream 500bp + 2kb)
  chrom_features <- subset(features, chrom == chr)
  
  # Create list of promoter coordinates. Start position is 500bp upstream of gene start site and end position is 1kb from the new start.
  promoter_coords <- data.frame(promoter = chrom_features$alias,
                                strand = chrom_features$strand, 
                                tss = chrom_features$tss,
                                start = ifelse(chrom_features$strand == '+', chrom_features$txStart - upstream, chrom_features$txEnd + upstream), 
                                end = ifelse(chrom_features$strand == '+', chrom_features$txStart + downstream, chrom_features$txEnd - downstream))
  
  # Round coordinates to the nearest 50bp
  promoter_coords$start_rounded <- promoter_coords$start - (promoter_coords$start %% 50)
  promoter_coords$end_rounded <- promoter_coords$end - (promoter_coords$end %% 50)
  
  # Throw out any coordinates > 500bp from the chromosome end
  promoter_coords <- promoter_coords[promoter_coords$start_rounded < (end - 2000), ]
  promoter_coords <- promoter_coords[promoter_coords$start_rounded >= (start + 500), ]
  
  return(promoter_coords)
}

get_steps <- function(chr) {
  # Detect OS and read in files from the correct paths
  if (Sys.info()[['sysname']] == "Windows") {
    # List the JS scores
    files <- list.files(paste0("D:/js_scores/", chr,"_scores/"), full.names = TRUE)
    
    # Load the chromosome coordinates
    coordinates <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = FALSE, sep = "\t")
    
    # Load genome features
    features <- read.csv("C:/Users/Kevin/Desktop/MacAlpine_KM/feature_files/sacCer3_genes_for_making_schematic.csv", header = TRUE)
    
    # Load mutants
    mutants_meta <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/Deletion_Mutants_Primers.txt", header = TRUE, sep = "\t")
    
    # Source the typhoon plot code
    source("D:/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
  } else {
    # List the JS scores
    files <- list.files(paste0("/Users/kmoyung/js_scores/", chr,"_scores_sliding/"), full.names = TRUE)
    
    # Load the chromosome coordinates
    coordinates <- read.table("/Users/kmoyung/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = FALSE, sep = "\t")
    
    # Load genome features
    features <- read.csv("/Users/kmoyung/MacAlpine_KM/feature_files/sacCer3_genes_for_making_schematic.csv", header = TRUE)
    
    # Load mutants
    mutants_meta <- read.table("/Users/kmoyung/MacAlpine_KM/Deletion_Mutants_Primers.txt", header = TRUE, sep = "\t")
    
    # Source the typhoon plot code
    source("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
    source('/Users/kmoyung/MacAlpine_KM/Code/get_midpoint_dataframe.R')
  }
  
  start <- coordinates[match(chr, coordinates$V1), 2]
  end <- coordinates[match(chr, coordinates$V1), 3]
  window_size <- 500
  offset <- 250 # Offset to start at the midpoints of the JS windows
  locus_window <- end - start
  steps <- seq(offset, (locus_window - window_size) + offset, by = 50)
  
  return(steps)
}