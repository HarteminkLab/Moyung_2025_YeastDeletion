## Get midpoint dataframes for sample and control
## Kevin Moyung
## MacAlpine Lab

library(dplyr)

chrs <- c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 
          'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 
          'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI')

if (Sys.info()[['sysname']] == "Windows") {
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/get_genes.R")
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/roman_conversion.R")
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/peakFinder.R")
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/ggdensity.R")
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/get_motifs.R")
 
  coordinates <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = T)
  
  # Load centromere sites for plotting
  centromeres <- read.table("C:/Users/Kevin/Desktop/MacAlpine_KM/feature_files/centromeres_sacCer3.txt", header = T)
  
  # Source the typhoon plot code
  source("D:/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
  
  # Source the js divergence functions
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/js_divergence_functions.R")
  
  # Source aneuploidy correction
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/aneuploidy_correction.R")
  
  # Get total read counts for each mutant (corrected for aneuploidies)
  total_reads.df <- readRDS("C:/Users/Kevin/Desktop/MacAlpine_KM/Metadata/total_read_counts_perchr.RDS")
  mut_reads.df <- aggregate(Count ~ Mutant, total_reads.df, sum)
  
  # Source the gene name conversion functions
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/convert_gene_names.R")
  
  gene.df = read.csv("D:/HARDAC_scripts/sacCer3_genes_for_making_schematic.csv",header = T)
  
  # Source typhoon plot functions
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/plot_typhoon_functions.R")
  
  # Source the fragment enrichment functions
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/get_fragments.R")
  
  # Source the nucleosome disorganization (entropy) function
  source("C:/Users/Kevin/Desktop/MacAlpine_KM/Code/calculate_nuc_disorganization.R")
} else {
  source("/Users/kmoyung/MacAlpine_KM/Code/get_genes.R")
  source("/Users/kmoyung/MacAlpine_KM/Code/roman_conversion.R")
  source("/Users/kmoyung/MacAlpine_KM/Code/peakFinder.R")
  source("/Users/kmoyung/MacAlpine_KM/Code/ggdensity.R")
  source("/Users/kmoyung/MacAlpine_KM/Code/get_motifs.R")
  
  coordinates <- read.table("/Users/kmoyung/MacAlpine_KM/chromosome_features/chromosome_coordinates.txt", header = T)
  
  # Load centromere sites for plotting
  centromeres <- read.table("/Users/kmoyung/MacAlpine_KM/feature_files/centromeres_sacCer3.txt", header = T)
  
  # Source the typhoon plot code
  source("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/surrounding_typhoon_plot_functions_local.R")
  source('/Users/kmoyung/MacAlpine_KM/Code/get_midpoint_dataframe.R')
  
  # Source the js divergence functions
  source("/Users/kmoyung/MacAlpine_KM/Code/js_divergence_functions.R")
  
  # Source aneuploidy correction
  source("/Users/kmoyung/MacAlpine_KM/Code/aneuploidy_correction.R")
  
  # Get total read counts for each mutant (corrected for aneuploidies)
  total_reads.df <- readRDS("/Users/kmoyung/MacAlpine_KM/Metadata/total_read_counts_perchr.RDS")
  mut_reads.df <- aggregate(Count ~ Mutant, total_reads.df, sum)
  
  # Source the gene name conversion functions
  source("/Users/kmoyung/MacAlpine_KM/Code/convert_gene_names.R")
  
  gene.df = read.csv("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/sacCer3_genes_for_making_schematic.csv", header = T)
  
  # Source typhoon plot functions
  source("/Users/kmoyung/MacAlpine_KM/Code/plot_typhoon_functions.R")
  
  # Source the fragment enrichment functions
  source("/Users/kmoyung/MacAlpine_KM/Code/get_fragments.R")
  
  # Source the nucleosome disorganization (entropy) function
  source("/Users/kmoyung/MacAlpine_KM/Code/calculate_nuc_disorganization.R")
}

## get_samples(mut): Returns a list of dataframes of reads, one df for every chromosome
get_samples <- function(mut, set_chr) {
  # Check if we want a specific chr or all the chrs
  if (set_chr != "all") {
    cur_chr <- set_chr
    chr_start <- coordinates[match(cur_chr, coordinates$Chr), "Start"]
    chr_end <- coordinates[match(cur_chr, coordinates$Chr), "End"]
    mid.df <- get_midpoint_dataframe(mut, cur_chr, chr_start, chr_end)
    
    # Truncate locations with more than n duplicate reads in the sample
    depth <- mut_reads.df[mut_reads.df$Mutant == mut, "Count"]
    n_dupes <- ceiling(depth / 4000000)
    mid.dt <- as.data.table(mid.df)
    freq_mid.dt <- mid.dt[, .N, by=names(mid.dt)]
    freq_mid.dt[freq_mid.dt$N > n_dupes, "N"] <- n_dupes
    filtered_mid.df <- data.frame(freq_mid.dt[rep(seq_len(dim(freq_mid.dt)[1]), freq_mid.dt$N), 1:2, drop = FALSE], row.names=NULL)
    
    return(filtered_mid.df)
  } else {
    samples <- list()
    for (cur_chr in chrs) {
      chr_start <- coordinates[match(cur_chr, coordinates$Chr), "Start"]
      chr_end <- coordinates[match(cur_chr, coordinates$Chr), "End"]
      mid.df <- get_midpoint_dataframe(mut, cur_chr, chr_start, chr_end)
      
      # Truncate locations with more than n duplicate reads in the sample
      depth <- mut_reads.df[mut_reads.df$Mutant == mut, "Count"]
      n_dupes <- ceiling(depth / 4000000)
      mid.dt <- as.data.table(mid.df)
      freq_mid.dt <- mid.dt[, .N, by=names(mid.dt)]
      freq_mid.dt[freq_mid.dt$N > n_dupes, "N"] <- n_dupes
      filtered_mid.df <- data.frame(freq_mid.dt[rep(seq_len(dim(freq_mid.dt)[1]), freq_mid.dt$N), 1:2, drop = FALSE], row.names=NULL)
      
      samples[[cur_chr]] <- filtered_mid.df
    }
    return(samples)
  }
  
}

## get_controls(): Returns a list of dataframes of reads for the control + downsampled, one df for every chromosome
# set_chr = either a specific chromosome of "all"
# down_mut = select the mutant to downsample the control to, otherwise downsample to the lowest mutant (HTL1)
get_controls <- function(set_chr, down_mut) {
  if (set_chr != "all") {
    cur_chr <- set_chr
    chr_start <- coordinates[match(cur_chr, coordinates$Chr), "Start"]
    chr_end <- coordinates[match(cur_chr, coordinates$Chr), "End"]
    # Read the control and subsample
    if (Sys.info()[['sysname']] == "Windows") {
      control.df <- readRDS(paste0("D:/aligned_mutants/Controls/", cur_chr, "_control.RDS"))
    } else {
      control.df <- readRDS(paste0("/Users/kmoyung/aligned_DM/Controls/", cur_chr, "_control.RDS"))
    }
    
    # Downsample to a mutant (or lowest depth mutant, aka HTL1)
    set.seed(999)
    if (down_mut != "HTL1") {
      min_chr_depths <- subset(total_reads.df, Mutant == down_mut)
      cur_depth = min_chr_depths[min_chr_depths$Chr == cur_chr, "Count"]
      sub_control <- sample_n(control.df, size = cur_depth, replace = T)
    } else {
      min_chr_depths <- subset(total_reads.df, Mutant == "HTL1")
      cur_depth = min_chr_depths[min_chr_depths$Chr == cur_chr, "Count"]
      sub_control <- sample_n(control.df, size = cur_depth, replace = T)
    }
    
    # Truncate locations with more than n duplicate reads in the control
    depth <- mut_reads.df[mut_reads.df$Mutant == down_mut, "Count"]
    n_dupes <- ceiling(depth / 4000000)
    control.dt <- as.data.table(sub_control)
    freq_control.dt <- control.dt[, .N, by=names(control.dt)]
    freq_control.dt[freq_control.dt$N > n_dupes, "N"] <- n_dupes
    filtered_control.df <- data.frame(freq_control.dt[rep(seq_len(dim(freq_control.dt)[1]), freq_control.dt$N), 1:2, drop = FALSE], row.names=NULL)
    
    return(filtered_control.df)
  } else {
    controls <- list()
    for (cur_chr in chrs) {
      chr_start <- coordinates[match(cur_chr, coordinates$Chr), "Start"]
      chr_end <- coordinates[match(cur_chr, coordinates$Chr), "End"]
      # Read the control and subsample
      if (Sys.info()[['sysname']] == "Windows") {
        control.df <- readRDS(paste0("D:/aligned_mutants/Controls/", cur_chr, "_control.RDS"))
      } else {
        control.df <- readRDS(paste0("/Users/kmoyung/aligned_DM/Controls/", cur_chr, "_control.RDS"))
      }
      
      # Downsample to a mutant (or lowest depth mutant, aka HTL1)
      set.seed(999)
      if (down_mut != "HTL1") {
        min_chr_depths <- subset(total_reads.df, Mutant == down_mut)
        cur_depth = min_chr_depths[min_chr_depths$Chr == cur_chr, "Count"]
        sub_control <- sample_n(control.df, size = cur_depth, replace = T)
      } else {
        min_chr_depths <- subset(total_reads.df, Mutant == "HTL1")
        cur_depth = min_chr_depths[min_chr_depths$Chr == cur_chr, "Count"]
        sub_control <- sample_n(control.df, size = cur_depth, replace = T)
      }
      
      # Truncate locations with more than n duplicate reads in the control
      depth <- mut_reads.df[mut_reads.df$Mutant == down_mut, "Count"]
      n_dupes <- ceiling(depth / 4000000)
      control.dt <- as.data.table(sub_control)
      freq_control.dt <- control.dt[, .N, by=names(control.dt)]
      freq_control.dt[freq_control.dt$N > n_dupes, "N"] <- n_dupes
      filtered_control.df <- data.frame(freq_control.dt[rep(seq_len(dim(freq_control.dt)[1]), freq_control.dt$N), 1:2, drop = FALSE], row.names=NULL)
      
      controls[[cur_chr]] <- filtered_control.df
    }
    return(controls)
  }
  
}




