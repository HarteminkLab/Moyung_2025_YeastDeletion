## Functions to plot genes and typhoon plots
## Kevin Moyung
## MacAlpine Lab

# Plot genes
plot_genes <- function(chr, start, end) {
  genes_sub <- get_genes(chr, start - 100, end)
  
  # Check if there are genes to be plotted; if not, plot a blank plot
  if (nrow(genes_sub) > 0) {
    # Plot genes
    group.colors <- c('+' = "#FF9999", '-' = magma(6)[6])
    # Add column to denote positive/negative strand numerically
    genes_sub$strandval <- as.numeric(paste0(genes_sub$strand, 1))
    # Assign tracks to each gene to avoid overlap
    genes_sub$index <- ((seq(1, nrow(genes_sub), by = 1) %% 2) + 1)* genes_sub$strandval
    plot_gene <- ggplot() +
      geom_segment(data = genes_sub,
                   aes(x = tss, xend = pas,
                       y = (((index + 0.5) * 0.4)), yend = (((index + 0.5) * 0.4))),
                   arrow = arrow(length=unit(0.3, "cm"), type = "closed")) +
      geom_segment(data = genes_sub,
                   aes(x = tss, xend = tss,
                       y = (index * 0.4) + 0.02, yend = ((index + 1) * 0.4) - 0.02)) +
      statebins:::geom_rrect(data = genes_sub,
                             aes(xmin = txStart, xmax = txEnd, ymin = (index * 0.4) + 0.02 - 0.2, ymax = ((index + 1) * 0.4) - 0.02 + 0.2, fill = strand), color = "black",
                             alpha = 1) +
      geom_text(data = subset(genes_sub, strand == "+"), aes(x = (txStart + 150), y = ((index + 0.5) * 0.4), label = alias)) +
      geom_text(data = subset(genes_sub, strand == "-"), aes(x = (txEnd - 150), y = ((index + 0.5) * 0.4), label = alias)) +
      scale_fill_manual(values = group.colors) +
      theme_bw() +
      theme(legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      labs(fill = "Strand") +
      ylim(-2, 2) +
      coord_cartesian(xlim = c(start, end), expand = c(0, 0))
  } else {
    plot_gene <- ggplot() +
      theme_bw() +
      theme(legend.position = "right",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      labs(fill = "Strand") +
      ylim(-2, 2) +
      coord_cartesian(xlim = c(start, end), expand = c(0,0))
  }
  
  return(plot_gene)
}

# Plot control typhoon plot
plot_control_typhoon <- function(chr, start, end, mut) {
  # Get the control df
  control.df <- get_controls(chr, mut)
  
  # Subset for the locus
  control.df <- subset(control.df, mid >= start & mid <= end)
  
  # Calculate smoothed density values
  control.df$density <- get_density(control.df$mid, control.df$length,
                                    n = 200,
                                    h = c(25, 30))
  
  # plot_control <- ggplot(control.df[order(control.df$density), ]) +
  #   geom_point(aes(x = mid, y = length), size = 0.55, shape = 20, fill = "gray") +
  #   geom_point(aes(x = mid, y = length, color = density), size = 0.45, shape = 20) +
  #   scale_color_gradientn(colors = rev(magma(6)[-1]), limits = c(0, 0.00003), oob = scales::squish) +
  #   ggtitle(paste0("Control ", chr, ":", start, "-", end)) +
  #   theme_bw() +
  #   theme(panel.grid = element_blank()) +
  #   labs(y = "Length", x = "", color = "Density") +
  #   coord_cartesian(xlim = c(start, end), expand = c(0, 0))
  
  plot_control <- ggplot(control.df[order(control.df$density), ]) +
    # geom_point(aes(x = mid, y = length), size = 0.45, shape = 20, fill = "gray") +
    geom_point(aes(x = mid, y = length, color = density), size = 0.45, shape = 20) +
    scale_color_gradientn(colors = rev(magma(6)[-1]), limits = c(0, 0.00003), oob = scales::squish) +
    ggtitle(paste0("Control")) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "none") +
    labs(y = "Length", x = chr, color = "Density") +
    coord_cartesian(xlim = c(start, end), expand = c(0, 0))
  
  return(plot_control)
}

# Plot sample typhoon plot
plot_sample_typhoon <- function(mut, chr, start, end) {
  
  # Get the control df
  mid.df <- get_samples(mut, chr)
  
  # Subset for the locus
  mid.df <- subset(mid.df, mid >= start & mid <= end)
  
  # Calculate smoothed density values
  mid.df$density <- get_density(mid.df$mid, mid.df$length, 
                                         n = 200, 
                                         h = c(25, 30))
  
  # plot_sample <- ggplot(mid.df[order(mid.df$density), ]) + 
  #   geom_point(aes(x = mid, y = length), size = 0.55, shape = 20, fill = "gray") +
  #   geom_point(aes(x = mid, y = length, color = density), size = 0.45, shape = 20) +
  #   scale_color_gradientn(colors = rev(magma(6)[-1]), limits = c(0, 0.00003), oob = scales::squish) + 
  #   ggtitle(paste0("\u0394", mut, " ", chr, ":", start, "-", end)) +
  #   theme_bw() +
  #   theme(panel.grid = element_blank()) + 
  #   labs(y = "Length", x = "", color = "Density") +
  #   coord_cartesian(xlim = c(start, end), expand = c(0, 0))
  
  plot_sample <- ggplot(mid.df[order(mid.df$density), ]) + 
    # geom_point(aes(x = mid, y = length), size = 0.45, shape = 20, fill = "gray") +
    geom_point(aes(x = mid, y = length, color = density), size = 0.45, shape = 20) +
    scale_color_gradientn(colors = rev(magma(6)[-1]), limits = c(0, 0.00003), oob = scales::squish) + 
    ggtitle(paste0(tolower(mut),"\u0394")) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(face = "italic")) + 
    labs(y = "Length", x = chr, color = "Density") +
    coord_cartesian(xlim = c(start, end), expand = c(0, 0))
  
  return(plot_sample)
}

# Plot TF binding features from Rossi, MacIsaac, and motifs
plot_tf_features <- function(mut){
  # Make a dataframe to store all factors
  factor_sites <- data.frame(position = 0, ypos = 0, Feature = NA)
  
  # Plot motif midpoint at this location
  # List of available mutants with motif information
  if (Sys.info()[['sysname']] == "Windows") {
    fimo_list <- list.files(path = "C:/Users/Kevin/Downloads/fimo_motifs/", full.names = T)
  } else {
    fimo_list <- list.files(path = "/Users/kmoyung/fimo_motifs/", full.names = T)
  }
  if (!isEmpty(grep(mut, fimo_list))) {
    fimo_results <- fread(fimo_list[grep(mut, fimo_list)])
    left <- start
    right <- end
    subset_motifs <- subset(fimo_results, seqnames == chr & midpoint >= left & midpoint <= right)
    if (nrow(subset_motifs) > 0) {
      # plot_TF <- plot_TF + geom_point(data = subset_motifs,
      #                                 aes(x = midpoint, y = 0.25, fill = "Motif"), shape = 18, size = 5, alpha = 0.75)
      factor_sites <- rbind(factor_sites, c(subset_motifs$midpoint[1], 0.75, "Motif"))
    }
  }
  
  # Plot all Rossi binding sites at this location
  rossi_peaks <- get_rossi(mut)
  sub_master_peaks <- subset(rossi_peaks, V1 == chr & V2 >= start & V3 <= end)
  sub_master_peaks$midpoint <- (sub_master_peaks$V2 + sub_master_peaks$V3) / 2
  if (nrow(sub_master_peaks) > 0) {
    # plot_TF <- plot_TF + geom_point(data = subset(sub_master_peaks, toupper(name) == mut),
    #                                  aes(x = (peakStart + peakEnd) / 2, y = 0.75, fill = "Rossi"), shape = 25, size = 5, alpha = 0.75)
    factor_sites <- rbind(factor_sites, c(sub_master_peaks$midpoint[1], 0.5, "Rossi"))
  }

  # Plot MacIsaac Sites
  macisaac_sites <- get_macisaac()
  sub_macisaac_sites <- subset(macisaac_sites, chrom == chr & left_coord >= start & right_coord <= end & alias == mut)
  sub_macisaac_sites$midpoint <- (sub_macisaac_sites$left_coord + sub_macisaac_sites$right_coord) / 2
  if (nrow(sub_macisaac_sites) > 0) {
    # plot_TF <- plot_TF + geom_point(data = subset(sub_macisaac_sites, alias == mut),
    #                                 aes(x = (left_coord + right_coord) / 2, y = 0.5, fill = "MacIsaac"), shape = 22, size = 5, alpha = 0.75)
    factor_sites <- rbind(factor_sites, c(sub_macisaac_sites$midpoint[1], 0.25, "MacIsaac"))
  }

  factor_sites <- factor_sites[-1, ]
  factor_sites$position <- as.numeric(factor_sites$position)
  factor_sites$ypos <- as.numeric(factor_sites$ypos)

  plot_TF <- ggplot(data = factor_sites, aes(x = position, y = ypos)) +
    geom_point(aes(shape = Feature, fill = Feature), size = 5, color = "black") +
    coord_cartesian(xlim = c(start, end), ylim = c(0, 1)) +
    theme_bw() +
    # theme(legend.position = c(0.5, 0.8),
    #       legend.direction = "vertical",
    #       legend.margin = margin(c(0, 0, 0, 0)),
    #       panel.background = element_blank(),
    #       panel.grid = element_blank(),
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.border = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.y = element_blank(),
    #       axis.ticks.y = element_blank(),
    #       axis.line.y = element_blank(),
    #       axis.title.x = element_blank(),
    #       axis.line.x = element_blank(),
    #       axis.text.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    #       plot.margin = unit(c(-5, 0, -5, 0), "cm")) +
    theme(legend.position = "none",
        legend.direction = "vertical",
        legend.margin = margin(c(0, 0, 0, 0)),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(-5, 0, -5, 0), "cm")) +
    scale_shape_manual(values = c(Rossi = 25, MacIsaac = 22, Motif = 18)) +
    scale_fill_manual(values = c(Rossi = "green", MacIsaac = "blue", Motif = "black"))
}
