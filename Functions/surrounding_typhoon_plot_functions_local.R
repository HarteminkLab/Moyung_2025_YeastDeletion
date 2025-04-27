## Typhoon plot functions (from Yulong Li)
## Kevin Moyung

# Name of functions contained in this functions file
# FUNCTION NAME			FUNCTION DEFINITION
# get_typhoon_plot_mat  	function(chr, start_pos, end_pos, bam_file_name)
# read_in_paired_bam 		function(bam_file_name, chr, type = "GR")
# dens_dot_plot 		function(dot.m, z_min = 0, z_max = 100, ...)
# make_gene_schematic 		function(feature_chr, feature_start, feature_end, ...)
# plot_gene 			function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 50){
# set_chromatin_schematic	function(x_start = 0, x_end = 1, y_start = 0, y_end = 1){

# Load the libraries
library(GenomicRanges)
library(Rsamtools)
library(gplots)
library(data.table)

# Ensure that the options(stringsAsFactors = FALSE)
options(stringsAsFactors = FALSE)

# Get mutant sample names
if (Sys.info()[['sysname']] == "Windows") {
  mutants <- read.table("D:/HARDAC_scripts/Mutants_List.txt", header = FALSE)
} else {
  mutants <- read.table("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/Mutants_List.txt", header = FALSE)
}

# Obtain density weights from smoothScatter based on the midpoint dataframe
densDM<-function (x, y = NULL, nbin, bandwidth,transformation = function(x) x^1,factor = 1) {
  xy <- xy.coords(x, y) #gets (x,y)  
  select <- is.finite(xy$x) & is.finite(xy$y)  
  x <- cbind(xy$x, xy$y)[select, ] #get rid of weird values  
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)  
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2  
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE) 
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)  
  dens <- map$fhat[cbind(xbin, ybin)]  
  dens[is.na(dens)] <- 0  
  dens[] <- transformation(dens)  
  
  dens
}

# Get the matrix
get_typhoon_plot_mat = function(chr, start_pos, end_pos, bam_file_name){
  
  # Load the bam_file
  chr.gr = read_in_paired_bam(bam_file_name, chr)
  
  # Subset on reads less than 250 bp
  chr.gr = chr.gr[width(chr.gr) <= 250]
  
  # Set the query GR
  query.gr = GRanges(seqnames = chr,
                     ranges = IRanges(start = start_pos, end = end_pos)
  )
  
  # Subset on the indices that overlap with the query.gr
  idx = subjectHits(findOverlaps(query.gr, chr.gr))
  chr.gr = chr.gr[idx]
  
  # Get a new typhoon_plot_gr that will contain the positions of the "half width" reads
  chr_tp.gr = chr.gr
  
  # Update the start and end coordinates of the "half width" reads
  ranges(chr_tp.gr) = IRanges(start = start(chr.gr) + round(width(chr.gr) / 4),
                              end = end(chr.gr) - round(width(chr.gr) / 4)
  )
  
  # Set up the matrix
  mat.m = matrix(0, nrow = 250, ncol = end_pos - start_pos + 1)
  colnames(mat.m) = start_pos:end_pos
  
  # Create the sub feature
  mat.gr = GRanges(seqnames = chr, ranges = IRanges(start = start_pos:end_pos, width = 1))
  
  # Iterate through each fragment width
  for(i in 20:250){
    
    # Get the reads that have a particular fragment width
    idx = which(width(chr.gr) == i)
    
    if(any(idx)){
      
      # Count the overlaps with the mat.gr
      mat.m[i,] = countOverlaps(mat.gr, chr_tp.gr[idx])
      
    }
    
  }
  
  # Return the mat.m
  return(mat.m)
  
}

# Function to convert a BAM file to a GR file
read_in_paired_bam = function(bam_file_name, chr, type = "GR"){
  
  # Get the yeast chr lengths
  yeast_chr = scanBamHeader(bam_file_name)[[1]]$targets[chr]
  
  # Make a GR file for the chromosome
  chr.gr = GRanges(seqnames = chr, ranges = IRanges(start = 1, end = yeast_chr), strand = "*")
  
  # Specify the scan bam paramaeters
  p = ScanBamParam(what = c("pos", "isize"), which = chr.gr, flag = scanBamFlag(isMinusStrand = FALSE))
  
  # Get the reads that meet these conditions
  reads.l = scanBam(bam_file_name, param = p)
  
  if(type == "GR"){
    
    chr.gr = GRanges(seqnames = chr,
                     ranges = IRanges(start = reads.l[[1]][["pos"]],
                                      end = reads.l[[1]][["pos"]] + reads.l[[1]][["isize"]] - 1
                     )
    )
    
    return(chr.gr)
    
  }else{
    
    # Convert this to a data frame
    chr.df = data.frame(chr = chr,
                        start = reads.l[[1]][["pos"]],
                        end = reads.l[[1]][["pos"]] + reads.l[[1]][["isize"]] - 1,
                        width = reads.l[[1]][["isize"]]
    )
    
    return(chr.df)
    
  }
  
  
}	

### Plot density dot typhoon instead of full-length typhoon
if (Sys.info()[['sysname']] == "Windows") {
  source('D:/HARDAC_scripts/get_midpoint_dataframe_hardac.R')
  source('D:/HARDAC_scripts/density_color_hardac.r')
} else {
  source('/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/get_midpoint_dataframe_hardac.R')
  source('/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/density_color_hardac.r')
}

control_density_dot_typhoon_plot<-function(mid.df, start_pos, end_pos, stream, region){
  if (!is.data.frame(mid.df)) {
    mid.df <- fread(mid.df, header = TRUE)
    # Subset using data table (faster)
    mid.df <- mid.df[mid >= start_pos & mid <= end_pos]
    # mid.df=subset(mid.df,mid>=start_pos&mid<=end_pos)
  }
  
  mycol <- colors()[1] # White
  mycol <- c(mycol,"#FEF0D9","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#990000")
  #mycol <- c("#FEF0D9","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#990000")
  densCols=densColors(x=mid.df$mid,y=mid.df$length,nbin=512,bandwidth=c(10,5),transformation = function(x) x^1,colramp = colorRampPalette(mycol)) 
  # Check if plotting upstream or downstream of cassette
  if (stream == "UP") {
    plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=0.5,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0("1kb Upstream Pileup of ", region),xlab="",ylab='Fragment size (bp)',yaxt='n')
    axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
    #abline(v = (start_pos + 1000), col = "blue", lwd = 3, lty = 2)
  }
  else if (stream == "DOWN") {
    plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=0.5,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0("1kb Downstream Pileup of ", region),xlab="",ylab='Fragment size (bp)',yaxt='n')
    axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
    #abline(v = (end_pos - 1000), col = "blue", lwd = 3, lty = 2)
  }
  else {
    plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=1,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0(region),xlab="",ylab='Fragment size (bp)',yaxt='n')
    axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
  }
}


density_dot_typhoon_plot<-function(DM_ID,chr,start_pos,end_pos, stream, region){
  mid.df=get_midpoint_dataframe(DM_ID,chr,start_pos,end_pos)
  mycol <- c("#FEF0D9","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#990000")
  densCols=densColors(x=mid.df$mid,y=mid.df$length,nbin=512,bandwidth=c(10,5),transformation = function(x) x^1,colramp = colorRampPalette(mycol))
  # Check if plotting upstream or downstream of cassette
  # if (stream == "UP") {
  #   plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=0.5,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0("1kb Upstream of ", region, " (", DM_ID, "-)"),xlab=chr,ylab='Fragment size (bp)',yaxt='n')
  #   axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
  #   #abline(v = (start_pos + 1000), col = "blue", lwd = 3, lty = 2)
  # }
  # else if (stream == "DOWN") {
  #   plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=0.5,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0("1kb Downstream of ", region, " (", DM_ID, "-)"),xlab=chr,ylab='Fragment size (bp)',yaxt='n')
  #   axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
  #   #abline(v = (end_pos - 1000), col = "blue", lwd = 3, lty = 2)
  # } 
  # else {
  #   plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=1,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0(region, " (", DM_ID, "-)"),xlab=chr,ylab='Fragment size (bp)',yaxt='n')
  #   axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
  # }
  
  plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=1,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0(region, " (", DM_ID, "-)"),xlab=chr,ylab='Fragment size (bp)',yaxt='n')
    axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
  
  print(nrow(mid.df))
  return(nrow(mid.df))
}

# Samples a subset from the pileup (control) matrix; sample_size is fraction of the total number of reads
sampled_dot_typhoon_plot<-function(mid.df,chr,start_pos,end_pos, sample_size){
  mid.df <- fread(mid.df, header = TRUE)
  print(nrow(mid.df))
  sampled_rows <- sample(1:nrow(mid.df), size = (sample_size))
  # Subset the mid.df by the sampled rows
  mid.df <- mid.df[sampled_rows, ]
  mycol <- c("#FEF0D9","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#990000")
  densCols=densColors(x=mid.df$mid,y=mid.df$length,nbin=512,bandwidth=c(10,5),transformation = function(x) x^1,colramp = colorRampPalette(mycol))
  
  plot(mid.df,bg=densCols,col=densCols,xaxs='i',yaxs='i',pch=21,cex=0.5,lwd=0,ylim=c(0,250),xlim=c(start_pos,end_pos),main=paste0("GAL1-GAL10 Region (", sample_size, " Sample)"),xlab=chr,ylab='Fragment size (bp)',yaxt='n')
  axis(2,at=c(0,100,200),labels = c(0,100,200),tick = T)
  
  print(nrow(mid.df))
}

### Plot gene bodies
# Set up the gene
make_gene_schematic = function(feature_chr, feature_start, feature_end, 
                               y_low = 0, y_high = 0.5, cex_title = 1, bg_type = "white",
                               proteinCoding = T, geneName = T, omit_genes = NA, x_pos_title = 50
){
  
  # Set up the plot
  plot(0, 0, type = "n", bty = "n", bg = bg_type,
       xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
       ylim = c(0, 1), yaxs = "i", yaxt = "n",
       ann = F
  )
  
  # Load the gene dataframe
  if (Sys.info()[['sysname']] == "Windows") {
    gene.df = read.csv("D:/HARDAC_scripts/sacCer3_genes_for_making_schematic.csv",header = T)
    origin.df=read.csv('D:/HARDAC_scripts/oridb_acs_feature_file_jab-curated-798-sites_sacCer3.csv',header = T)
  } else {
    gene.df = read.csv("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/sacCer3_genes_for_making_schematic.csv",header = T)
    origin.df=read.csv('/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/oridb_acs_feature_file_jab-curated-798-sites_sacCer3.csv',header = T)
  }
  origin.df=subset(origin.df,oridb_status=='Confirmed')
  # Convert to a GenomicRanges object
  gene.gr = GRanges(seqnames = gene.df$chr,
                    ranges = IRanges(start = gene.df$left_coord, end = gene.df$right_coord),
                    strand = gene.df$strand
  )
  names(gene.gr) = gene.df$alias
  
  origin.gr=GRanges(seqnames = origin.df$chr,
                    ranges = IRanges(start=origin.df$pos,width = 1),
                    strand = origin.df$strand)
  names(origin.gr)=origin.df$ars_name
  # Create the feature gr
  feature.gr = GRanges(seqnames = feature_chr,	
                       ranges = IRanges(start = feature_start, end = feature_end)
  )
  
  # Find the overlaps
  gene_overlaps.df = as.data.frame(as.matrix(findOverlaps(feature.gr, gene.gr)))
  
  if(any(nrow(gene_overlaps.df))){
    
    # Enter in the genes
    for(i in 1:nrow(gene_overlaps.df)){
      plot_gene(gene.df[gene_overlaps.df$subjectHits[i],], y_low, y_high, 
                feature_start, feature_end, cex_title, geneName, x_pos_title)
    }
    
  }
  ori_overlaps.df = as.data.frame(as.matrix(findOverlaps(feature.gr, origin.gr)))
  
  if(any(nrow(ori_overlaps.df))){
    
    # Enter in the origins
    for(i in 1:nrow(ori_overlaps.df)){
      plot_origin(origin.df[ori_overlaps.df$subjectHits[i],], y_low, y_high, 
                  feature_start, feature_end, cex_title, x_pos_title)
    }
    
  }
  
}

# Plot gene
plot_gene = function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 50){
  
  # Get y_mid
  ymid = (yhigh + ylow) / 2
  
  # Add in the text
  if(gene.v$strand == "+"){
    
    # Make the rectangle
    rect(gene.v$left_coord, ymid, gene.v$right_coord, yhigh - 0.3, col = rgb(1,0,0,0.5))
    
    if(geneName){
      if(gene.v$left_coord >= xstart){
        text(x = gene.v$left_coord + x_pos_title, y = yhigh - 0.2, adj = c(0, 1),
             labels = gene.v$alias, font = 3, cex = cex_title)
      }else{
        text(x = gene.v$end - x_pos_title, y = yhigh - 0.2, adj = c(1, 1),
             labels = gene.v$alias, font = 3, cex = cex_title)
      }
    }
  }else{
    
    # Make the rectangle
    #rect(gene.v$left_coord, ylow + 0.1, gene.v$right_coord, ymid - 0.3, col = rgb(0,0,1,0.5))
    # Make the rectangle
    rect(gene.v$left_coord, ymid, gene.v$right_coord, yhigh - 0.3, col = rgb(0,0,1,0.5))
    
    if(geneName){
      if(gene.v$right_coord <= xend){
        text(x = gene.v$right_coord - x_pos_title, y = yhigh - 0.2, adj = c(0, 1),
             labels = gene.v$alias, srt = 180, font = 3, cex = cex_title)
      }else{
        text(x = gene.v$left_coord + x_pos_title, y = yhigh - 0.2, adj = c(1, 1),
             labels = gene.v$alias, srt = 180, font = 3, cex = cex_title)
      }
    }
  }
  
}

# Plot origin
plot_origin = function(gene.v,ylow, yhigh, xstart, xend, cex_title, x_pos_title = 50){
  
  # Get y_mid
  ymid = (yhigh + ylow) / 2
  
  # Add in the text
  # Add in the text
  if(gene.v$strand == "+"){
    
    # Make the rectangle
    rect(gene.v$pos-10, ymid + 0.1, gene.v$pos+10, yhigh - 0.1, col = rgb(1,0,0,0.5))
    
    #if(geneName){
    # if(gene.v$pos-10 >= xstart){
    #  text(x = gene.v$pos-10 + x_pos_title, y = yhigh - 0.2, adj = c(0, 1),
    #       labels = gene.v$ars_name, font = 3, cex = cex_title)
    #}else{
    # text(x = gene.v$pos+10 - x_pos_title, y = yhigh - 0.2, adj = c(1, 1),
    #     labels = gene.v$ars_name, font = 3, cex = cex_title)
    #}
    # }
  }else{
    
    # Make the rectangle
    rect(gene.v$pos-10, ylow + 0.1, gene.v$pos+10, ymid - 0.1, col = rgb(0,0,1,0.5))
    
    #if(geneName){
    # if(gene.v$pos+10 <= xend){
    #  text(x = gene.v$pos+10 - x_pos_title, y = ylow + 0.2, adj = c(0, 1),
    #      labels = gene.v$ars_name, srt = 180, font = 3, cex = cex_title)
    #}else{
    # text(x = gene.v$pos-10 + x_pos_title, y = ylow + 0.2, adj = c(1, 1),
    #     labels = gene.v$ars_name, srt = 180, font = 3, cex = cex_title)
    #}
    
    
    
    
    
    # }
  }
  
}

# TODO: Set up the TFs
make_TF_schematic = function(feature_chr, feature_start, feature_end, 
                               y_low = 0, y_high = 0.5, cex_title = 0.7, bg_type = "white",
                               proteinCoding = T, geneName = T, omit_genes = NA, x_pos_title = 50
){
  
  # Set up the plot
  plot(0, 0, type = "n", bty = "n", bg = bg_type,
       xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
       ylim = c(0, 1), yaxs = "i", yaxt = "n",
       ann = F
  )
  
  # Load the gene dataframe
  if (Sys.info()[['sysname']] == "Windows") {
    tf.df = read.table("D:/HARDAC_scripts/MacIsaac_p005_c1_V64_SGD.gff3",header = T)
  } else {
    tf.df = read.table("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/MacIsaac_p005_c1_V64_SGD.gff3",header = T)
  }
    
  # Rename columns for consistency
  names(tf.df)[names(tf.df) == "start"] <- "left_coord"
  names(tf.df)[names(tf.df) == "end"] <- "right_coord"
  tf_names <- data.frame(do.call('rbind', strsplit(as.character(tf.df$Name.TF),'=',fixed=TRUE)))
  tf.df <- cbind(tf.df, tf_names$X2)
  names(tf.df)[names(tf.df) == "tf_names$X2"] <- "alias"
  
  # Convert to a GenomicRanges object
  tf.gr = GRanges(seqnames = tf.df$chr,
                    ranges = IRanges(start = tf.df$left_coord, end = tf.df$right_coord),
                    strand = tf.df$strand
  )
  names(tf.gr) = tf_names$alias
  

  # Create the feature gr
  feature.gr = GRanges(seqnames = feature_chr,	
                       ranges = IRanges(start = feature_start, end = feature_end)
  )
  
  # Find the overlaps
  tf_overlaps.df = as.data.frame(as.matrix(findOverlaps(feature.gr, tf.gr)))
  
  
  if(any(nrow(tf_overlaps.df))){
  
    # Enter in the genes
    for(i in 1:nrow(tf_overlaps.df)){
      plot_TF(tf.df[tf_overlaps.df$subjectHits[i],], y_low, y_high, 
                feature_start, feature_end, cex_title, geneName, x_pos_title)
    }
    
  }
}

# Plot TF
plot_TF = function(gene.v, ylow, yhigh, xstart, xend, cex_title, geneName, x_pos_title = 0){
  
  # Get y_mid
  ymid = (yhigh + ylow) / 2
  
  # Add in the text
  if(gene.v$strand == "+"){
    
    # Make the rectangle
    rect(gene.v$left_coord, ymid, gene.v$right_coord, yhigh - 0.3, col = rgb(0,0.75,0,0.5))
    
    if(geneName){
      if(gene.v$left_coord >= xstart){
        text(x = gene.v$left_coord + x_pos_title, y = yhigh - 0.2, adj = c(0, 1),
             labels = gene.v$alias, font = 3, cex = cex_title, srt = 90)
      }else{
        text(x = gene.v$end - x_pos_title, y = yhigh - 0.2, adj = c(1, 1),
             labels = gene.v$alias, font = 3, cex = cex_title, srt = 90)
      }
    }
  }else{
    
    # Make the rectangle
    #rect(gene.v$left_coord, ylow + 0.1, gene.v$right_coord, ymid - 0.3, col = rgb(0,0,1,0.5))
    # Make the rectangle
    rect(gene.v$left_coord, ymid, gene.v$right_coord, yhigh - 0.3, col = rgb(0,0.75,0,0.5))
    
    if(geneName){
      if(gene.v$right_coord <= xend){
        text(x = gene.v$right_coord - x_pos_title, y = yhigh - 0.2, adj = c(0, 1),
             labels = gene.v$alias, font = 3, cex = cex_title, srt = 90)
      }else{
        text(x = gene.v$left_coord + x_pos_title, y = yhigh - 0.2, adj = c(1, 1),
             labels = gene.v$alias, font = 3, cex = cex_title, srt = 90)
      }
    }
  }
  
}