get_midpoint_dataframe<-function(Gene,chr,start_pos,end_pos){
  # mutants <- read.table("/Users/kmoyung/MacAlpine_KM/HARDAC_scripts/Mutants_List.txt", header = F)
  dm_id <- paste0("DM", mutants$V2[mutants$V1 == Gene])
  filename=list.files(path=paste0('/Users/kmoyung/aligned_DM/', dm_id),pattern = glob2rx("*.bam"),full.names = T)
  range.gr=GRanges(seqnames = chr,
                   ranges = IRanges(start=start_pos-500,end=end_pos+500))
  p = ScanBamParam(what = c("rname","pos", "strand", "isize"),which=range.gr)
  reads.l=scanBam(filename,param = p)
  reads.gr=GRanges(seqnames = reads.l[[1]][['rname']],
                   ranges = IRanges(start = reads.l[[1]][['pos']],width = reads.l[[1]][['isize']]))
  reads.gr=reads.gr[width(reads.gr)<=250]
  mnase.df=data.frame(mid=start(reads.gr)+floor((width(reads.gr)-1)/2),length=width(reads.gr))
  mnase_mid.df=subset(mnase.df,mid>=start_pos&mid<=end_pos)
  return(mnase_mid.df)
  }
