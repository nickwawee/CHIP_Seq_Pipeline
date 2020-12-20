##Peak combine
#This script will combine and merge all peaks called for a specific histone mark
rm(list=ls())
library(GenomicRanges)
library(parallel)

#parameters
#p="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/DNAmapping2/output/bupMACS2/"#path to bed3 files
p="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/DNAmapping2/output/bupMACS2/"
op="/Volumes/projects_secondary/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/Peak_Counting_Compare/"#output path where merged bed3 is going to output
#op="/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/Peak_Counting_Compare/"
nc=4#number of cores
qvalthresh=0.05
#nc=4
mo=1#minimum overlap, (bp)

hms=c("H3K27me3","H3K27ac","H3K9me3","H3K9ac")

mclapply(1:length(hms), function(j){
   hm=hms[j]
   #loading
   bfs=list.files(path=p, pattern=paste(hm,".filtered.BAMPE_peaks.broadPeak", sep=""))
   #bfs=bfs[grep("lean_wt", bfs)]
   dflist=mclapply(1:length(bfs), function(i){
      bf=read.table(file=paste(p, bfs[i], sep=""), sep="\t")
      return(data.frame(chr=as.character(bf$V1), st=as.integer(bf$V2), en=as.integer(bf$V3),qval=10^(-1*bf$V9)))
}, mc.cores = nc/2)
   combined=do.call(rbind, dflist)
   if (is.numeric(qvalthresh)&qvalthresh!=0){
      combined<-combined[which(combined$qval<=qvalthresh),]#fiitering in qval threshold bed regions
   }  
   gr<-GRanges(seqnames = combined$chr, ranges=IRanges(start=combined$st, end=combined$en, names=combined$name))#creating GRanges object
   gr.new<-reduce(gr, min.gapwidth=mo)#merging regions to one bed (or saf) file 
   mergedbed=data.frame(chr=as.character(gr.new@seqnames), st=as.integer(gr.new@ranges@start), en=as.integer(gr.new@ranges@start+gr.new@ranges@width-1))
   write.table(mergedbed, file=paste(op,hm,"_merged.bed", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
}, mc.cores=nc)
