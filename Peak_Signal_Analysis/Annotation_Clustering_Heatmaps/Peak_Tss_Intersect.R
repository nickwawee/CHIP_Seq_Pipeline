## This script is a copy of the TSS intersection script that can be ran on the server.
#Its inputs are the bed files that contain peak AUCs, the number of cores, kmeans cluster numbers for the promoter and enhancer, and the intracluster sorting method.
#It outputs sorted bed files in separate kmeans clusters that are used to compute binned matrices and generate heatmaps
#libraries
library(R.methodsS3, lib.loc = "/secondary/projects/mnp/nick/software/R/lib/")
library(R.oo, lib.loc = "/secondary/projects/mnp/nick/software/R/lib/")
library(R.utils, lib.loc = "/secondary/projects/mnp/nick/software/R/lib/")
library(data.table, lib.loc = "/secondary/projects/mnp/nick/software/R/lib/")
library(parallel)

#functions
bednames<-function(bedmat){#this function takes a BED3 file and returns a single character vector of names for regions
  bedmat[,1]=paste(bedmat[,1],rep(":", nrow(bedmat)),sep="")
  bedmat[,2]=paste(as.character(bedmat[,2]), as.character(bedmat[,3]),sep="-")
  bednames<-paste(bedmat[,1],bedmat[,2],sep="")      
}

#parameters
ncores=20
prom_cluster_num=9
enhancer_cluster_num=9
intracluster="mean"
filterd9='yes'
hm="H3K27ac"
datatype="seqdepthnorm"
##loading and combining
sp=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",datatype,"/bedfiles/", sep="")

allfiles=list.files(path = sp, pattern="seq_depth_norm")
if (length(grep(pattern="sortedregions", x=allfiles))>0){
     allfiles=allfiles[-grep(pattern="sortedregions", x=allfiles)]
}

tb=read.csv(paste(sp,allfiles[1],sep=""),sep="\t",header=FALSE, stringsAsFactors = FALSE)
combined<-list()
combined[["chr"]]<-tb$V1
combined[["start"]]<-tb$V2
combined[["end"]]<-tb$V3
for (filen in allfiles){
  tb=read.csv(paste(sp,filen,sep=""),sep="\t",header=FALSE, stringsAsFactors = FALSE)
  n=paste(strsplit(filen, split="_")[[1]][1:3], collapse ="_")
  combined[[n]]=tb$V4
}
save(combined,file=paste('/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/',hm,"_",datatype,'/matrixfiles/combinedlist.RData', sep=""))

enrichmentmat<-do.call(cbind,combined[4:length(combined)])
rownames(enrichmentmat)=paste(combined$chr,":",combined$start,"-",combined$end,sep="")
if (filterd9=='yes'){
   enrichmentmat<-enrichmentmat[,-grep(x=colnames(enrichmentmat), pattern='nnatd9')]
}
#loading TSS to identify promoters
load("/primary/projects/mnp/tools/data/repository/organisms/GRCm38.p6_GENCODE/annotation/gene.filter/GRCm36.p6.GENCODE.RefSeq.filt.RData")

Type=unlist(mclapply(1:length(combined$chr),function(i){
  ch=combined$chr[i]; st=combined$start[i]; en=combined$end[i]
  id=which(genes$chr==ch & genes$tss>=st & genes$tss<=en)
  ret=NULL
  if (length(id)==0){ret="Enhancer"}
  if (length(id)!=0){ret="Promoter"}
  return(ret)
},mc.cores=ncores))

enrichmentdf<-as.data.frame(cbind(as.data.frame(enrichmentmat), Type))
save(enrichmentdf, file=paste("/secondary/projects/mnp/nick/Leslie.Yang/adipocyte_data/ChIPseq2/",hm,"_",datatype,"/matrixfiles/enrichmentdf.RData", sep=""))
##clustering
promotermat<-enrichmentdf[which(enrichmentdf$Type=="Promoter"),-ncol(enrichmentdf)]
set.seed(123)
kmeans_res<-kmeans(scale(promotermat), centers=prom_cluster_num, iter.max = 1000)
promoterdf<-as.data.frame(promotermat)
promoterdf$clusternumber<-kmeans_res$cluster

if (intracluster=="ward.D2"){
  #below clusters the kmeans clustering
  ordervec.all<-data.frame(names=as.character(),order=as.integer())
  for (i in 1:prom_cluster_num){
    newdf<-promoterdf[which(promoterdf$clusternumber==i),-ncol(promoterdf)]
    newregions<-rownames(newdf)
    distmat<-dist(newdf)
    ordervec<-hclust(distmat, method="ward.D2")$order
    ordervec.all=rbind(ordervec.all, cbind(newregions,ordervec))
  }
  row.names(ordervec.all)<-ordervec.all$newregions
  promoterdf$hierclust<-as.numeric(as.character(ordervec.all$ordervec))
  promoterdf<-promoterdf[order(promoterdf$clusternumber, -promoterdf$hierclust),]
  
  #enhancer clustering
  enhancerdf<-enrichmentdf[which(enrichmentdf$Type=="Enhancer"),-ncol(enrichmentdf)]
  set.seed(123)
  kmeans_res<-kmeans(scale(enhancerdf), centers=enhancer_cluster_num)
  enhancerdf$clusternumber<-kmeans_res$cluster
  enhancerdf<-enhancerdf[order(enhancerdf$clusternumber),]
  
  #enhancer cluster clustering
  ordervec.all<-data.frame(names=as.character(),order=as.integer())
  for (i in 1:enhancer_cluster_num){
    newdf<-enhancerdf[which(enhancerdf$clusternumber==i),-ncol(enhancerdf)]
    newregions<-rownames(newdf)
    distmat<-dist(newdf)
    ordervec<-hclust(distmat, method="ward.D2")$order
    ordervec.all=rbind(ordervec.all, cbind(newregions,ordervec))
  }
  row.names(ordervec.all)<-ordervec.all$newregions
  enhancerdf$hierclust<-as.numeric(as.character(ordervec.all$ordervec))
  enhancerdf<-enhancerdf[order(enhancerdf$clusternumber, -enhancerdf$hierclust),]
}else{#the clusters will be sorted by the mean AUC value
  ordervec.all<-data.frame(names=as.character(),order=as.integer())
  for (i in 1:prom_cluster_num){
    newdf<-promoterdf[which(promoterdf$clusternumber==i),-ncol(promoterdf)]
    newregions<-rownames(newdf)
    ordervec<-rowMeans(newdf)
    ordervec.all=rbind(ordervec.all, cbind(newregions,ordervec))
  }
  row.names(ordervec.all)<-ordervec.all$newregions
  promoterdf$meanclust<-as.numeric(as.character(ordervec.all$ordervec))
  promoterdf<-promoterdf[order(promoterdf$clusternumber, -promoterdf$meanclust),]
  
  #enhancer clustering
  enhancerdf<-enrichmentdf[which(enrichmentdf$Type=="Enhancer"),-ncol(enrichmentdf)]
  set.seed(123)
  kmeans_res<-kmeans(scale(enhancerdf), centers=enhancer_cluster_num)
  enhancerdf$clusternumber<-kmeans_res$cluster
  enhancerdf<-enhancerdf[order(enhancerdf$clusternumber),]
  
  #enhancer cluster clustering
  ordervec.all<-data.frame(names=as.character(),order=as.integer())
  for (i in 1:enhancer_cluster_num){
    newdf<-enhancerdf[which(enhancerdf$clusternumber==i),-ncol(enhancerdf)]
    newregions<-rownames(newdf)
    ordervec<-rowMeans(newdf)
    ordervec.all=rbind(ordervec.all, cbind(newregions,ordervec))
  }
  row.names(ordervec.all)<-ordervec.all$newregions
  enhancerdf$meanclust<-as.numeric(as.character(ordervec.all$ordervec))
  enhancerdf<-enhancerdf[order(enhancerdf$clusternumber, -enhancerdf$meanclust),]
}

for (h in 1:prom_cluster_num){
  #subsetting per cluster
  promdf<-promoterdf[which(promoterdf$clusternumber==h),]
  #creating promoter bed
  prom_chr<-unlist(mclapply(1:nrow(promdf), function(i){
    strsplit(rownames(promdf),split=":")[[i]][1]
  }, mc.cores=ncores))
  prom_start<-as.integer(unlist(mclapply(1:nrow(promdf), function(i){
    strsplit(strsplit(rownames(promdf),split=":")[[i]][2],split="-")[[1]][1]
  }, mc.cores = ncores)))
  prom_end<-as.integer(unlist(mclapply(1:nrow(promdf), function(i){
    strsplit(rownames(promdf),split="-")[[i]][2]
  }, mc.cores=ncores)))
  prom_bed<-data.frame(chr=as.character(prom_chr), start=as.integer(prom_start), end=as.integer(prom_end))
  
  tssindices2=unlist(mclapply(1:nrow(prom_bed), function(i){
    chr=prom_bed$chr[i]; peakstart=prom_bed$start[i]; peakend=prom_bed$end[i];
    indices=which(chr==genes$chr & peakstart<=genes$tss & peakend>=genes$tss)
  }, mc.cores = ncores))
  prom_bed_tss=cbind(genes$chr[tssindices2], genes$tss[tssindices2], genes$tss[tssindices2]+1)
  #prom_bed_tss=cbind(prom_bed_tss,bednames(prom_bed_tss), genes$strand[tssindices2], rep(1, nrow(prom_bed_tss)))
  write.table(prom_bed_tss, file=paste(sp,"promoter_cluster",as.character(h),".bed",sep=""), quote = F, row.names = F, col.names = F, sep="\t")
}  

#subsetting gene matrix for ps and pe for enhancer heatmap
for (h in 1:enhancer_cluster_num){
  endf<-enhancerdf[which(enhancerdf$clusternumber==h),]
  en_chr<-unlist(mclapply(1:nrow(endf), function(i){
    strsplit(rownames(endf),split=":")[[i]][1]
  }, mc.cores=ncores))
  en_start<-as.integer(unlist(mclapply(1:nrow(endf), function(i){
    strsplit(strsplit(rownames(endf),split=":")[[i]][2],split="-")[[1]][1]
  }, mc.cores = ncores)))
  en_end<-as.integer(unlist(mclapply(1:nrow(endf), function(i){
    strsplit(rownames(endf),split="-")[[i]][2]
  }, mc.cores=ncores)))
  en_bed<-data.frame(chr=as.character(en_chr), start=as.integer(en_start), end=as.integer(en_end))
  write.table(en_bed, file=paste(sp,"enhancer_cluster",as.character(h),".bed",sep=""), quote = F, row.names = F, col.names = F, sep="\t")
}  
